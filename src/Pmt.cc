////////////////////////////////////////////////////////////////////////////////
// PMT Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////
#include "Utils.h"
#include "Pmt.h"

#include "IdealChargeFunction.h"

#include "TMath.h"
#include "TF1.h"

PMT::PMT( )
{ };

//------------------------------------------------------------------------------

PMT::PMT( int run, int board, int channel )
  : m_run(run)
  , m_board(board)
  , m_channel(channel)
{
  m_pmt_number=m_channel+geo::nchannels*m_board;
  this->initHist();
};

//------------------------------------------------------------------------------

PMT::PMT( int run, int pmt_number)
  : m_run(run)
  , m_pmt_number(pmt_number)
{
  m_board=0; m_channel=0;
  this->getBoardAndChannel(m_pmt_number, m_board, m_channel);
  this->initHist();
};

//------------------------------------------------------------------------------

PMT::~PMT()
{
    //delete h_amplitude; //amplitude of pulse
    //delete h_amplitude_low; // amplitude of pulse in a smaller range
    //delete h_charge; // charge of pulse
    //delete h_pulsetime; //start time of pulse 
    //delete h_pulsepeaktime; //peak time of pulse
    //delete h_NbOfPulse; // number of pulses in a waveform
};

//------------------------------------------------------------------------------

void PMT::initHist()
{
  char hname[100];

  // Amplitude histograms
  sprintf(hname, "hist_run%d_board%d_channel%d_amplitude", m_run, m_board,  m_channel);
  h_amplitude = new TH1D(hname, hname, 100, 0, 2000);
  h_amplitude->SetXTitle("Amplitude (ADC counts)");
  h_amplitude->SetYTitle("Counts");

  sprintf(hname, "hist_run%d_board%d_channel%d_amplitude_low", m_run, m_board,  m_channel);
  h_amplitude_low = new TH1D(hname, hname, 100, 0, 150);
  h_amplitude_low->SetXTitle("Amplitude (ADC counts)");
  h_amplitude_low->SetYTitle("Counts");

  // for charge histograms
  sprintf(hname, "hQ_run%d_board%d_ch%d", m_run , m_board, m_channel);
  h_charge = new TH1D(hname, hname, 200,0,100);
  h_charge->SetXTitle("Charge (pC)");
  h_charge->SetYTitle("Counts");
  // for pulse time
  sprintf(hname, "hPulseTime_run%d_board%d_ch%d", m_run, m_board, m_channel);
  h_pulsetime = new TH1D(hname, hname, 500,400,500);
  h_pulsetime->SetXTitle("Pulse start time (ns)");
  h_pulsetime->SetYTitle("Counts");
  sprintf(hname, "hPulsePeakTime_run%d_board%d_ch%d", m_run, m_board, m_channel);
  h_pulsepeaktime = new TH1D(hname, hname, 500,400,500);
  h_pulsepeaktime->SetXTitle("Pulse start time (ns)");
  h_pulsepeaktime->SetYTitle("Counts");
  // for number of pulse couting
  sprintf(hname, "hPulseCounting_run%d_board%d_ch%d", m_run, m_board, m_channel);
  h_NbOfPulse = new TH1I(hname, hname, 15,0,15);
  h_NbOfPulse->SetXTitle("Number of pulses in the time set window");
  h_NbOfPulse->SetYTitle("Counts");
  
};

//------------------------------------------------------------------------------

void PMT::getBoardAndChannel( int pmt_number, int & board, int & channel )
{
  //Return the board and channel number given the absolute pmt number
  board = pmt_number/geo::nchannels;
  channel = pmt_number%geo::nchannels;
};

//------------------------------------------------------------------------------

bool PMT::isIlluminated( int optical_channel )
{
  // True if a pmt is illuminated by a given optical channel with direct light
  bool m_isIlluminated=false;

  int laser_groups=10;

  // We compare the optical channel with the absolute pmt number.
  // NB: This assumption works only in the case all boards are taken and optical
  // channel number is changed in progressive increasing order.
  if( m_pmt_number>=optical_channel*laser_groups &
                                m_pmt_number<(optical_channel+1)*laser_groups ){
      m_isIlluminated=true;
  }

  return m_isIlluminated;
};

//------------------------------------------------------------------------------

void PMT::loadWaveform( Waveform *waveform )
{

  // Fetch the relevant quantities from the waveform
  double amplitude = abs( waveform->getAmplitude() );
  double charge = abs( waveform->getCharge() );
  double pulsestarttime = waveform->getStartTime();
  double pulsepeaktime = waveform->getPeakTime();
  int npulse = waveform->getNPulses();

  // Fill the histograms
  h_amplitude->Fill(amplitude);
  h_amplitude_low->Fill(amplitude);

  h_charge->Fill(charge);
  h_pulsetime->Fill(pulsestarttime);
  h_pulsepeaktime->Fill(pulsepeaktime);
  h_NbOfPulse->Fill(npulse);

  // Fill the arrays
  // m_amplitude_array.push_back(amplitude);

};

//------------------------------------------------------------------------------
/*
double PMT::ChargeIdealResponse() (double* x, double* par){
  double mu = par[0];
  double q = par[1];
  double sigma = par[2];
  double amplitude = par[3];
  double sum=0;
  for(int n=1; n<50; n++){
    sum += (TMath::Power(mu,n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n)*TMath::Exp(-1.0*(x[0]-q*n)*(x[0]-q*n)/(2.0*n*sigma*sigma))/(sigma*TMath::Sqrt(2.0*TMath::Pi()*n)));
  }
  return amplitude * sum;
}
*/

//------------------------------------------------------------------------------
void PMT::FitChargeWithIdeal(){
  
  IdealChargeFunction idealchargeF_obj;

  TF1* func = new TF1("idealfunc",idealchargeF_obj,0, 100, 4);
  func->SetParNames("mu","q0","sigma","amp.");
  func->SetParameters(20,1.0,0.1,h_charge->Integral()/2.0);
  func->SetParLimits(0, 0, 50);
  func->SetParLimits(1, 0, 5);
  func->SetParLimits(2, 0, 5);
  func->SetParLimits(3, 0, 10.0*h_charge->Integral());
  h_charge->Fit("idealfunc","RQ","",0,80);
  cout<< h_charge->GetName()<<"\t"
      << getHV() << "\t"
      << func->GetParameter(0) <<"\t"
      << func->GetParameter(1) <<"\t"
      << func->GetParameter(2) <<"\t"
      << func->GetParameter(3) <<"\n";
}

//------------------------------------------------------------------------------

void PMT::writeHist(TFile* thefile, TDirectory* ampDir, TDirectory* chargeDir)
{
  
  // Write the histograms to file
  ampDir->cd();
  h_amplitude->Write();
  h_amplitude_low->Write();

  chargeDir->cd();
  h_charge->Write();

};

//------------------------------------------------------------------------------

void PMT::clean()
{
  // m_amplitude_array.clear();
};

//------------------------------------------------------------------------------
