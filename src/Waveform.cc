////////////////////////////////////////////////////////////////////////////////
// Waveform methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"

#include "PulseShapeFunction.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "TCanvas.h"

//------------------------------------------------------------------------------

Waveform::Waveform(){};

//------------------------------------------------------------------------------

Waveform::Waveform(int run, int subrun, int event, int board, int channel)
  : m_run(run)
  , m_subrun(subrun)
  , m_event(event)
  , m_board(board)
  , m_channel(channel)
{ };

//------------------------------------------------------------------------------

Waveform::Waveform( int run, int subrun, int event ,int board, int channel,
                                                      Rawdigits_t raw_waveform )
  : Waveform(run, subrun, event, board, channel)
{
  this->loadData(raw_waveform);
};


//------------------------------------------------------------------------------

Waveform::Waveform(int run, int subrun, int event, int board, int channel, double n_sigma)
  : m_run(run)
  , m_subrun(subrun)
  , m_event(event)
  , m_board(board)
  , m_channel(channel)
  , m_nsigma(n_sigma)
{ };

//------------------------------------------------------------------------------

Waveform::Waveform( int run, int subrun, int event ,int board, int channel, double n_sigma,
                                                      Rawdigits_t raw_waveform )
  : Waveform(run, subrun, event, board, channel, n_sigma)
{
  this->loadData(raw_waveform);
};

//------------------------------------------------------------------------------

Waveform::~Waveform(){};

//------------------------------------------------------------------------------

void Waveform::loadData( Rawdigits_t raw_waveform )
{
  m_nsamples=int(raw_waveform.size());
  m_raw_waveform = raw_waveform;
  for( double w : raw_waveform ){ m_waveform.push_back( double(w) ); }
  BaselineSubtraction();
  //if(hasPulse(m_nsigma)){
    ComputePulseCharacteristics();
  //}
};

//------------------------------------------------------------------------------

void Waveform::BaselineSubtraction()
{
  // Calculate the baseline as the mean values on the first part of the spectrum
  // n_sample_baseline = m_nsamples;

  for(int t=0; t<n_sample_baseline; t++)
  {
    baseline_mean += m_waveform.at(t);
  }
  baseline_mean /= n_sample_baseline;

  // Calculate the stdev of the baseline
  for(int t=0; t<n_sample_baseline; t++)
  {
    baseline_width += pow(m_waveform.at(t)-baseline_mean, 2);
  }
  baseline_width = sqrt( baseline_width / (n_sample_baseline-1) );

  // Subtract the baseline from the signal waveform
  std::transform(m_waveform.begin(), m_waveform.end(), m_waveform.begin(),
                                [ & ] (double x) { return x - baseline_mean; });

};

//------------------------------------------------------------------------------

bool Waveform::find(int run, int subrun, int event, int board, int channel )
{
  return (run==m_run && subrun==m_subrun && event==m_event
                                       && board==m_board && channel==m_channel);
};

//------------------------------------------------------------------------------

bool Waveform::hasSignal(double n_sigma)
{
  // Find if the waveform has n consecutive counts above a threshold expressed
  // in number of sigmas

  bool has_signal=false;
  int n_counts=5.0;
  int counts=0;

  for( double value : m_waveform )
  {
    if( abs(value) > n_sigma*baseline_width ){ counts++; }
    else{ counts=0; } // Reset the counts

    if( counts > n_counts ){ has_signal=true; break; }
  }

  return has_signal;
};

//------------------------------------------------------------------------------

bool Waveform::hasPulse( double n_sigma )
{
  // Define the Pulse region of the signal as the region with >n consecutive
  // counts above a threshold expressed in number of sigmas.

  // this will find the first pulse
  //    (1) with n_counts consecutive samples above threshold
  //        the threshold is defined as n_sigma*baseline_width

  bool has_pulse=false;
  int counts=0;
  int n_counts=5.0;

  //int t_start=0, t_end=0;
  m_start_time=0.0;
  m_end_time=0.0;

  for(int t=0; t<m_nsamples; t++)
  {
    double value = m_waveform.at(t);

    if( abs(value) > n_sigma*baseline_width )
    {
      counts++;
      if( counts==1 ){ m_start_time=t-2; }
    }
    else
    {
      // if it is a pulse region, then this is the ending clause
      if( counts>n_counts && m_start_time>0 )
      {
        m_end_time=t+3;
        has_pulse=true;
        break; // Force to have only one pulse
      }
      else
      {
        // reset everything
        counts=0; m_start_time=0; m_end_time=0;
      }
    }
  }

//  // If no pulse is found, we end the games here..
//  if( !has_pulse ) { return has_pulse; }

  // Here we define the characteristics of the pulse

  //cout<<m_start_time<<endl;

  return has_pulse;
};

//------------------------------------------------------------------------------
void Waveform::ComputePulseCharacteristics(){

  // this will give amplitude, charge, pulse width of the first pulse
  double amp = 0.0, charge=0.0;
  
  //for( int t=m_start_time; t<m_end_time; t++ )
  // here I force the charge to be integrated in a 100 ns window, 
  // starting from a little earlier than the pulse time.
  int startbin = 210, nbins=50;
  for( int t=startbin; t<startbin+nbins; t++ )
  {
    if( m_waveform.at(t) < amp ){ amp = m_waveform.at(t); m_peak_time = t; }
    //charge += abs(m_waveform.at(t));
    charge += m_waveform.at(t); // should not take abs()
  }

  m_width = (m_end_time-m_start_time);
  m_amplitude = -0.122*amp; // convert to mV and inverse to get positive number
  m_integral = -1.0*charge * 0.122*2.0*0.02; // the Q converting factor

  // we count for number of pulses in a separate function
  //CountingPulses();
}

// -- fit pulse shape using the LogNormal function
// not correct yet
void Waveform::FitPulseTime_LogNormalFunc(){
  // convert the waveform in time 
  // and invert so I have a positive pulse
  if(!hasPulse(m_nsigma)){
    cout<<"no pulse, skip .. "<<endl;
  }
  TH1D* hwave = getWaveformHistInverted();

  PulseShapeFunction_LogNormal func_LogNormal_obj;
  TF1* func1 = new TF1("exp_mod_LogNormal",func_LogNormal_obj,0,10000,5);
  func1->SetParNames("t0","mu","w","a","baseline");
  double t_peak = 2.0*hwave->GetMaximumBin();
  double t1 = t_peak - 5;
  //cout<<t_peak<<"\t"<<hwave->Integral( (t_peak-10)/2,(t_peak+20)/2)<<endl;
  func1->SetParameters(t_peak-5, 2, 1.0, 2.0*hwave->Integral( (t_peak-10)/2, (t_peak+20)/2 ), 0);
  hwave->Fit("exp_mod_LogNormal","RQ","",t_peak-15,t_peak+25);
  double par[4];
  for(int k=0; k<2; k++){
  for(int i=0; i<4; i++)
    {
      par[i] = func1->GetParameter(i);
      func1->SetParameter(i,par[i]);
    }
    hwave->Fit("exp_mod_LogNormal","RQ","",t_peak-15,t_peak+25);
  }
  //cout<<func1->GetChisquare()<<"\t"<<endl;
  hwave->GetXaxis()->SetRangeUser(420,520);

  m_pulse_shape_fit_chi2 = func1->GetChisquare();
  m_pulse_shape_fit_ndf = func1->GetNDF();
  m_pulse_shape_time = func1->GetParameter(0);
  m_pulse_shape_time_uncer = func1->GetParError(0);
  m_pulse_shape_w = func1->GetParameter(1);
  m_pulse_shape_w_uncer = func1->GetParError(1);
  m_pulse_shape_c = func1->GetParameter(2);
  m_pulse_shape_c_uncer = func1->GetParError(2);
  m_pulse_shape_a = func1->GetParameter(3);
  m_pulse_shape_a_uncer = func1->GetParError(3);

}

// -- fit pulse shape using the two exponential terms
void Waveform::FitPulseTime_TwoExpoentials(TFile* ofile){
  // convert the waveform in time 
  // and invert so I have a positive pulse
  if(!hasPulse(m_nsigma)){
    cout<<"no pulse, skip .. "<<endl;
  }
  TH1D* hwave = getWaveformHistInverted();

  PulseShapeFunction_TwoExpoentials func_obj;
  TF1* func1 = new TF1("pulse_func",func_obj,0,10000,5);
  func1->SetParNames("t0","#tau","#tau_{S}","a","baseline");
  double t_peak = 2.0*hwave->GetMaximumBin();
  double t1 = t_peak - 5;
  //cout<<t_peak<<"\t"<<hwave->Integral( (t_peak-10)/2,(t_peak+20)/2)<<endl;
  func1->SetParameters(t_peak-5, 3, 7, 2.0*hwave->Integral( (t_peak-10)/2, (t_peak+20)/2 ), 0);
  hwave->Fit("pulse_func","RQ","",t_peak-15,t_peak+25);
  double par[5];
  for(int k=0; k<2; k++){
  for(int i=0; i<5; i++)
    {
      par[i] = func1->GetParameter(i);
      func1->SetParameter(i,par[i]);
    }
    hwave->Fit("pulse_func","RQ","",t_peak-15,t_peak+25);
  }
  //cout<<func1->GetChisquare()<<"\t"<<endl;
  hwave->GetXaxis()->SetRangeUser(420,520);
  //ofile->cd();
  //hwave->Write();

  /*
  char name[60]; sprintf(name,"waveform_fit_twoExp/%s.png",hwave->GetName());
  TCanvas* c = new TCanvas();
  hwave->Draw();
  c->SaveAs(name);
  c->Close();
  */

  m_pulse_shape_fit_chi2 = func1->GetChisquare();
  m_pulse_shape_fit_ndf = func1->GetNDF();
  m_pulse_shape_time = func1->GetParameter(0);
  m_pulse_shape_time_uncer = func1->GetParError(0);
  // the parameters below are not correct for the moment,
  // and they are not used yet.
  m_pulse_shape_w = func1->GetParameter(1);
  m_pulse_shape_w_uncer = func1->GetParError(1);
  m_pulse_shape_c = func1->GetParameter(2);
  m_pulse_shape_c_uncer = func1->GetParError(2);
  m_pulse_shape_a = func1->GetParameter(3);
  m_pulse_shape_a_uncer = func1->GetParError(3);

}

// -- fit pulse shape using the exp-Gaus function
void Waveform::FitPulseTime_Laser(){
  // convert the waveform in time 
  // and invert so I have a positive pulse
  if(!hasPulse(m_nsigma)){
    cout<<"no pulse, skip .. "<<endl;
  }
  TH1D* hwave = getWaveformHistInverted();

  PulseShapeFunction_ExpGaus func_expGaus_obj;
  TF1* func1 = new TF1("exp_mod_Gaus",func_expGaus_obj,0,10000,4);
  func1->SetParNames("t0","w","c","a");
  double t_peak = 2.0*hwave->GetMaximumBin();
  double t1 = t_peak - 5;
  //cout<<t_peak<<"\t"<<hwave->Integral( (t_peak-10)/2,(t_peak+20)/2)<<endl;
  func1->SetParameters(t_peak-5, 2, 0.1, 2.0*hwave->Integral( (t_peak-10)/2, (t_peak+20)/2 ));
  hwave->Fit("exp_mod_Gaus","RQ","",t_peak-15,t_peak+25);
  double par[4];
  for(int k=0; k<2; k++){
  for(int i=0; i<4; i++)
    {
      par[i] = func1->GetParameter(i);
      func1->SetParameter(i,par[i]);
    }
    hwave->Fit("exp_mod_Gaus","RQ","",t_peak-15,t_peak+25);
  }
  //cout<<func1->GetChisquare()<<"\t"<<endl;
  hwave->GetXaxis()->SetRangeUser(420,520);

  m_pulse_shape_fit_chi2 = func1->GetChisquare();
  m_pulse_shape_fit_ndf = func1->GetNDF();
  m_pulse_shape_time = func1->GetParameter(0);
  m_pulse_shape_time_uncer = func1->GetParError(0);
  m_pulse_shape_w = func1->GetParameter(1);
  m_pulse_shape_w_uncer = func1->GetParError(1);
  m_pulse_shape_c = func1->GetParameter(2);
  m_pulse_shape_c_uncer = func1->GetParError(2);
  m_pulse_shape_a = func1->GetParameter(3);
  m_pulse_shape_a_uncer = func1->GetParError(3);
  //cout<< func1->GetMaximum()<<"\t"<<func1->GetMaximumX()<<"\t"<<func1->GetX( 0.1*func1->GetMaximum(), 400, m_pulse_shape_time )<<endl;
  //if(m_pulse_shape_time>550) cout<<"not a good fit "<<endl;
  //m_pulse_shape_time = func1->GetX( 0.1*func1->GetMaximum(), 400, m_pulse_shape_time  );
}

//------------------------------------------------------------------------------
void Waveform::CountingPulses(TFile* ofile){
  // to be developped
  // ideally we want to save time, amplitude, charge of each pulse in a waveform
  m_pulseChargeVec.resize(0); // store each pulse's charge
  m_pulseTimeVec.resize(0);   // store each pulse's start time
  m_pulseAmpVec.resize(0);    // store each pulse's amplitude
  m_pulseWidthVec.resize(0);  // store each pulse's width  
  m_npulse = 0;
  m_peakX.resize(0);

  TH1D* hwave = getWaveformHistInverted();
  hwave->SetXTitle("sample bins");
  double threshold =  baseline_width*m_nsigma;
  char title[100];
  sprintf(title, "%s threshold %.2f", hwave->GetName(), threshold);
  hwave->SetTitle(title);

  TSpectrum* aspec = new TSpectrum();

  if(hwave->GetMaximum() > threshold){
      double sigma = 5.0;
      int n = aspec->Search(hwave, sigma, "",  threshold*3.0/5.0/hwave->GetMaximum());
      m_npulse = aspec->GetNPeaks();
      double* peakX = aspec->GetPositionX();
      double* peakY = aspec->GetPositionY();
      TList *functions = hwave->GetListOfFunctions();
      //TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
      TPolyMarker* pm = new TPolyMarker();
      pm->SetMarkerColor(4);
      pm->SetMarkerStyle(8);
      functions->Add(pm);
      //cout<<"pulse finding "<<title<<"\t"
      //    <<hwave->GetMaximum()<<"\t"<<n<<"\t"<<m_npulse<<"..."<<endl;
      //for(int i=0; i<m_npulse; i++){
      //    cout<<m_peakX[i]<<"\t";
      // }
      //cout<<endl;

      // try to remove those single-bin spikes
      //double * newVec = new double [0];
      int newcnt=0;
      cout<<hwave->GetName()<<"\t"<<m_npulse<<endl;
      for(int i=0; i<m_npulse; i++){
          //cout<<i<<"\t"<<m_peakX[i]
          //    <<"\t"<<hwave->GetBinContent( int(peakX[i]/2) )
          //    <<"\t"<<hwave->GetBinContent( int(peakX[i]/2+1) )  
          //    <<"\t"<<hwave->GetBinContent( int(peakX[i]/2+2) )
          //    <<endl; 
          if( hwave->GetBinContent( int(peakX[i]/2) ) < threshold
             && hwave->GetBinContent( int(peakX[i]/2+2) ) < threshold ){
            pm->SetNextPoint(peakX[i], peakY[i]);
            //cout<<"remove this pulse"<<endl;
            continue;
          }
          newcnt++;
          m_peakX.push_back( peakX[i] );
      }
      m_npulse = newcnt;

      ofile->cd();
      hwave->Write();
  }
  else{
    m_npulse = 0;
  }
  
  //ofile->cd();
  //hwave->Write();
}

//------------------------------------------------------------------------------

bool Waveform::isValidWaveform()
{
  bool isValid;

  if(m_raw_waveform.size()==0){
    isValid=false;
  }
  else {
    isValid=true;
  }

  return isValid;
}

//------------------------------------------------------------------------------

  // Use fit to find the pulse instead

  // Better way to find the pulse

  // Noise filter ?

//------------------------------------------------------------------------------

Waveform::Complex_t Waveform::doFFT(Waveform::Waveform_t m_time_domain)
{
    Eigen::FFT<double> fft;
    Waveform::Complex_t  m_frequency_domain;
    fft.fwd(m_frequency_domain, m_time_domain);
    return m_frequency_domain;
}

TH1D* Waveform::getPowerSpectrum()
{
  Waveform::Complex_t m_spectrum = this->Waveform::doFFT(m_waveform);
  int m_fft_size = int(m_spectrum.size())/2;
  double max_sampling = m_sampling_freq/2;
  double freq_res = m_sampling_freq/m_spectrum.size();

  TH1D *h_power = new TH1D("", ";Frequency [MHz];Power",
                             m_spectrum.size(), 0, m_spectrum.size()*freq_res );

  for(int f=0; f<int(m_spectrum.size()); f++)
  {
    double ampls = pow(m_spectrum.at(f).real(), 2)
                                              + pow(m_spectrum.at(f).imag(), 2);
    h_power->Fill(f*freq_res, ampls );
  }

  return h_power;
}

//------------------------------------------------------------------------------

TH1D* Waveform::getWaveformHist()
{
  char hname[100];
  sprintf(hname, "Run%d-Subrun%d-Evt%d-Bd%d-Ch%d_hist", m_run,
                                         m_subrun, m_event, m_board, m_channel);

  TH1D *hist = new TH1D(hname, ";Time (ns);ADC", m_nsamples,
                                               0, m_nsamples*m_sampling_period);

  for(int t=0; t<m_nsamples; t++){ hist->Fill( t*m_sampling_period,
                                                           m_waveform.at(t) ); }

  return hist;
};

//------------------------------------------------------------------------------

TH1D* Waveform::getWaveformHistInverted()
{
  char hname[100];
  sprintf(hname, "Run%d-Subrun%d-Evt%d-Bd%d-Ch%d_hist", m_run,
                                         m_subrun, m_event, m_board, m_channel);

  TH1D *hist = new TH1D(hname, ";Time (ns);ADC", m_nsamples,
                                               0, m_nsamples*m_sampling_period);

  for(int t=0; t<m_nsamples; t++){ 
      //hist->Fill( t*m_sampling_period, -1.0*m_waveform.at(t) );
    hist->SetBinContent( t+1, -1.0*m_waveform.at(t) );  
  }

  return hist;
};

//------------------------------------------------------------------------------

TH1D* Waveform::getRawWaveformHist()
{
  char hname[100];
  sprintf(hname, "Run%d-Subrun%d-Evt%d-Bd%d-Ch%d_raw_hist", m_run,
                                         m_subrun, m_event, m_board, m_channel);

  TH1D *hist = new TH1D(hname, ";Time [ns];ADC", m_nsamples,
                                               0, m_nsamples*m_sampling_period);

  for(int t=0; t<m_nsamples; t++){ hist->Fill( t*m_sampling_period,
                                                       m_raw_waveform.at(t) ); }

  return hist;
};
