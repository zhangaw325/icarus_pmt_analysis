////////////////////////////////////////////////////////////////////////////////
// Macro to analyze the PMT noise
// Process the noise analyis for the filename given as argument. Edit the input
// filename with the file you prefer. Note: the output will always be in the format
// "noise_run%d_00%d.root". you must change manually the variables run and subrun for now.
//
// Usage: root -l loadLib.cc noise_ana.cc
// Out of the box, it will produce a root file in the macro folder called "noise_run1264_1.root"
// with the ADC distribution of all the pmts in data_dl1_run1264_1_20200227T235326-decoded.root
// after baseline subtraction.
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////
#include "Run.h"
#include "Waveform.h"
#include "Pmt.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TGraph.h"

using namespace std;

void charge_ana( string filename="run211_filenamelist.txt", string hvtablefile="run211_HV.txt" )
{

  // CHANGE THE FILENAME BELOW IF WANT TO USE ONE OTHER FILE FOR YOUR ANALYSIS
  fstream fin(filename.c_str(), ios::in);
  string firstfilename = "";
  getline(fin, firstfilename);
  fin.close();
  RUN my_run(firstfilename);
  int run=my_run.getRun();
  int subrun=my_run.getSubrun();

  const int nboards=1;
  const int nchannels=16;
  const double nsigma = 5.0;

  // set a hv table to the PMTs
  double hvtable[nboards][nchannels];
  fstream finHV(hvtablefile.c_str(), ios::in);
  double hv_value;
  for(int board=0; board<nboards; board++){
    for(int ch=0; ch<nchannels; ch++){
      finHV>>hv_value;
      hvtable[board][ch]=hv_value;
    }
  }
  finHV.close();

  //****************************************************************************
  // Input
  //

  TChain *tchain =  new TChain("caenv1730dump/events");

  // HERE you filename is loaded inside the macro to be processed!
  fstream fin1(filename.c_str());
  int nbOfFileCnt=0;
  while(fin1>>firstfilename){
    tchain->Add(firstfilename.c_str());
    nbOfFileCnt++;
  }
  fin1.close();
  cout<<"Processing run "<< run << ", with "<< nbOfFileCnt <<" files ..."<<endl;

  std::vector<std::vector<uint16_t> > *data=0; //unsigned short
  tchain->SetBranchAddress("fWvfmsVec", &data);


  //****************************************************************************
  // Now we create and initialize an histogram per each board and channel
  // This is an example, you can create more histograms using this model

  // baseline is really the basics for further processing the waveforms
  TH1D *h_pmt_rms[nboards][nchannels];
  TH1D *h_pmt_baseline[nboards][nchannels];

  // then we will want to look at waveforms characteristics:
  //     charge, amplitude, pulse time, number of pulse (for rate), etc.
  //TH1D *h_pmt_charge[nboards][nchannels]; 
  //TH1D* h_pmt_amplitude[nboards][nchannels];
 // TH1D* h_pmt_pulsetime[nboards][nchannels];
  //TH1D* h_pmt_NbOfPulses[nboards][nchannels];

  for(int board=0; board<nboards; board++)
  {
    for(int channel=0; channel<nchannels; channel++)
    {
      //DEFINE HERE THE HISTOGRAMS. EACH HISTOGRAM SHOULD HAVE A DIFFERENT NAME
      //TO AVOID CONFUSION WHEN YOU SAVE THEM TO FILE.
      char hname[100];
      char htitle[100];
      sprintf(hname, "hAllWaveBaseline_run%d_board%d_channel%d_rms",
                                                   run, board, channel);
      sprintf(htitle, "Waveform sample value distribution run%d board%d ch%d",
                                                   run, board, channel);
      h_pmt_rms[board][channel]= new TH1D(hname, htitle, 40, -20, 20);      
      h_pmt_rms[board][channel]->SetXTitle("Values (ADC count)");
      h_pmt_rms[board][channel]->SetYTitle("Counts");

      sprintf(hname, "hMeanBaseline_run%d_board%d_channel%d_baseline",
                                                   run, board, channel);
      sprintf(htitle, "Mean baselie of waveforms run%d board%d ch%d",
                                                   run, board, channel);
      h_pmt_baseline[board][channel]= new TH1D(hname, htitle, 40, 0, -1);
      h_pmt_baseline[board][channel]->SetBuffer(1);
      h_pmt_baseline[board][channel]->SetXTitle("Mean baseline of waveforms (ADC count)");
      h_pmt_baseline[board][channel]->SetYTitle("Count");
    }
  }

  //****************************************************************************
  // Here we create and initialize a vector with the class PMT. If you want you
  // can define also other objects (like TH1D classes) in here. See macro
  // noise_ana.cc with an example
  //
  PMT* pmts[nboards][nchannels];
  for( int board=0; board<nboards; board++ )
  {
    for( int channel=0; channel<nchannels; channel++ )
    {
      pmts[board][channel] = new PMT( run, board, channel );
      pmts[board][channel]->setHV(hvtable[board][channel]);
    }
  }

  //****************************************************************************
  // Now it is time to loop over the events
  //

  for(int e=0; e<tchain->GetEntries(); e++)
  {
    if(e%1000 ==0)
      cout << "         evt: " << e << endl;

    //if(e==3000) break;

    // WE TAKE THE EVENT
    tchain->GetEvent(e);

    for(int board=0; board<nboards; board++)
    {
      for(int channel=0; channel<nchannels; channel++)
      {

        // Create the WAVEFORM OBJECT
        Waveform *waveform = new Waveform(run, subrun ,e, board, channel, nsigma);

        // ONCE THE DATA ARE LOADED IN THE WAVEFORM OBJECT, BASELINE IS
        // AUTOMATICALLY SUBSTRACTED. YOU CAN USE THE FUNCTION "getWaveform()"
        // TO EXTRACT THE WAVEFORM VECTOR AFTER BASELINE SUBTRACTION OR
        // "getRawWaveform()" BEFORE BASELINE SUBTRACTION.
        waveform->loadData((*data).at(channel+nchannels*board));

        pmts[board][channel]->loadWaveform(waveform);

        if(!waveform->isValidWaveform()){ delete waveform; continue;  }

        // Loop over the entries of the waveform after the baseline subtraction
        // and fill the entries of the histogram h_pmt_rms

        for( float entry : waveform->getWaveform() )
        {
          h_pmt_rms[board][channel]->Fill( entry ); // this is distribution of baseline
        }

        // this is distribution of mean baseline from waveforms
        h_pmt_baseline[board][channel]->Fill( waveform->getBaselineMean() ); 

        delete waveform;

      } // channel
    } // boards
  } // event
  cout<<"      total events "<< tchain->GetEntries() <<endl;


  //****************************************************************************
  // Output: create an output TFile and write the histogram in it
  //

  char ofilename[100]; sprintf(ofilename, "Result_run%d_charge.root", run);
  TFile* ofile=new TFile(ofilename, "RECREATE"); //
  ofile->cd();

  TDirectory* baselineDir = ofile->mkdir("baselineDir");
  TDirectory* noiseDir = ofile->mkdir("noiseDir");
//  TDirectory* chargeDir = new TDirectory("chargeDir","chargeDirT");
//  TDirectory* timeDir = new TDirectory("timeDir","timeDirT");

  TGraph *g_rms = new TGraph( nboards*nchannels );
  char gname[100]; sprintf(gname, "g_rmsBaselineDist_run%d_", run);
  // this is the rms of the distribution of baseline, 
  // it is a measure of noise in the waveform
  g_rms->SetTitle(gname); g_rms->SetName(gname);

  for(int board=0; board<nboards; board++)
  {
    for(int channel=0; channel<nchannels; channel++)
    {
      //Write the histogram to the output TFILE

      //g_rms->SetPoint(channel+nchannels*board, channel+nchannels*board,
      //                                     h_pmt_rms[board][channel]->GetRMS());
      //baselineDir->cd();
      //h_pmt_rms[board][channel]->Write();
      //h_pmt_baseline[board][channel]->Write();
    }
  }

  noiseDir->cd();
  TCanvas* cNoise = new TCanvas();
  g_rms->Draw("ALP"); g_rms->SetMarkerStyle(7);
  cNoise->Write();
  cNoise->Close();

  TDirectory* ampDir = ofile->mkdir("AmplitudeDir");
  TDirectory* chargeDir = ofile->mkdir("ChargeDir");
  for( int board=0; board<nboards; board++ )
  {
    for( int channel=0; channel<nchannels; channel++ )
    {
      if( pmts[board][channel]->isHVon() ){
         
         //double this_hv_value = pmts[board][channel]->getHV();
         //cout<< "HV on for board "<<board<<", ch "<<channel<<": "<<this_hv_value<<endl;
         // Fit the charge distribution with the ideal response
         pmts[board][channel]->FitChargeWithIdeal();
      }
      pmts[board][channel]->writeHist(ofile, ampDir, chargeDir);
    }
  }  

  ofile->Close();

  cout << "All done" << endl;
  exit(0);

} //end macro
