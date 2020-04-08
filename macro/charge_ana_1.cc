#include "Run.h"
#include "Waveform.h"
#include "Pmt.h"
#include "CaliGainNovData.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TGraph.h"


void charge_ana_1(string metadatafile="metadata1.txt"){

  string path = "../../process20191121data/";

  const int nboards=1;
  const int nchannels=16;
  const double nsigma = 5.0;

  int run, nbsubrun;
  double hvtable[nchannels];

  fstream fin_meta(metadatafile.c_str(), ios::in);

  CaliGainNovData* allChargeHist[nboards][nchannels];// = new CaliGainNovData();  
  for(int bd=0; bd<nboards; bd++){
    for(int ch=0; ch<nchannels; ch++){
      allChargeHist[bd][ch] = new CaliGainNovData(bd,ch);
    }
  }

  while(fin_meta>>run){ // process one run
    
    fin_meta>>nbsubrun;
    // set the HV table
    for(int i=0; i<nchannels; i++){
      fin_meta>>hvtable[i]; // here hard coded already by assuming only one board 
    }

    // add all subruns into one 
    TChain *tchain =  new TChain("caenv1730dump/events");
    char filenamecontainer[100];
    for(int i=0; i<nbsubrun; i++){
      if(i+1 <10)
        sprintf(filenamecontainer, "%srun%d_0%d.root", path.c_str(), run, i+1) ;
      else sprintf(filenamecontainer, "%srun%d_%d.root", path.c_str(), run, i+1) ;
      //cout << filenamecontainer << endl;

      tchain->Add(filenamecontainer);
    }

    cout<<"Processing run "<< run << ", with "<< nbsubrun <<" files ..."<<endl;

    std::vector<std::vector<uint16_t> > *data=0; //unsigned short
    tchain->SetBranchAddress("fWvfmsVec", &data);    

    // prepare PMT containers
    PMT* pmts[nboards][nchannels];
    for( int board=0; board<nboards; board++ )
    {
      for( int channel=0; channel<nchannels; channel++ )
      {
        pmts[board][channel] = new PMT( run, board, channel );
        pmts[board][channel]->setHV(hvtable[channel]);
      }
    }

    // process the waveforms
    for(int e=0; e<tchain->GetEntries(); e++)
    {
      int subrun=1;
      if(e%1000 ==0){
        cout << "         evt: " << e << endl;
        subrun++;
      }
      if(e==2500) break;

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
          /*
          for( float entry : waveform->getWaveform() )
          {
            h_pmt_rms[board][channel]->Fill( entry ); // this is distribution of baseline
          }

          // this is distribution of mean baseline from waveforms
          h_pmt_baseline[board][channel]->Fill( waveform->getBaselineMean() ); 
          */
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

           allChargeHist[board][channel]->loadOneHist( pmts[board][channel]->getChargeHist(),  hvtable[channel] );

        }
        pmts[board][channel]->writeHist(ofile, ampDir, chargeDir);
      }
    }  

    ofile->Close();

    cout << "All done for run "<< run << endl;

  }//end run loop
  
  fin_meta.close();


  // test
  //TCanvas* cTest = new TCanvas();
  //allChargeHist[0][0]->plotHist();
  //cTest->SaveAs("testplot.png");

  // okay, let's do our global fit
  cout<<"Now we will fit all charge distributions \n"
      <<" and plot gain vs hv " <<endl;
  for(int board=0; board<nboards; board++){
    for(int ch=0; ch<nchannels; ch++){
      if( allChargeHist[board][ch]->getNbOfDataPoints() == 0 ) continue;
      cout <<" .. fitting board "<<board<<", ch "<<ch<<endl;

      allChargeHist[board][ch]->fitGainCurve( );
    }
  }

  cout<<"Success, exiting ..."<<endl;
  exit(0);
}