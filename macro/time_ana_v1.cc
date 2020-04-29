#include "Run.h"
#include "Waveform.h"
#include "Pmt.h"
#include "CaliGainNovData.h"
#include "IdealChargeFunction.h"
#include "ChargeFunctionWithPedestal.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMinuit.h"

const int Sum(const int number){
  int sum = 0;
  for(int i=0; i<number; i++)
    sum += i;
  return sum;
}

double myExpGaus(double t, double t0, double w, double c){
  return c/2.0*TMath::Exp(c*c*w*w/2.0)*TMath::Exp(-1.0*c*(t-t0))
          * TMath::Erfc( 1.0/1.414* (c*w-(t-t0)/w) );
}

double func_poissonExpGaus(double * x, double * par){
  double t0 = par[0];
  double w = par[1];
  double c = par[2];
  double a = par[3];
  double npe = par[4];

  double t = x[0];

  double value = 0;
  for(int k=1; k<100; k++){
    value += myExpGaus(t,t0,w*k,c)*TMath::Exp(-1.0*npe)*TMath::Power(npe,k)/TMath::Factorial(k);
  }
  return value*a;
}

double func_TwoExpOnExpGaus(double * x, double * par){
  double t0 = par[0];
  double 
}

double func_lognormal(double * x, double * par){

  double t0 = par[0];
  double mu = par[1];
  double w = par[2];
  double a = par[3];

  double t = x[0];
  if(t<=mu) return par[4]; //baseline
  else
  return par[4]+a/TMath::Sqrt(2.0*TMath::Pi()) / (w*(t-t0))
         * TMath::Exp( 1.0/2.0 * TMath::Power( (TMath::Log10(t-t0)-mu)/w ,2) );

}

double func_expGaus(double * x, double * par){
  double t0 = par[0];
  //double mu = par[1];
  double w = par[1];
  double c = par[2];
  double a = par[3];

  double t = x[0];

  return a*c/2.0*TMath::Exp(c*c*w*w/2.0)*TMath::Exp(-1.0*c*(t-t0))
          * TMath::Erfc( 1.0/1.414* (c*w-(t-t0)/w) );
}

void FittingPulseLogNormal(TFile* ofile, TH1D* hwave){
  ofile->cd();

  gStyle->SetOptFit(1111);

  TF1* func1 = new TF1("exp_mod_logNormal",func_lognormal,0,10000,4);
  func1->SetParNames("t0","mu","w","a");
  double t_peak = 2.0*hwave->GetMaximumBin();
  double t1 = t_peak - 5;
  //cout<<t_peak<<"\t"<<hwave->Integral( (t_peak-10)/2,(t_peak+20)/2)<<endl;
  func1->SetParameters(t_peak-5, 3, 1, 2.0*hwave->Integral( (t_peak-10)/2, (t_peak+20)/2 ));
  hwave->Fit("exp_mod_logNormal","RQ","",t_peak-15,t_peak+25);
  double par[4];
  for(int k=0; k<2; k++){
  for(int i=0; i<4; i++)
    {
      par[i] = func1->GetParameter(i);
      func1->SetParameter(i,par[i]);
    }
    hwave->Fit("exp_mod_logNormal","RQ","",t_peak-15,t_peak+25);
  }

  hwave->GetXaxis()->SetRangeUser(420,520);
  /*
  char name[50]; sprintf(name,"waveform_fit/%s.png",hwave->GetName());
  TCanvas* c = new TCanvas();
  hwave->Draw();
  c->SaveAs(name);
  c->Close();
  */
  hwave->Write();
}

void FittingPulseExpGaus(TFile* ofile, TH1D* hwave){
  ofile->cd();

  gStyle->SetOptFit(1111);

  TF1* func1 = new TF1("exp_mod_Gaus",func_expGaus,0,10000,4);
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

  hwave->GetXaxis()->SetRangeUser(420,520);
  /*
  char name[50]; sprintf(name,"waveform_fit/%s.png",hwave->GetName());
  TCanvas* c = new TCanvas();
  hwave->Draw();
  c->SaveAs(name);
  c->Close();
  */
  hwave->Write();
}

void FittingPulsePoissonExpaus(TFile* ofile, TH1D* hwave){
  ofile->cd();

  gStyle->SetOptFit(1111);

  TF1* func1 = new TF1("pulse_shape_func",func_poissonExpGaus,0,10000,5);
  func1->SetParNames("t0","w","c","a","npe");
  double t_peak = 2.0*hwave->GetMaximumBin();
  double npe = hwave->GetBinContent(t_peak/2)/40.0;
  double t1 = t_peak - 5;
  //cout<<t_peak<<"\t"<<hwave->Integral( (t_peak-10)/2,(t_peak+20)/2)<<endl;
  func1->SetParameters(t_peak-5, 2, 0.1, 2.0*hwave->Integral( (t_peak-10)/2, (t_peak+20)/2 ), npe);
  hwave->Fit("pulse_shape_func","RQ","",t_peak-15,t_peak+25);
  double par[5];
  for(int k=0; k<2; k++){
  for(int i=0; i<5; i++)
    {
      par[i] = func1->GetParameter(i);
      func1->SetParameter(i,par[i]);
    }
    hwave->Fit("pulse_shape_func","RQ","",t_peak-15,t_peak+25);
  }

  hwave->GetXaxis()->SetRangeUser(420,520);
  /*
  char name[50]; sprintf(name,"waveform_fit/%s.png",hwave->GetName());
  TCanvas* c = new TCanvas();
  hwave->Draw();
  c->SaveAs(name);
  c->Close();
  */
  hwave->Write();
}


//int main(int argc, char* argv[] ) {
void time_ana_v1(string metadatafile="metadata1.txt"){

  //string metadatafile="metadata1.txt"
  //string metadatafile (argv[1]);
  cout<<metadatafile<<endl;

  // turn off canvas display during running.
  // also stop showing statbox on histograms
  gROOT->SetBatch(kTRUE);
  //gStyle->SetOptStat(0);

  string path = "../../process20191121data/";

  const int nboards=1;
  const int nchannels=16;
  const double nsigma = 5.0;

  int run, nbsubrun;
  double hvtable[nchannels];

  fstream fin_meta(metadatafile.c_str(), ios::in);

  /*
  CaliGainNovData* allChargeHist[nboards][nchannels];// = new CaliGainNovData();  
  for(int bd=0; bd<nboards; bd++){
    for(int ch=0; ch<nchannels; ch++){
      allChargeHist[bd][ch] = new CaliGainNovData(bd,ch);
    }
  }
  */

  double pulsetime[nchannels];
  for(int i=0;i<nchannels;i++){
    pulsetime[i]=0.0;
  }
  const int nDiff = Sum(nchannels);//1+2+3+4;
  TH1D* hPulseTimeDiff[nDiff];
  char name[50];
  //for(int i=0; i<4; i++){
  //  sprintf(name,"hTimeDiff_ch%d_to_ch0",i+1);
  //  hPulseTimeDiff[i] = new TH1D(name,name,100,-5,5);
  //}
  int seqcnt=0;
  for(int i=0; i<nchannels-1; i++){
    for(int j=i+1; j<nchannels; j++){
      sprintf(name,"hTimeDiff_ch%d_to_ch%d",j, i);
    
      hPulseTimeDiff[seqcnt] = new TH1D(name,name,200,-10,10);
      hPulseTimeDiff[seqcnt]->SetXTitle("Time difference (ns)");
      hPulseTimeDiff[seqcnt]->SetYTitle("Counts");
      seqcnt++;
      //cout<<j+i-1<<endl;
    }
  } 

  while(fin_meta>>run){ // process one run

    //

    
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

    // prepare a root file to save the waveforms
    sprintf(filenamecontainer,"Waveforms_run%d.root", run);
    TFile* of_wave = new TFile(filenamecontainer,"recreate");

    // process the waveforms
    for(int e=0; e<tchain->GetEntries(); e++)
    {
      int subrun=1;
      if(e%1000 ==0){
        cout << "         evt: " << e << endl;
        subrun++;
      }
      if(e==10) break;

      // WE TAKE THE EVENT
      tchain->GetEvent(e);

      for(int board=0; board<nboards; board++)
      {
        for(int channel=0; channel<nchannels; channel++)
        {

          //if(channel!=8 && channel!=7) continue;

          // Create the WAVEFORM OBJECT
          Waveform *waveform = new Waveform(run, subrun ,e, board, channel, nsigma);

          // ONCE THE DATA ARE LOADED IN THE WAVEFORM OBJECT, BASELINE IS
          // AUTOMATICALLY SUBSTRACTED. YOU CAN USE THE FUNCTION "getWaveform()"
          // TO EXTRACT THE WAVEFORM VECTOR AFTER BASELINE SUBTRACTION OR
          // "getRawWaveform()" BEFORE BASELINE SUBTRACTION.
          waveform->loadData((*data).at(channel+nchannels*board));

          if(waveform->hasPulse(nsigma) ){
              TH1D* onewave = waveform->getWaveformHistInverted();
              FittingPulsePoissonExpaus(of_wave, onewave);
              // do pulse fitting here
              //FittingPulseExpGaus(of_wave, onewave);
              //FittingPulseLogNormal(of_wave, onewave);
              //waveform->FitPulseTime_Laser();
              //waveform->FitPulseTime_TwoExpoentials(of_wave);
              //pulsetime[channel] = waveform->GetPulseTimeFromFit();
              pmts[board][channel]->FillPulseTimeFromFit(waveform);
          }

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

        // since I have just one board, it is ok to process the time difference calculation here
        /*
        for(int ch=0; ch<4; ch++){
          hPulseTimeDiff[ch]->Fill(pulsetime[ch+1] - pulsetime[0] );
        }
        */

      } // boards
      seqcnt = 0;
      for(int ch=0; ch<nchannels-1; ch++) {
        for(int j=ch+1; j<nchannels; j++){
            hPulseTimeDiff[seqcnt]->Fill( pulsetime[j] - pulsetime[ch] );
            seqcnt++;
        }
      }
      for(int ch=0; ch<nchannels; ch++) pulsetime[ch] = 0.0;
    } // event
    cout<<"      total events "<< tchain->GetEntries() <<endl;

    of_wave->Close();

    //****************************************************************************
    // Output: create an output TFile and write the histogram in it
    //

    char ofilename[100]; sprintf(ofilename, "Result_run%d_charge.root", run);
    TFile* ofile=new TFile(ofilename, "RECREATE"); //
    ofile->cd();
    TDirectory* ampDir = ofile->mkdir("AmplitudeDir");
    TDirectory* chargeDir = ofile->mkdir("ChargeDir");
    TDirectory* timeDir = ofile->mkdir("timeDir");
    for( int board=0; board<nboards; board++ )
    {
      for( int channel=0; channel<nchannels; channel++ )
      {
        pmts[board][channel]->writeHist(ofile, ampDir, chargeDir, timeDir);
      }
    } 

    TDirectory* timeDiffDir = ofile->mkdir("TimeDiffDir"); 
    timeDiffDir->cd();
    seqcnt = 0;
    for(int ch=0; ch<nchannels-1; ch++){
      for(int j=ch+1; j<nchannels; j++){
          hPulseTimeDiff[seqcnt]->Write();
          seqcnt++;
      }
    }

    ofile->Close();

    cout << "All done for run "<< run << endl;

  }//end run loop
  
  fin_meta.close();

  cout<<"Success, exiting ..."<<endl;
  cout<<"what? can't stop"<<endl;
  //exit(0);
  //return 0;
}
