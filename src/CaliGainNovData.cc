#include "CaliGainNovData.h"
#include "IdealChargeFunction.h"
#include "ChargeFunctionWithPedestal.h"
//#include "Utils.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TSystem.h"

//using namespace std;

//global variables for hosting charge distribution
struct ds_global {
  vector<TH1D*> ghist;
  vector<double> gv;
  //static const int npars = 22;
  bool isLowChargeFlag = false;
  const int npars_highcharge = 10;
  const int npars_lowcharge = 22;
  int ndf=0;
} ds;

static CaliGainNovData* fitobj;

CaliGainNovData::CaliGainNovData(){
  fitobj=this;
  this->chargehist.resize(0);
  this->voltages.resize(0);
}

CaliGainNovData::CaliGainNovData(int board, int channel):m_board(board),m_ch(channel)
{
      chargehist.resize(0);
      lowchargehist.resize(0);
      voltages.resize(0);
}

CaliGainNovData::~CaliGainNovData(){
  for(int i=0; i<(int)chargehist.size(); i++){
    delete chargehist[i]; 
    delete lowchargehist[i];
  }
}

void CaliGainNovData::Initialize(vector<TH1D*>  hist, vector<double>  hv){
  cout<<"Initialize called "<< hist.size()<<"\t"<< hv.size() <<endl;
  for(int i=0; i< (int)hist.size(); i++){
    this->loadOneHist(hist[i], hv[i]);
    cout<<"testing static loading "<<hv[i]<<endl;
  }
  cout<<""<<endl;
}

void CaliGainNovData::loadOneHist(TH1D* hist, double hv){
    chargehist.push_back(hist);
    voltages.push_back(hv);
}

void CaliGainNovData::loadOneLowChargeHist(TH1D* hist){
    lowchargehist.push_back(hist);
}

void CaliGainNovData::plotHist(){
      cout<<chargehist.size()<<endl;
      
      char text[100]; 
      TCanvas* cTest = new TCanvas();
      gStyle->SetOptStat(0);
      for(int i=0; i<(int)chargehist.size(); i++){
        if(! chargehist[i]) continue;
        if(i==0) chargehist[i]->Draw();
        else chargehist[i]->Draw("same");
        chargehist[i]->SetLineColor(i+1);
        cout<<chargehist[i]->GetNbinsX()<<"\t"<<voltages[i]<<endl;
      }
      TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
      leg->SetHeader("HV values");
      for(int i=0; i<(int)chargehist.size(); i++){
        sprintf(text, "%.1f", voltages[i]);
        leg->AddEntry(chargehist[i],text,"lp");
      }
      leg->Draw();
      cTest->Modified(); cTest->Update();
      sprintf(text,"chargeDist_bd%d_ch%d.png",m_board, m_ch);
      cTest->SaveAs(text);
      cTest->Close();
      
}

// -------------------------------------------------------------------------------
// the hook for starting the global fit to charge distributions

void CaliGainNovData::fitGainCurve(TFile* ofile){
    if(chargehist.size()==0){ cout<<"fitGainCurve will not work if no data points"<<endl; exit(0); }

    // set the global variables to be used by the fitter
    // try to determinie whether the normal charge or lowcharge distributions 
    // to be used
    int bin1 = lowchargehist[0]->GetXaxis()->FindBin(-0.1);
    int bin2 = lowchargehist[0]->GetXaxis()->FindBin(0.1);
    double ratio  = lowchargehist[2]->Integral( bin1, bin2 )/lowchargehist[2]->Integral();
    cout<<"Ratio of events in pedestal (approx.) "<<ratio<<endl;
    bool flag_lowcharge = false;
    // if near 0 pC there are entries, use low charge
    // otherwise use normal charge
    if(ratio > 0.01){
      cout<<"fitting low charge distribution"<<endl;
      // use low charge hist in the fit
      flag_lowcharge = true;
      ds.isLowChargeFlag = true;
      for(int i=0; i< (int)lowchargehist.size(); i++){
          ds.ghist.push_back(lowchargehist[i]);
          ds.gv.push_back(voltages[i]);
      }
      fitLowCharge_fullExpr();
      plotGainCurveLowCharge(ofile,flag_lowcharge);
    }
    else{
      cout<<"fitting high charge distribution"<<endl;
      // use normal charge hist in the fit
      ds.isLowChargeFlag = false;
      for(int i=0; i< (int)chargehist.size(); i++){
          ds.ghist.push_back(chargehist[i]);
          ds.gv.push_back(voltages[i]);
      }
      fitHighCharge();
      plotGainCurveHighCharge(ofile,flag_lowcharge);
    }

    //fitChargeDistributions();

    // in the end of the fit function, need to reset the global variables
    ds.ghist.resize(0);
    ds.gv.resize(0);

    // plot gain curve
    //plotGainCurve(ofile, flag_lowcharge);
}

// plot gain curve (high charge)
void CaliGainNovData::plotGainCurveHighCharge(TFile* ofile, bool flagLowCharge){

    char gname[100]; 
    TCanvas* cGain = new TCanvas();
    sprintf(gname, "c_gain_bd%d_ch%d",m_board,m_ch);
    cGain->SetName(gname);
    TGraphErrors* gGain = new TGraphErrors();
    sprintf(gname, "gain_bd%d_ch%d",m_board,m_ch);
    gGain->SetName(gname);    gGain->SetTitle(gname);
    gGain->GetXaxis()->SetTitle("HV (V)");
    gGain->GetYaxis()->SetTitle("Gain (10^{7})");
    for(int i=0; i< (int)voltages.size(); i++){
        gGain->SetPoint(i, voltages[i], vresult[voltages.size()*i+1]/1.6);
        gGain->SetPointError(i, 0, verrors[voltages.size()*i+1]/1.6);
    }
    gGain->Draw("ALP");
    gGain->SetMarkerStyle(8);
    // fit the gain curve 
    TF1* fGain = new TF1("fGain","[0]*TMath::Power(x,[1])",1,2000);
    fGain->SetParNames("a","b");
    fGain->SetParameters(1e-24,7.5);
    //gGain->Fit("fGain","QR","",voltages[0]-10, voltages[2]+10);
    for(int i=0; i<5;i++){
      double a=fGain->GetParameter(0);
      double b=fGain->GetParameter(1);
      fGain->SetParameters(a,b);
      gGain->Fit("fGain","QR","",voltages[0]-10, voltages[2]+10);
    }
    ofile->cd();
    gGain->Write();
    cGain->Close();
    // plot charge distributions and the fit curves
    TCanvas* cCharge = new TCanvas();
    sprintf(gname,"cChargeDis_bd%d_ch%d",m_board,m_ch);
    cCharge->SetName(gname);
    IdealChargeFunction idealchargePlot_obj; // hard coded
    TF1* func[3];
    for(int i=0; i<(int)chargehist.size(); i++){
      if(i==0) {
          chargehist[i]->Draw(); 
      }
      else {
          chargehist[i]->Draw("same");
      }
      chargehist[i]->SetLineColor(i+1);
      chargehist[i]->Write();

      sprintf(gname,"idealfunc%d",i);
      func[i] = new TF1(gname,idealchargePlot_obj,0, 100, 4);
      //func[i]->SetParNames("mu","q0","sigma","amp.");
      func[i]->SetParameters(vresult[0],
                            vresult[i*3+1],
                            vresult[i*3+2],
                            vresult[i*3+3]);

      func[i]->SetNpx(600);
      func[i]->SetLineColor(1); func[i]->SetLineStyle(2);
      func[i]->Draw("same");
    }
    TLegend* leg = new TLegend(0.75,0.3,0.9,0.5);
    leg->SetHeader("HV values");
    for(int i=0; i<(int)chargehist.size(); i++){
      sprintf(gname, "%.1f", voltages[i]);
      leg->AddEntry(chargehist[i],gname,"lp");
    }
    leg->Draw(); 
    TPaveText* pt = new TPaveText(0.5,0.5,0.9,0.9, "brNDC");   
    sprintf(gname, "Fit status: %d, #chi^{2}/ndf: %.1f/%d", m_fitstatus, minChi2, ds.ndf-ds.npars_highcharge);
    pt->AddText(gname);
    sprintf(gname,"npe: %.2f #pm %.2f", vresult[0], verrors[0]);
    pt->AddText(gname);
    pt->AddText("# , q , w , a:");
    for(int i=0; i<(int)chargehist.size();i++){
      sprintf(gname, "%d, %.2f #pm %.2f, %.2f #pm %.2f, %.0f #pm %.0f", i+1,
                 vresult[i*chargehist.size()+1], verrors[i*chargehist.size()+1],
                 vresult[i*chargehist.size()+2], verrors[i*chargehist.size()+2],
                 vresult[i*chargehist.size()+3], verrors[i*chargehist.size()+3]);
      pt->AddText(gname);
    }
    pt->Draw();
    pt->SetBorderSize(1);
    pt->SetFillColor(0);
    pt->SetTextFont(42);
    pt->SetTextAlign(12);

    cCharge->Modified(); cCharge->Update();
    cCharge->Write();
    cCharge->Close();
  
}

// plot gain curve (high charge)
void CaliGainNovData::plotGainCurveLowCharge(TFile* ofile, bool flagLowCharge){

    char gname[100]; 
    TCanvas* cGain = new TCanvas();
    sprintf(gname, "c_gain_bd%d_ch%d",m_board,m_ch);
    cGain->SetName(gname);
    TGraphErrors* gGain = new TGraphErrors();
    sprintf(gname, "gain_bd%d_ch%d",m_board,m_ch);
    gGain->SetName(gname);    gGain->SetTitle(gname);
    gGain->GetXaxis()->SetTitle("HV (V)");
    gGain->GetYaxis()->SetTitle("Gain (10^{7})");
    for(int i=0; i< (int)voltages.size(); i++){
        gGain->SetPoint(i, voltages[i], vresult[i*7+1]/1.6);
        gGain->SetPointError(i, 0, verrors[i*7+2]/1.6);
    }
    gGain->Draw("ALP");
    gGain->SetMarkerStyle(8);
    // fit the gain curve 
    TF1* fGain = new TF1("fGain","[0]*TMath::Power(x,[1])",1,2000);
    fGain->SetParNames("a","b");
    fGain->SetParameters(1e-24,7.5);
    //gGain->Fit("fGain","QR","",voltages[0]-10, voltages[2]+10);
    for(int i=0; i<5;i++){
      double a=fGain->GetParameter(0);
      double b=fGain->GetParameter(1);
      fGain->SetParameters(a,b);
      gGain->Fit("fGain","QR","",voltages[0]-10, voltages[2]+10);
    }
    ofile->cd();
    gGain->Write();
    cGain->Close();
    // plot charge distributions and the fit curves
    TCanvas* cCharge = new TCanvas();
    sprintf(gname,"cChargeDis_bd%d_ch%d",m_board,m_ch);
    cCharge->SetName(gname);
    myChargeFunctionWithRealPedestal myChargeFunctionWithRealPedestal_obj; // hard coded  
    TF1* func[3];
    for(int i=0; i<(int)lowchargehist.size(); i++){
      if(i==0) {
        lowchargehist[i]->Draw();
      }
      else {
        lowchargehist[i]->Draw("same");
      }
      lowchargehist[i]->SetLineColor(i+1);
      lowchargehist[i]->Write();

      sprintf(gname,"myresponsefunc%d",i);
      func[i] = new TF1(gname,myChargeFunctionWithRealPedestal_obj,-1, 50, 8);
      func[i]->SetParameter(0, vresult[0]);
      for(int m=0; m<7; m++){
        func[i]->SetParameter(m+1, vresult[i*7+m+1]);
      }
      func[i]->SetNpx(600);
      func[i]->SetLineColor(1); func[i]->SetLineStyle(2);
      func[i]->Draw("same");
    }
    TLegend* leg = new TLegend(0.8,0.3,0.9,0.5);
    leg->SetHeader("HV values");
    for(int i=0; i<(int)lowchargehist.size(); i++){
      sprintf(gname, "%.1f", voltages[i]);
      leg->AddEntry(lowchargehist[i],gname,"lp");
    }
    leg->Draw(); 
    TPaveText* pt = new TPaveText(0.5,0.5,0.9,0.9, "brNDC");   
    sprintf(gname, "Fit status: %d, #chi^{2}/ndf: %.1f/%d", m_fitstatus, minChi2, ds.ndf-ds.npars_lowcharge);
    pt->AddText(gname);
    sprintf(gname,"npe: %.2f #pm %.2f", vresult[0], verrors[0]);
    pt->AddText(gname);
    pt->AddText("# , q , #sigma :");
    for(int i=0; i<(int)lowchargehist.size();i++){
      sprintf(gname, "%d, %.2f #pm %.2f, %.2f #pm %.2f", i+1,
                 vresult[i*7+1], verrors[i*7+1],
                 vresult[i*7+2], verrors[i*7+2]);
      pt->AddText(gname);
    }
    pt->Draw();
    pt->SetBorderSize(1);
    pt->SetFillColor(0);
    pt->SetTextFont(42);
    pt->SetTextAlign(12);

    TPaveText* pt1 = new TPaveText(0.5, 0.3, 0.8, 0.5, "brNDC");
    pt1->AddText("other parameters");
    for(int i=0; i<(int)chargehist.size();i++){
        sprintf(gname, "%.2f #pm %.2f, %.2f #pm %.2f, %.2f #pm %.2f, %.2f #pm %.2f, %.2f #pm %.2f ",
                   vresult[i*7+3], verrors[i*7+3],
                   vresult[i*7+4], verrors[i*7+4],
                   vresult[i*7+5], verrors[i*7+5],
                   vresult[i*7+6], verrors[i*7+6],
                   vresult[i*7+7], verrors[i*7+7]);
        pt1->AddText(gname);
    }
    pt1->Draw();
    pt1->SetBorderSize(1);
    pt1->SetFillColor(0);
    pt1->SetTextFont(42);
    pt1->SetTextAlign(12);

    cCharge->Modified(); cCharge->Update();
    cCharge->Write();
    cCharge->Close();
  
}

// fit charge distributions
void CaliGainNovData::fitLowCharge_fullExpr(){
    initializeParVectors(ds.npars_lowcharge);
    const int nVs = 3;
    const int npars=22;
    double vstart[npars];
    double step[npars];
    double minVal[npars]={0.01, 
                           0.01, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0,
                           0.01, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0,
                           0.01, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0};
    double maxVal[npars]={50, 
                           10, 1, 0.5, 1, 1e5, 1, 1e7,
                           10, 1, 0.5, 1, 1e5, 1, 1e7,
                           10, 1, 0.5, 1, 1e5, 1, 1e7};
    string parName[npars]={"npe",
          "q11", "#sigma_{11}", "q01", "#sigma_{01}", "c01", "w1","a1", 
          "q12", "#sigma_{12}", "q02", "#sigma_{02}", "c02", "w2","a2", 
          "q13", "#sigma_{13}", "q03", "#sigma_{03}", "c03", "w3","a3"};

    // initial values to the parameters
    vstart[0] = ds.ghist[2]->GetMean()/1.6;
    if(vstart[0]<0.5) vstart[0]=0.5;
    step[0] = 0.05;
    for(int i=0; i<nVs; i++){
        step[i*7+1] = 0.05; // step for q1
        step[i*7+2] = 0.02; // step for sigma1
        step[i*7+3] = 0.01; // step for q0
        step[i*7+4] = 0.01; // step for sigma0
        step[i*7+5] = 0.01; // step for c0;
        step[i*7+6] = 0.05; // step for w;
        step[i*7+7] = 20; // step for a

        vstart[i*7+1] = 1; // start for q1
        vstart[i*7+2] = 0.5; // start for sigma1
        vstart[i*7+3] = 0; // start for q0
        vstart[i*7+4] = 0.1; // start for sigma0
        vstart[i*7+5] = 10; // start for c0
        vstart[i*7+6] = 0.3; // start for w
        vstart[i*7+7] = ds.ghist[i]->Integral()*TMath::Exp(-1.0*vstart[0]); // start for a
    }
    
    // set the fitter: use TMinuit
    //gSystem->Load("libMinuit");
    TMinuit* minimizer = new TMinuit(npars); // was ds.npars here
    // the fit function is set
    //cout<<"Set the minimizer with "<<npar<<" parameters "<<endl;
    minimizer->SetFCN(CaliGainNovData::mygainfcn);
    //  
    for(int i=0; i<npars; i++){
      minimizer->DefineParameter(i, parName[i].c_str(), vstart[i], step[i], 
                             minVal[i], maxVal[i]);
    }
    // do the fit
    minimizer->Command("SET PRINT -1");//
    Int_t ierrs=0;
    minimizer->mnexcm("SET NOWarnings", 0, 0, ierrs);
    int fitstatus= minimizer->Migrad();

    // get parameters
    for(int i=0; i<npars; i++){
      minimizer->GetParameter(i, vresult[i], verrors[i]);
    }
    cout<<"check fit status "<<fitstatus<<endl;
    // if fit is unsuccessful, try again with updated initial values
    if(fitstatus==4){
      for(int i=0; i<npars; i++){
        minimizer->DefineParameter(i, parName[i].c_str(), vresult[i], step[i]/2.0, 
                               minVal[i], maxVal[i]);
      }
      fitstatus = minimizer->Migrad();
      cout<<"status of the second fit "<<fitstatus<<endl;
      for(int i=0; i<npars; i++){
        minimizer->GetParameter(i, vresult[i], verrors[i]);
      }
    }  

    // here we get the result from the fit
    double amin,edm,errdef;
    int nvpar,nparx,icstat;
    minimizer->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    cout<<amin<<"\t"<<edm<<"\t"<<errdef<<"\t"<<nvpar<<"\t"<<nparx<<"\t"<<icstat<<endl;  
    cout<<"---fit results "<<endl;
    for(int i=0; i<npars; i++){
      cout<<parName[i]<<"\t"<<vresult[i]<<"\t"<<verrors[i]<<endl;
    }

    m_fitstatus = fitstatus;
    minChi2 = amin;
}


/*
// fit high charge distributions
void CaliGainNovData::fitLowCharge(){


    cout<<"fitting low charge distributions..."<<endl;

    // set the vectors to hold fit results
    initializeParVectors(ds.npars_lowcharge);

    // define number of parameters and their initial values

    int bin1[3];
    int bin2[3];
    for(int i=0;i<3;i++){
      bin1[i] = lowchargehist[i]->GetXaxis()->FindBin(-0.2);
      bin2[i] = lowchargehist[i]->GetXaxis()->FindBin(0.2);
    }
    //int bin3 = lowchargehist[0]->GetXaxis()->FindBin(0.5);
    //int bin4 = lowchargehist[0]->GetNbinsX();


    int n = getNbOfDataPoints() ;
    const int n_pars = 19;//ds.npars_lowcharge;//10;//1 + 3*n ;
    double vstart[n_pars]={20, 
                       1.0, 0.5, lowchargehist[0]->Integral()/400.0*40,
                       1.0, 0.5, lowchargehist[1]->Integral()/400.0*40,
                       1.0, 0.5, lowchargehist[2]->Integral()/400.0*40,
                       0.0, 0.1, lowchargehist[0]->Integral(bin1[0],bin2[0])/10.0,
                       0.0, 0.1, lowchargehist[1]->Integral(bin1[1],bin2[1])/10.0,
                       0.0, 0.1, lowchargehist[2]->Integral(bin1[2],bin2[2])/10.0 } ;
    double step[n_pars]={0.02,
                          0.05, 0.1, 5, 0.05, 0.1, 5, 0.05, 0.1, 5, 
                          0.005, 0.01, 10, 0.005, 0.01, 10, 0.005, 0.01, 10};
    double minVal[n_pars]={0,
                           0.01, 0.001, 0, 0.01, 0.001, 0, 0.01, 0.001, 0,
                           -0.5, 0.001, 0, -0.5, 0.001, 0, -0.5, 0.001, 0};
    double maxVal[n_pars]={10, 10, 10,lowchargehist[0]->Integral()/4.0, 
                               10, 10,lowchargehist[1]->Integral()/4.0, 
                               10, 10,lowchargehist[2]->Integral()/4.0,
                               0.5, 0.5, lowchargehist[0]->Integral(bin1[0],bin2[0])*10.0,
                               0.5, 0.5, lowchargehist[1]->Integral(bin1[1],bin2[1])*10.0,
                               0.5, 0.5, lowchargehist[2]->Integral(bin1[2],bin2[2])*10.0 };
    string parName[n_pars]={"npe",
                            "q1", "w1", "a1",
                            "q2", "w2", "a2",
                            "q3", "w3", "a3",
                            "pedm1", "pedw1", "peda1",
                            "pedm2", "pedw2", "peda2",
                            "pedm3", "pedw3", "peda3"
                           };

    // reset some of the initial parameters
    for(int i=0; i< (int)lowchargehist.size(); i++){
      //step[i*3+1] = 0.05;
      vstart[0] = lowchargehist[i]->GetMean()/1.6;
      //if(vstart[0]<1) vstart[0]=1;
    }

    // set the fitter: use TMinuit
    TMinuit* minimizer = new TMinuit(n_pars);
    // the fit function is set
    minimizer->SetFCN(CaliGainNovData::myfcnWithPedestal);
    //  
    for(int i=0; i<n_pars; i++){
      minimizer->DefineParameter(i, parName[i].c_str(), vstart[i], step[i], 
                             minVal[i], maxVal[i]);
    }
    // do the fit
    minimizer->Command("SET PRINT -1");//
    //minimizer->SetPrintLevel(0);
    Int_t ierrs=0;
    minimizer->mnexcm("SET NOWarnings", 0, 0, ierrs);
    int fitstatus= minimizer->Migrad();

    // get parameters
    for(int i=0; i<n_pars; i++){
      minimizer->GetParameter(i, vresult[i], verrors[i]);
    }
    cout<<"check fit status "<<fitstatus<<endl;
    if(fitstatus==4){
      for(int i=0; i<n_pars; i++){
        minimizer->DefineParameter(i, parName[i].c_str(), vresult[i], step[i]/2.0, 
                               minVal[i], maxVal[i]);
      }
      fitstatus = minimizer->Migrad();
      cout<<"status of the second fit "<<fitstatus<<endl;
      for(int i=0; i<n_pars; i++){
        minimizer->GetParameter(i, vresult[i], verrors[i]);
      }
    }  

    // here we get the result from the fit
    double amin,edm,errdef;
    int nvpar,nparx,icstat;
    minimizer->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    cout<<amin<<"\t"<<edm<<"\t"<<errdef<<"\t"<<nvpar<<"\t"<<nparx<<"\t"<<icstat<<endl;  
    cout<<"---fit results "<<endl;
    for(int i=0; i<n_pars; i++){
      cout<<parName[i]<<"\t"<<vresult[i]<<"\t"<<verrors[i]<<endl;
    }

    m_fitstatus = fitstatus;
    minChi2 = amin;

}
*/

/*
// fit high charge distributions
void CaliGainNovData::fitLowCharge_expGaus(){ // pedestal with expoentially-modified Guassian


    cout<<"fitting low charge distributions (ped w/ exp-mod Gaus)..."<<endl;

    // set the vectors to hold fit results
    initializeParVectors(ds.npars_lowcharge);

    // define number of parameters and their initial values

    int bin1[3];
    int bin2[3];
    for(int i=0;i<3;i++){
      bin1[i] = lowchargehist[i]->GetXaxis()->FindBin(-0.2);
      bin2[i] = lowchargehist[i]->GetXaxis()->FindBin(0.2);
    }
    //int bin3 = lowchargehist[0]->GetXaxis()->FindBin(0.5);
    //int bin4 = lowchargehist[0]->GetNbinsX();


    int n = getNbOfDataPoints() ;
    const int n_pars = 22;//ds.npars_lowcharge;//10;//1 + 3*n ;
    double vstart[n_pars]={20, 
                       1.0, 0.5, lowchargehist[0]->Integral()/400.0*40,
                       1.0, 0.5, lowchargehist[1]->Integral()/400.0*40,
                       1.0, 0.5, lowchargehist[2]->Integral()/400.0*40,
                       0.0, 0.1, lowchargehist[0]->Integral(bin1[0],bin2[0])/10.0*4, 0.1,
                       0.0, 0.1, lowchargehist[1]->Integral(bin1[1],bin2[1])/10.0*4, 0.1,
                       0.0, 0.1, lowchargehist[2]->Integral(bin1[2],bin2[2])/10.0*4, 0.1};
    double step[n_pars]={0.02,
                          0.05, 0.1, 5, 
                          0.05, 0.1, 5, 
                          0.05, 0.1, 5, 
                          0.005, 0.01, 10, 0.01,
                          0.005, 0.01, 10, 0.01,
                          0.005, 0.01, 10, 0.01};
    double minVal[n_pars]={0,
                           0.01, 0.001, 0, 
                           0.01, 0.001, 0, 
                           0.01, 0.001, 0,
                           -0.5, 0.001, 0, 0,
                           -0.5, 0.001, 0, 0,
                           -0.5, 0.001, 0, 0};
    double maxVal[n_pars]={10, 10, 10,lowchargehist[0]->Integral()/4.0, 
                               10, 10,lowchargehist[1]->Integral()/4.0, 
                               10, 10,lowchargehist[2]->Integral()/4.0,
                               0.5, 0.5, lowchargehist[0]->Integral(bin1[0],bin2[0])*10.0, 0.5,
                               0.5, 0.5, lowchargehist[1]->Integral(bin1[1],bin2[1])*10.0, 0.5,
                               0.5, 0.5, lowchargehist[2]->Integral(bin1[2],bin2[2])*10.0, 0.5 };
    string parName[n_pars]={"npe",
                            "q1", "w1", "a1",
                            "q2", "w2", "a2",
                            "q3", "w3", "a3",
                            "pedm1", "pedw1", "peda1", "pedt1",
                            "pedm2", "pedw2", "peda2", "pedt2",
                            "pedm3", "pedw3", "peda3", "pedt3"
                           };

    // reset some of the initial parameters
    for(int i=0; i< (int)lowchargehist.size(); i++){
      //step[i*3+1] = 0.05;
      vstart[0] = lowchargehist[i]->GetMean()/1.6;
      //if(vstart[0]<1) vstart[0]=1;
    }

    // set the fitter: use TMinuit
    TMinuit* minimizer = new TMinuit(n_pars);
    // the fit function is set
    minimizer->SetFCN(CaliGainNovData::myfcnWithPedestal_expGaus);
    //  
    for(int i=0; i<n_pars; i++){
      minimizer->DefineParameter(i, parName[i].c_str(), vstart[i], step[i], 
                             minVal[i], maxVal[i]);
    }
    // do the fit
    minimizer->Command("SET PRINT -1");//
    //minimizer->SetPrintLevel(0);
    Int_t ierrs=0;
    minimizer->mnexcm("SET NOWarnings", 0, 0, ierrs);
    int fitstatus= minimizer->Migrad();

    // get parameters
    for(int i=0; i<n_pars; i++){
      minimizer->GetParameter(i, vresult[i], verrors[i]);
    }
    cout<<"check fit status "<<fitstatus<<endl;
    if(fitstatus==4){
      for(int i=0; i<n_pars; i++){
        minimizer->DefineParameter(i, parName[i].c_str(), vresult[i], step[i]/2.0, 
                               minVal[i], maxVal[i]);
      }
      fitstatus = minimizer->Migrad();
      cout<<"status of the second fit "<<fitstatus<<endl;
      for(int i=0; i<n_pars; i++){
        minimizer->GetParameter(i, vresult[i], verrors[i]);
      }
    }  

    // here we get the result from the fit
    double amin,edm,errdef;
    int nvpar,nparx,icstat;
    minimizer->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    cout<<amin<<"\t"<<edm<<"\t"<<errdef<<"\t"<<nvpar<<"\t"<<nparx<<"\t"<<icstat<<endl;  
    cout<<"---fit results "<<endl;
    for(int i=0; i<n_pars; i++){
      cout<<parName[i]<<"\t"<<vresult[i]<<"\t"<<verrors[i]<<endl;
    }

    m_fitstatus = fitstatus;
    minChi2 = amin;

}
*/

// fit high charge distributions
void CaliGainNovData::fitHighCharge(){

    //cout<<"fitting high charge distributions..."<<endl;

    // set the vectors to hold fit results
    initializeParVectors(ds.npars_highcharge);

    // define number of parameters and their initial values
    int n = getNbOfDataPoints() ;
    const int n_pars = 10;//ds.npars_highcharge;//10;//1 + 3*n ;
    double vstart[n_pars]={20, 
                       1.0, 0.5, chargehist[0]->Integral()/2.0,
                       1.0, 0.5, chargehist[1]->Integral()/2.0,
                       1.0, 0.5, chargehist[2]->Integral()/2.0} ;
    double step[n_pars]={0.1, 0.1, 0.1, 10, 0.1, 0.1, 10, 0.1, 0.1, 10};
    double minVal[n_pars]={0, 0.01, 0.001, 0, 0.01, 0.001, 0, 0.01, 0.001, 0};
    double maxVal[n_pars]={50, 10, 10,chargehist[0]->Integral()*2.0, 
                          10, 10,chargehist[1]->Integral()*2.0, 
                          10, 10,chargehist[2]->Integral()*2.0};
    string parName[n_pars]={"npe",
                            "q1", "w1", "a1",
                            "q2", "w2", "a2",
                            "q3", "w3", "a3"
                           };

    //for(int i=0; i<(int)chargehist.size(); i++){
      vstart[0] = chargehist[2]->GetMean()/1.6;
    //}

    // set the fitter: use TMinuit
    //gSystem->Load("libMinuit");
    TMinuit* minimizer = new TMinuit(n_pars);
    // the fit function is set
    minimizer->SetFCN(CaliGainNovData::myfcn);
    //  
    for(int i=0; i<n_pars; i++){
      minimizer->DefineParameter(i, parName[i].c_str(), vstart[i], step[i], 
                             minVal[i], maxVal[i]);
    }
    // do the fit
    minimizer->Command("SET PRINT -1");//SetPrintLevel(-1);
    int fitstatus= minimizer->Migrad();

    // get parameters
    for(int i=0; i<n_pars; i++){
      minimizer->GetParameter(i, vresult[i], verrors[i]);
    }
    cout<<"check fit status "<<fitstatus<<endl;
    if(fitstatus==4){
      for(int i=0; i<n_pars; i++){
        minimizer->DefineParameter(i, parName[i].c_str(), vresult[i], step[i]/2.0, 
                               minVal[i], maxVal[i]);
      }
      fitstatus = minimizer->Migrad();
      cout<<"status of the second fit "<<fitstatus<<endl;
      for(int i=0; i<n_pars; i++){
        minimizer->GetParameter(i, vresult[i], verrors[i]);
      }
    }  

    // here we get the result from the fit
    double amin,edm,errdef;
    int nvpar,nparx,icstat;
    minimizer->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    cout<<amin<<"\t"<<edm<<"\t"<<errdef<<"\t"<<nvpar<<"\t"<<nparx<<"\t"<<icstat<<endl;  

    cout<<"---fit results "<<endl;
    for(int i=0; i<n_pars; i++){
      cout<<parName[i]<<"\t"<<vresult[i]<<"\t"<<verrors[i]<<endl;
    }

    m_fitstatus = fitstatus;
    minChi2 = amin;

}

double CaliGainNovData::myGausPDF(double x, double mean, double sigma){
  // test output
  double value = 1.0/TMath::Sqrt(2.0*TMath::Pi()*sigma*sigma)*TMath::Exp(-1.0/2.0*TMath::Power( (x-mean)/sigma ,2));

  //if( TMath::IsNaN(value) ){
  //  cout<<"nan in myGausPDF"<<endl;
  //  cout<<x<<"\t"<<mean<<"\t"<<sigma<<endl;
  //}

  return value;
}

double CaliGainNovData::myExpGausPDF(double x, double mean, double sigma, double c0){
  double value = sigma/c0 * TMath::Sqrt(TMath::Pi()/2.0)
         * TMath::Exp( 1./2.*(sigma/c0)*(sigma/c0)- (x-mean)/c0 )
         * TMath::Erfc( 1./TMath::Sqrt(2.) * (sigma/c0 - (x-mean)/sigma) );
  // test output
  //if( TMath::IsNaN(value) ){
  //  cout<<"nan in myExpGausPDF"<<endl;
  //  cout<<x<<"\t"<<mean<<"\t"<<sigma<<"\t"<<c0<<endl;
  //}
  return value;
  //return c0*0.5*TMath::Exp(0.5*c0*(2.0*mean+c0*sigma*sigma-2.0*x))
  //       * TMath::Erfc( (mean+c0*sigma*sigma-x)/sigma/1.414);
}

// --------------------------------------------------------------------
// using the full expression derived by Milind
//double CaliGainNovData::mygainResponse(double * par){
double CaliGainNovData::mygainResponse(double p0,
                                       double p1, double p2, double p3, double p4, double p5, double p6, double p7,
                                       double p8, double p9, double p10, double p11, double p12, double p13, double p14,
                                       double p15, double p16, double p17, double p18, double p19, double p20, double p21){
  // following Milind's paper, arXiv: 1909.05373
  const int nVs = 3; // number of voltage points
  int n=7; // for each charge distribution, number of pars (not counting npe)
  double mu = p0;//par[0];  // mean npe
  double q[nVs]; //= par[1];   // charge of spe
  double sigma[nVs];// = par[2]; // width of spe
  double q0[nVs]; //= par[3];   // mean of backgroud
  double sigma0[nVs]; //= par[4]; // width of backgroud
  double c0[nVs]; //= par[5];    // exp. component in the backgroud
  double w[nVs]; //= par[6]; // the weight: probability of dark pulse 
  double a[nVs]; //= par[7]; // final scale factor or amplitude
  /*
  for(int i=0; i<nVs;i++){
    q[i] = 0;//par[i*n+1];
    sigma[i] = 0;//par[i*n+2];
    q0[i] = 0;//par[i*n+3];
    sigma0[i] = 0;//par[i*n+4];
    c0[i] = 0;//par[i*n+5];
    w[i] = 0;//par[i*n+6];
    a[i] = 0;//par[i*n+7];
  }
  */

  q[0] = p1;     q[1] = p8;     q[2] = p15;
  sigma[0]=p2;   sigma[1]=p9;   sigma[2]=p16;
  q0[0]=p3;      q0[1]=p10;     q0[2]=p17;
  sigma0[0]=p4;  sigma0[1]=p11; sigma0[2]=p18;
  c0[0]=p5;      c0[1]=p12;     c0[2]=p19;
  w[0]=p6;       w[1]=p13;      w[2]=p20;
  a[0]=p7;       a[1]=p14;      a[2]=p21;

//  cout<<""<<endl;

//  if(TMath::IsNaN(mu)){
//    cout<<"iniital parameters to use (nan?): (note isLowChargeFlag)"<<": "<<ds.isLowChargeFlag<<endl;
//    cout<<p0<<endl;
//    cout<<p1<<"\t"<<p2<<"\t"<<p3<<"\t"<<p4<<"\t"<<p5<<"\t"<<p6<<"\t"<<p7<<endl;
//    cout<<p8<<"\t"<<p9<<"\t"<<p10<<"\t"<<p11<<"\t"<<p12<<"\t"<<p13<<"\t"<<p14<<endl;
//    cout<<p15<<"\t"<<p16<<"\t"<<p17<<"\t"<<p18<<"\t"<<p19<<"\t"<<p20<<"\t"<<p21<<endl;
    
//  }

  double pdf1[nVs], pdf2[nVs], pdf3[nVs];
  double chi2[nVs];
  double prediction = 0;

  int ndf = 0;
  for(int i=0; i<nVs; i++){
      pdf1[i]=0; pdf2[i]=0; pdf3[i]=0;
      chi2[i] =0;
      int nbins = (int)ds.ghist[i]->GetNbinsX();

      prediction = 0;
      for(int j=1; j<nbins-1; j++){
          double x = ds.ghist[i]->GetBinCenter(j+1);
          double y = ds.ghist[i]->GetBinContent(j+1);

          if(x>50) continue;

          pdf1[i] = TMath::Exp(-1.0*mu) * ( (1.0-w[i])*myGausPDF(x,q0[i],sigma0[i])
                                      +w[i]*myExpGausPDF(x,q0[i],sigma0[i],c0[i]) );
          // test nan value
          //if(TMath::IsNaN(q0[i]) || TMath::IsNaN(sigma0[i]) )
          //  cout<<"weried place 1 for nan value. "<<x<<"\t"<<q0[i] <<"\t"<<sigma0[i]<<"\t"<<c0[i]<<"\t"<<c0[i]<<endl;

          pdf2[i] = mu*TMath::Exp(-1.0*mu) * ( (1.0-w[i])*myGausPDF(x,q[i]+q0[i],TMath::Sqrt(sigma[i]*sigma[i]+sigma0[i]*sigma0[i]))
                           + w[i]*myExpGausPDF(x, q[i]+q0[i], TMath::Sqrt(sigma[i]*sigma[i]+sigma0[i]*sigma0[i]), c0[i]) );
          // test nan value
          //if( TMath::IsNaN(q[i]+q0[i]) || TMath::IsNaN(TMath::Sqrt(sigma[i]*sigma[i]+sigma0[i]*sigma0[i])) ){
          //   cout<<"wiried place 2 for nan: "<<x<<"\t"<<mu<<"\t"<<q[i]<<"\t"<<q0[i]<<"\t"
          //                                   <<sigma[i]<<"\t"<<sigma0[i]<<"\t"<<c0[i]<<"\t"<<w[i]<<endl;
          //}

          for(int k=2; k<=15; k++){
              pdf3[i] = TMath::Power(mu, k)*TMath::Exp(-1.0*mu)/TMath::Factorial(k)
                      * ( (1.0-w[i])*myGausPDF(x,q[i]*k+q0[i], TMath::Sqrt(sigma[i]*sigma[i]*k+sigma0[i]*sigma0[i]))
                       + w[i]*myExpGausPDF(x,q[i]*k+q0[i], TMath::Sqrt(sigma[i]*sigma[i]*k+sigma0[i]*sigma0[i]), c0[i]) ) ;

              // test nan value
             //if( TMath::IsNaN(q[i]*k+q0[i]) || TMath::IsNaN(TMath::Sqrt(sigma[i]*sigma[i]*k+sigma0[i]*sigma0[i])) ){
             //  cout<<"wiried place 3 for nan: "<<k<<"\t"<<x<<"\t"<<mu<<"\t"<<"\t"<<q[i]<<"\t"<<q0[i]<<"\t"<<sigma[i]<<"\t"<<sigma0[i]<<"\t"<<c0[i]<<"\t"<<w[i]<<endl;
             //}

          }

          prediction = a[i]*(pdf1[i]+pdf2[i]+pdf3[i]);

          // apply the poisson likelihood
          //if(prediction!=0){
            if(y!=0)
              chi2[i] += TMath::Power( y-prediction, 2)/y;
              ndf++;
            //else chi2[i] += 2.0*prediction;
            //cout<<"chi2 "<<i<<chi2[i]<<"\t"<<prediction<<"\t"<<y<<endl;
          //}
          //else{
            // just for test
            //cout<<"in fit p=0 : "<<"bin "<<j<<" , x= "<<x<<" , y= "<<y<<endl;
          //}
          //cout<<"thoise pdfs "<<i<<"\t"<<pdf1[i]<<"\t"<<pdf2[i]<<"\t"<<pdf3[i]<<endl;
      }
  }
  ds.ndf=ndf;
  //cout<<"return value "<<chi2[0]+chi2[1]+chi2[2]<<endl;
  return chi2[0]+chi2[1]+chi2[2];

}

/*
// --------------------------------------------------------------------
// this is the function that is used at this moment
// for computing the chi2, using the ideal response of a PMT
double CaliGainNovData::ResponseWithPedestal(double npe,
                                      double Q1, double w1, double a1,            
                                      double Q2, double w2, double a2,
                                      double Q3, double w3, double a3,
                                      double pedm1, double pedw1, double pedamp1,
                                      double pedm2, double pedw2, double pedamp2,
                                      double pedm3, double pedw3, double pedamp3)
{
  //cout<<"idealResponse is called ..."<<endl;
  const int npoints = (int) ds.ghist.size(); //3;//this->getNbOfDataPoints();
  //cout<<"check "<<global_hist.size()<<endl;
  if(npoints==0) {cout<<"fitting no data in CaliGainNovData:myfcn"<<endl; exit(0);}

  // I have these parameters:
  //double npe = mu; // mean npe
  double q[npoints]; q[0]=Q1; q[1]=Q2; q[2]=Q3;
  double sigma[npoints]; sigma[0]=w1; sigma[1]=w2; sigma[2]=w3;
  double amp[npoints]; amp[0]=a1; amp[1]=a2; amp[2]=a3;
  double q0[npoints]; q0[0]=pedm1; q0[1]=pedm2; q0[2]=pedm3;
  double w0[npoints]; w0[0]=pedw1; w0[1]=pedw2; w0[2]=pedw3;
  double a0[npoints]; a0[0]=pedamp1; a0[1]=pedamp2; a0[2]=pedamp3;

  // compute the chi2
  double chi2[3]; 
  for(int i=0; i<npoints; i++){
    chi2[i]=0.0;
    int nbins = (int)ds.ghist[i]->GetNbinsX();//200;

    for(int j=0; j<nbins-1; j++){
      double x = ds.ghist[i]->GetBinCenter(j+1);
      double y = ds.ghist[i]->GetBinContent(j+1);
      double prediction = 0;

      // ideal response term
      double p2 = 0;
      for(int k=1;k<100;k++){
          p2 +=    (TMath::Power(npe,k)*TMath::Exp(-1.0*npe)/TMath::Factorial(k)
            *TMath::Exp(-1.0*(x-q[i]*k)*(x-q[i]*k) 
              /(2.0*k*sigma[i]*sigma[i]))/(sigma[i]*TMath::Sqrt(2.0*TMath::Pi()*k)) ) ;
      }

      prediction = p2*amp[i] 
                         + a0[i]/TMath::Sqrt(2.0*TMath::Pi()*w0[i]*w0[i])
                           *TMath::Exp(-1.0*(x-q0[i])*(x-q0[i])/(2.0*w0[i]*w0[i])) ;

      if(y!=0)
        chi2[i] += TMath::Power( y-prediction , 2)/y; // the chi2 function
    }
  }

  return chi2[0]+chi2[1]+chi2[2];

}
*/

/*
// --------------------------------------------------------------------
// this is the function that is used at this moment
// for computing the chi2, using the response with pedestal, expoentially modified guassian
double CaliGainNovData::ResponseWithPedestal_expGaus(double npe,
                                      double Q1, double w1, double a1,            
                                      double Q2, double w2, double a2,
                                      double Q3, double w3, double a3,
                                      double pedm1, double pedw1, double pedamp1, double pedtau1,
                                      double pedm2, double pedw2, double pedamp2, double pedtau2,
                                      double pedm3, double pedw3, double pedamp3, double pedtau3)
{
  //cout<<"idealResponse is called ..."<<endl;
  const int npoints = (int) ds.ghist.size(); //3;//this->getNbOfDataPoints();
  //cout<<"check "<<global_hist.size()<<endl;
  if(npoints==0) {cout<<"fitting no data in CaliGainNovData:myfcn"<<endl; exit(0);}

  // I have these parameters:
  //double npe = mu; // mean npe
  double q[npoints]; q[0]=Q1; q[1]=Q2; q[2]=Q3;
  double sigma[npoints]; sigma[0]=w1; sigma[1]=w2; sigma[2]=w3;
  double amp[npoints]; amp[0]=a1; amp[1]=a2; amp[2]=a3;
  double q0[npoints]; q0[0]=pedm1; q0[1]=pedm2; q0[2]=pedm3; 
  double w0[npoints]; w0[0]=pedw1; w0[1]=pedw2; w0[2]=pedw3;
  double a0[npoints]; a0[0]=pedamp1; a0[1]=pedamp2; a0[2]=pedamp3;
  double tau0[npoints]; tau0[0]=pedtau1; tau0[1]=pedtau2; tau0[2]=pedtau3;

  // compute the chi2
  int ndf=0;
  double chi2[3]; 
  for(int i=0; i<npoints; i++){
    chi2[i]=0.0;
    int nbins = (int)ds.ghist[i]->GetNbinsX();//200;

    for(int j=0; j<nbins-1; j++){
      double x = ds.ghist[i]->GetBinCenter(j+1);
      double y = ds.ghist[i]->GetBinContent(j+1);
      double prediction = 0;

      // ideal response term
      double p2 = 0;
      for(int k=1;k<100;k++){
          p2 +=    (TMath::Power(npe,k)*TMath::Exp(-1.0*npe)/TMath::Factorial(k)
            *TMath::Exp(-1.0*(x-q[i]*k)*(x-q[i]*k) 
              /(2.0*k*sigma[i]*sigma[i]))/(sigma[i]*TMath::Sqrt(2.0*TMath::Pi()*k)) ) ;
      }

      prediction = p2*amp[i] 
                         + a0[i]*w0[i]/tau0[i]*TMath::Sqrt(TMath::Pi()/2.0)
                           *TMath::Exp(0.5*(w0[i]*w0[i])/(tau0[i]*tau0[i])-(x-q0[i])/tau0[i])
                           *TMath::Erfc(TMath::Sqrt(1./2.)*(w0[i]/tau0[i]-(x-q0[i])/w0[i]));

      if(y!=0){
        chi2[i] += TMath::Power( y-prediction , 2)/y; // the chi2 function
        ndf++;
      }
    }
  }

  ds.ndf=ndf;
  return chi2[0]+chi2[1]+chi2[2];

}
*/

// --------------------------------------------------------------------
// this is the real function that is used at this moment
// for computing the chi2, using the ideal response of a PMT
double CaliGainNovData::idealResponse(double npe,
                                      double Q1, double w1, double a1,            
                                      double Q2, double w2, double a2,
                                      double Q3, double w3, double a3)
{
  //cout<<"idealResponse is called ..."<<endl;
  const int npoints = (int) ds.ghist.size(); //3;//this->getNbOfDataPoints();
  //cout<<"check "<<global_hist.size()<<endl;
  if(npoints==0) {cout<<"fitting no data in CaliGainNovData:myfcn"<<endl; exit(0);}

  // I have these parameters:
  //double npe = mu; // mean npe
  double q[npoints]; q[0]=Q1; q[1]=Q2; q[2]=Q3;
  double sigma[npoints]; sigma[0]=w1; sigma[1]=w2; sigma[2]=w3;
  double amp[npoints]; amp[0]=a1; amp[1]=a2; amp[2]=a3;

  // compute the chi2
  double chi2 = 0;
  int ndf = 0;
  for(int i=0; i<npoints; i++){
    int nbins = (int)ds.ghist[i]->GetNbinsX();//200;
    for(int j=0; j<nbins-1; j++){
      double x = ds.ghist[i]->GetBinCenter(j+1);
      //if(x<0.25) continue;
      double y = ds.ghist[i]->GetBinContent(j+1);
      //fstream fouttest("testout.txt",ios::out|ios::app);
      //if(i==0) {fouttest<<j<<x<<"\t"<<y<<endl;}
      double prediction = 0;
      for(int k=1;k<100;k++){
          prediction += (TMath::Power(npe,k)*TMath::Exp(-1.0*npe)/TMath::Factorial(k)
            *TMath::Exp(-1.0*(x-q[i]*k)*(x-q[i]*k)
              /(2.0*k*sigma[i]*sigma[i]))/(sigma[i]*TMath::Sqrt(2.0*TMath::Pi()*k)));
      }      

      //if(prediction!=0){
        if(y!=0)
          chi2 += TMath::Power( y-prediction*amp[i] , 2)/y; // the chi2 formula Neyman
          ndf++;
          //chi2 += 2.0*(prediction*amp[i] - y + y*TMath::Log(y/prediction*amp[i]) );  // Poisson chi2
        //if(y==0) chi2 += 2.0*prediction*amp[i];
      //}
    }
  }

  ds.ndf = ndf;

  return chi2;

}


// --------------------------------------------------------------------------------
// this is the function to be called by the fitter
// inside it will call another function to do real calculations
void CaliGainNovData::myfcn(int &npar, double * deriv, double &f, double * par, int iflag){
  f = fitobj->idealResponse(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9]);
}

/*
// --------------------------------------------------------------------------------
// this is the function to be called by the fitter
// inside it will call another function to do real calculations
void CaliGainNovData::myfcnWithPedestal(int &npar, double * deriv, double &f, double * par, int iflag){
  f = fitobj->ResponseWithPedestal(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], par[10], par[11], par[12], par[13], par[14], par[15], par[16], par[17], par[18]);
}
*/

/*
// --------------------------------------------------------------------------------
// this is the function to be called by the fitter
// inside it will call another function to do real calculations
void CaliGainNovData::myfcnWithPedestal_expGaus(int &npar, double * deriv, double &f, double * par, int iflag){
  f = fitobj->ResponseWithPedestal_expGaus(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], par[10], par[11], par[12], par[13], par[14], par[15], par[16], par[17], par[18], par[19], par[20], par[21]);
}
*/

// --------------------------------------------------------------------------------
// this is the function to be called by the fitter
// inside it will call another function to do real calculations
void CaliGainNovData::mygainfcn(int &npar, double * deriv, double &f, double * par, int iflag){
  //cout<<"mygainfcn called with "<<npar<<" parameters : " <<endl;
  //for(int i=0; i<22;i++) cout<<par[i]<<"\t";
  //cout<<endl;
  f = fitobj->mygainResponse(par[0],
                             par[1],par[2],par[3],par[4],par[5],par[6],par[7],
                             par[8],par[9],par[10],par[11],par[12],par[13],par[14],
                             par[15],par[16],par[17],par[18],par[19],par[20],par[21]);
  //cout<<"f="<<f<<endl;
}

void CaliGainNovData::initializeParVectors(int npars){
  vresult.resize(npars);
  verrors.resize(npars);
  for(int i=0; i<npars; i++){vresult[i]=0.0; verrors[i]=0.0;}
  minChi2=0.0;      
  m_fitstatus = 0.0;
}
