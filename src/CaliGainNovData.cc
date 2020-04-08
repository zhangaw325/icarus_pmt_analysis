#include "CaliGainNovData.h"
#include "IdealChargeFunction.h"
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

//using namespace std;

//global variables for hosting charge distribution
struct ds_global {
  vector<TH1D*> ghist;
  vector<double> gv;
} ds;

//ds.ghist.resize(0);
//ds.gv.resize(0);

static CaliGainNovData* fitobj;

CaliGainNovData::CaliGainNovData(){
  fitobj=this;
  this->chargehist.resize(0);
  this->voltages.resize(0);
}

CaliGainNovData::CaliGainNovData(int board, int channel):m_board(board),m_ch(channel)
{
      chargehist.resize(0);
      voltages.resize(0);
}

CaliGainNovData::~CaliGainNovData(){
  for(int i=0; i<(int)chargehist.size(); i++)
    delete chargehist[i]; 
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
  for(int i=0; i<npoints; i++){
    int nbins = (int)ds.ghist[i]->GetNbinsX();//200;

    for(int j=0; j<nbins-1; j++){
      double x = ds.ghist[i]->GetBinCenter(j+1);
      double y = ds.ghist[i]->GetBinContent(j+1);
      //fstream fouttest("testout.txt",ios::out|ios::app);
      //if(i==0) {fouttest<<j<<x<<"\t"<<y<<endl;}
      double prediction = 0;
      for(int k=1;k<100;k++){
          prediction += (TMath::Power(npe,k)*TMath::Exp(-1.0*npe)/TMath::Factorial(k)
            *TMath::Exp(-1.0*(x-q[i]*k)*(x-q[i]*k)
              /(2.0*k*sigma[i]*sigma[i]))/(sigma[i]*TMath::Sqrt(2.0*TMath::Pi()*k)));
      }      

      if(y!=0)
        chi2 += TMath::Power( y-prediction*amp[i] , 2)/y; // the chi2 function
    }
  }

  return chi2;

}

// --------------------------------------------------------------------------------
// this is the function to be called by the fitter
// inside it will call another function to do real calculations
void CaliGainNovData::myfcn(int &npar, double * deriv, double &f, double * par, int iflag){
  f = fitobj->idealResponse(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9]);
}

// -------------------------------------------------------------------------------
// the hook for starting the global fit to charge distributions

void CaliGainNovData::fitGainCurve(TFile* ofile){
    if(chargehist.size()==0){ cout<<"fitGainCurve will not work if no data points"<<endl; exit(0); }

    // set the global variables to be used by the fitter
    for(int i=0; i< (int)chargehist.size(); i++){
        ds.ghist.push_back(chargehist[i]);
        ds.gv.push_back(voltages[i]);
    }

    // define number of parameters and their initial values
    int n = getNbOfDataPoints() ;
    const int n_pars = 10;//1 + 3*n ;
    double vstart[n_pars]={20, 
                       1.0, 0.5, chargehist[0]->Integral()/2.0,
                       1.0, 0.5, chargehist[1]->Integral()/2.0,
                       1.0, 0.5, chargehist[2]->Integral()/2.0} ;
    double step[n_pars]={0.5, 0.1, 0.1, 20, 0.1, 0.1, 20, 0.1, 0.1, 20};
    double minVal[n_pars]={0, 0.01, 0.001, 0, 0.01, 0.001, 0, 0.01, 0.001, 0};
    double maxVal[n_pars]={50, 10, 10,chargehist[0]->Integral()*2.0, 
                          10, 10,chargehist[1]->Integral()*2.0, 
                          10, 10,chargehist[2]->Integral()*2.0};
    string parName[n_pars]={"npe",
                            "q1", "w1", "a1",
                            "q2", "w2", "a2",
                            "q3", "w3", "a3"
                           };
    // set the fitter: use TMinuit
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
    double outpar[n_pars], errors[n_pars];
    for(int i=0; i<n_pars; i++){
      minimizer->GetParameter(i, outpar[i], errors[i]);
    }
    cout<<"check fit status "<<fitstatus<<endl;
    if(fitstatus==4){
      for(int i=0; i<n_pars; i++){
        minimizer->DefineParameter(i, parName[i].c_str(), outpar[i], step[i]/2.0, 
                               minVal[i], maxVal[i]);
      }
      fitstatus = minimizer->Migrad();
      cout<<"status of the second fit "<<fitstatus<<endl;
      for(int i=0; i<n_pars; i++){
        minimizer->GetParameter(i, outpar[i], errors[i]);
      }
    }

    // here we get the result from the fit
    double amin,edm,errdef;
    int nvpar,nparx,icstat;
    minimizer->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    cout<<amin<<"\t"<<edm<<"\t"<<errdef<<"\t"<<nvpar<<"\t"<<nparx<<"\t"<<icstat<<endl;

    // plot gain curves and save them in the root file
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
        gGain->SetPoint(i, voltages[i], outpar[voltages.size()*i+1]/1.6);
        gGain->SetPointError(i, 0, errors[voltages.size()*i+1]/1.6);
    }
    gGain->Draw("ALP");
    gGain->SetMarkerStyle(8);
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
      if(i==0) chargehist[i]->Draw();
      else chargehist[i]->Draw("same");
      chargehist[i]->SetLineColor(i+1);
      sprintf(gname,"idealfunc%d",i);
      func[i] = new TF1(gname,idealchargePlot_obj,0, 100, 4);
      func[i]->SetParNames("mu","q0","sigma","amp.");
      func[i]->SetParameters(outpar[0],
                          outpar[i*chargehist.size()+1],
                          outpar[i*chargehist.size()+2],
                          outpar[i*chargehist.size()+3]);
      func[i]->SetLineColor(1); func[i]->SetLineStyle(2);
      func[i]->Draw("same");
    }
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->SetHeader("HV values");
    for(int i=0; i<(int)chargehist.size(); i++){
      sprintf(gname, "%.1f", voltages[i]);
      leg->AddEntry(chargehist[i],gname,"lp");
    }
    leg->Draw(); 
    TPaveText* pt = new TPaveText(0.6,0.3,0.9,0.7, "brNDC");   
    pt->AddText("Fit parameters:");
    sprintf(gname,"npe: %.2f #pm %.2f", outpar[0], errors[0]);
    pt->AddText(gname);
    for(int i=0; i<(int)chargehist.size();i++){
      sprintf(gname, "q%d_%.2f/w%d_%.2f/a%d_%.1f",
                 i+1, outpar[i*chargehist.size()+1], 
                 i+1, outpar[i*chargehist.size()+2],
                 i+1, outpar[i*chargehist.size()+3]);
      pt->AddText(gname);

    }
    pt->Draw();
    pt->SetBorderSize(1);
    pt->SetFillColor(0);
    pt->SetTextFont(42);
    //pt->SetTextAlign(1);
    cCharge->Modified(); cCharge->Update();
    cCharge->Write();
    cCharge->Close();
    // in the end of the fit function, need to reset the global variables
    ds.ghist.resize(0);
    ds.gv.resize(0);
}
