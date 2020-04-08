#include "CaliGainNovData.h"
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

#include "TMinuit.h"
#include "TH1D.h"

//using namespace std;

double global_hist_y[3][200];
double global_hist_x[3][200];
//global_hist.resize(0);
/*
for(int i=0; i<3; i++){
  for(int j=0;j<200;j++){
    global_hist_y[i][j]=0;
    global_hist_x[i][j]=0;
  }
}
*/

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
  //for(int i=0; i<(int)chargehist.size(); i++)
  //  delete chargehist[i]; 
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

double CaliGainNovData::idealResponse(double npe,
                                      double Q1, double w1, double a1,            
                                      double Q2, double w2, double a2,
                                      double Q3, double w3, double a3)
{
  //cout<<"idealResponse is called ..."<<endl;
  const int npoints = 3;//this->getNbOfDataPoints();
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
    int nbins = 200;//global_hist[i]->GetNbinsX();

    for(int j=1; j<nbins-1; j++){
      double x = global_hist_x[i][j];//->GetBinCenter(j+1);
      double y = global_hist_y[i][j];//->GetBinContent(j+1);
      //fstream fouttest("testout.txt",ios::out|ios::app);
      //if(i==0) {fouttest<<j<<x<<"\t"<<y<<endl;}
      double prediction = 0;
      for(int k=1;k<75;k++){
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

//I want to make a global FNC function
void CaliGainNovData::myfcn(int &npar, double * deriv, double &f, double * par, int iflag){
  //cout<<"myfcn is called"<<endl;
  //for(int i=0; i<10; i++){
  //   cout<<"par"<<i+1<<"\t"<<par[i]<<endl;
  //}
  f = fitobj->idealResponse(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9]);
}

void CaliGainNovData::fitGainCurve(){
    if(chargehist.size()==0){ cout<<"fitGainCurve will not work if no data points"<<endl; exit(0); }

    cout<<"fitGainCurve called " << voltages[0]<<"\t"
                                <<voltages[1]<<"\t"<<voltages[2]<<endl;

    //cout<<"test output "<<chargehist[0]->Integral()<<"\t"
    //                    <<chargehist[1]->Integral()<<"\t"
    //                    <<chargehist[2]->Integral()<<"\n";

    //global_hist.resize(0);
    //for(int i=0; i< (int)chargehist.size(); i++){
    //    global_hist.push_back(chargehist[i]);
    //}

    for(int i=0; i<3; i++){
      for(int j=0;j<200;j++){
        global_hist_y[i][j]=chargehist[i]->GetBinContent(j+1);
        global_hist_x[i][j]=chargehist[i]->GetBinCenter(j+1);
      }
    }

    //fitobj->Initialize(chargehist, voltages);

    //int n = getNbOfDataPoints() ;
    const int n_pars = 10;//1 + 3*n ;

    double vstart[10]={20, 
                       1.0, 0.5, chargehist[0]->Integral()/2.0,
                       1.0, 0.5, chargehist[1]->Integral()/2.0,
                       1.0, 0.5, chargehist[2]->Integral()/2.0} ;
    double step[10]={0.5, 0.1, 0.1, 20, 0.1, 0.1, 20, 0.1, 0.1, 20};
    double minVal[10]={0, 0.01, 0.001, 0, 0.01, 0.001, 0, 0.01, 0.001, 0};
    double maxVal[10]={50, 10, 10,chargehist[0]->Integral()*2.0, 
                          10, 10,chargehist[1]->Integral()*2.0, 
                          10, 10,chargehist[2]->Integral()*2.0};
    string parName[10]={"npe",
                            "q1", "w1", "a1",
                            "q2", "w2", "a2",
                            "q3", "w3", "a3"
                           };

    //TFitter* minimizer = new TFitter(10);
    TMinuit* minimizer = new TMinuit(n_pars);
    minimizer->SetFCN(CaliGainNovData::myfcn);
    for(int i=0; i<n_pars; i++){
      minimizer->DefineParameter(i, parName[i].c_str(), vstart[i], step[i], 
                             minVal[i], maxVal[i]);
    }
    //minimizer->ExecuteCommand("SIMPLEX",0,0);
    //minimizer->ExecuteCommand("MIGRAD", vstart ,n_pars);
    minimizer->Migrad();

    double amin,edm,errdef;
    int nvpar,nparx,icstat;
    minimizer->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    cout<<amin<<"\t"<<edm<<"\t"<<errdef<<"\t"<<nvpar<<"\t"<<nparx<<"\t"<<icstat<<endl;
}