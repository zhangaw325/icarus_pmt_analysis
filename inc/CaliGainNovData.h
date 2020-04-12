#ifndef  __CaliGainNovData_H
#define __CaliGainNovData_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "TH1D.h"
#include "TFile.h"

#include "TFitter.h"

using namespace std;

class CaliGainNovData
{

  public:

    CaliGainNovData();
    CaliGainNovData(int board, int channel);

    ~CaliGainNovData();
    void Initialize(vector<TH1D*>  hist, vector<double> v);

    void loadOneHist(TH1D* hist, double hv);
    void loadOneLowChargeHist(TH1D* hist);
    void plotHist();

    int getNbOfDataPoints(){return (int)voltages.size();}
    vector<TH1D*> getHistograms(){return chargehist; };
    vector<TH1D*> getLowChargeHist(){return lowchargehist;}

    void fitGainCurve(TFile* ofile);
    void fitHighCharge();
    void fitLowCharge();
    void plotGainCurve(TFile* ofile, bool flag);
    void initializeParVectors(const int npars);

    //void myfcn(int &npar, double * deriv, double &f, double * par, int iflag);
    static void myfcn(int &npar, double * deriv, double &f, double * par, int iflag);
    double idealResponse(double mu, double Q1, double w1, double a1, double Q2, double w2, double a2, double Q3, double w3, double a3);

    static void myfcnWithPedestal(int &npar, double * deriv, double &f, double * par, int iflag);
    double ResponseWithPedestal(double npe,
                                      double Q1, double w1, double a1,            
                                      double Q2, double w2, double a2,
                                      double Q3, double w3, double a3,
                                      double pedm1, double pedw1, double pedamp1,
                                      double pedm2, double pedw2, double pedamp2,
                                      double pedm3, double pedw3, double pedamp3 );

  private:

    vector<TH1D*> chargehist;
    vector<TH1D*> lowchargehist;
    vector<double> voltages;
    int m_board;
    int m_ch;
    
    // vectors to store fit parameters
/*
    vector<double> vstart;
    vector<double> vstep;
    vector<double> minVal;
    vector<double> maxVal;
    vector<string> parName;
*/
    vector<double> vresult;
    vector<double> verrors;
    double minChi2;
    int m_fitstatus;
};

#endif
