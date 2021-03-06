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
    void fitHighCharge_improve();
    //void fitLowCharge();
    void fitLowCharge_fullExpr();
    void fitLowCharge_fullExpr_Single();
    void fitLowCharge_expGaus();
    //void fitChargeDistributions();
    //void plotGainCurve(TFile* ofile, bool flag);
    void plotGainCurveHighCharge(TFile* ofile, bool flag);
    void plotGainCurveHighCharge_improve(TFile* ofile, bool flag);
    void plotGainCurveLowCharge(TFile* ofile, bool flag);
    void plotGainCurveLowCharge_expGaus(TFile* ofile, bool flag);
    void plotGainCurveLowCharge_fullExpr_Single(TFile* ofile, bool flag);
    void initializeParVectors(int npars);
    void SaveJointFitParameters(bool flag);

    //void myfcn(int &npar, double * deriv, double &f, double * par, int iflag);
    static void myfcn(int &npar, double * deriv, double &f, double * par, int iflag);
    double idealResponse(double mu, double Q1, double w1, double a1, double Q2, double w2, double a2, double Q3, double w3, double a3);

    static void myfcn_improve(int &npar, double * deriv, double &f, double * par, int iflag);
    double idealResponse_improve(double mu, double Q1, double w1, double a1, double Q2, double w2, double a2, double Q3, double w3, double a3);

    //static void myfcnWithPedestal(int &npar, double * deriv, double &f, double * par, int iflag);
    //double ResponseWithPedestal(double npe,
    //                                  double Q1, double w1, double a1,            
    //                                  double Q2, double w2, double a2,
    //                                  double Q3, double w3, double a3,
    //                                  double pedm1, double pedw1, double pedamp1,
    //                                  double pedm2, double pedw2, double pedamp2,
    //                                  double pedm3, double pedw3, double pedamp3 );

    static void myfcnWithPedestal_expGaus(int &npar, double * deriv, double &f, double * par, int iflag);
    double ResponseWithPedestal_expGaus(double npe,
                                      double Q1, double w1, double a1,            
                                      double Q2, double w2, double a2,
                                      double Q3, double w3, double a3,
                                      double pedm1, double pedw1, double pedamp1, double pedtau1,
                                      double pedm2, double pedw2, double pedamp2, double pedtau2,
                                      double pedm3, double pedw3, double pedamp3, double pedtau3 );

    void fitLowCharge_expGaus_improve();
    void plotGainCurveLowCharge_expGaus_improve(TFile* ofile, bool flag);
    double ResponseWithPedestal_expGaus_improve(double npe,
                                      double Q1, double w1, double a1,            
                                      double Q2, double w2, double a2,
                                      double Q3, double w3, double a3,
                                      double pedm1, double pedw1, double pedamp1, double pedtau1,
                                      double pedm2, double pedw2, double pedamp2, double pedtau2,
                                      double pedm3, double pedw3, double pedamp3, double pedtau3 );
    static void myfcnWithPedestal_expGaus_improve(int &npar, double * deriv, double &f, double * par, int iflag);

    // these two are my final charge distribution fitting function.
    // this is based on 3 voltage data points;
    static void mygainfcn(int &npar, double * deriv, double &f, double * par, int iflag);
    double mygainResponse(double p0,
                          double p1, double p2, double p3, double p4, double p5, double p6, double p7,
                          double p8, double p9, double p10, double p11, double p12, double p13, double p14,
                          double p15, double p16, double p17, double p18, double p19, double p20, double p21);

    double myGausPDF(double x, double mean, double sigma); // Gaus distribution
    double myExpGausPDF(double x, double mean, double sigma, double c0);//exp-mod Gauss

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
