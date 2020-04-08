#ifndef  __CaliGainNovData_H
#define __CaliGainNovData_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "TH1D.h"

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
    void plotHist();

    int getNbOfDataPoints(){return (int)voltages.size();}
    vector<TH1D*> getHistograms(){return chargehist; };

    void fitGainCurve();

    //void myfcn(int &npar, double * deriv, double &f, double * par, int iflag);
    static void myfcn(int &npar, double * deriv, double &f, double * par, int iflag);
    double idealResponse(double mu, double Q1, double w1, double a1, double Q2, double w2, double a2, double Q3, double w3, double a3);

  private:

    vector<TH1D*> chargehist;
    vector<double> voltages;
    int m_board;
    int m_ch;
};

#endif
