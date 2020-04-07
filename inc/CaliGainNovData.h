#ifndef  __CaliGainNovData_H
#define __CaliGainNovData_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "TH1D.h"

using namespace std;

class CaliGainNovData
{

  public:

    CaliGainNovData();
    CaliGainNovData(int board, int channel);

    ~CaliGainNovData();

    void loadOneHist(TH1D* hist, double hv);
    void plotHist();

  private:

    vector<TH1D*> chargehist;
    vector<double> voltages;
    int m_board;
    int m_ch;
};

#endif
