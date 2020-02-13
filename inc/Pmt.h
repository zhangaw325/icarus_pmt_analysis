////////////////////////////////////////////////////////////////////////////////
// PMT Class definition
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef  __PMT_H
#define __PMT_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "Waveform.h"

#include "TH1D.h"

using namespace std;

class PMT
{

  public:

    PMT(int run, int event, int board, int channel, int n_samples,
                                                Waveform::Rawdigits_t waveform);
    ~PMT();

    // Getters
    int getRun(){ return fRun; };
    int getEvent(){ return fEvent; };
    int getBoard(){ return fBoard; };
    int getChannel(){ return fChannel; }
    Waveform::Rawdigits_t getRawWaveform(){ return fRawWaveform; };
    Waveform::Waveform_t getWaveform(){ return fWaveform.signal; };
    void getBaselineParams( double &mean, double &stdev );

    TH1D* getWaveformHist();

    // Helper to identify a PMT in list
    bool find(int run, int event, int board, int channel );

  private:

    int fRun;
    int fEvent;
    int fBoard;
    int fChannel;
    int fNSamples;

    double m_sampling_period;

    vector<unsigned short> fRawWaveform;
    Waveform fWaveform;

};

#endif //__PMT_H
