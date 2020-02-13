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

#include "TH1D.h"

using namespace std;

class PMT
{

  public:

    typedef vector<double> Waveform_t;

    PMT(int run, int event, int board, int channel, int n_samples);
    ~PMT();

    // Setters
    void setWaveform( vector<unsigned short> waveform  );

    // Getters
    int getRun(){ return fRun; };
    int getEvent(){ return fEvent; };
    int getBoard(){ return fBoard; };
    int getChannel(){ return fChannel; }
    vector<unsigned short> getRawWaveform(){ return fRawWaveform; };
    Waveform_t getWaveform(){ return fWaveform; };
    void getBaselineParams( double &mean, double &stdev );

    TH1D* getWaveformHist();

    // Basic signal processing functions
    void removeBaseline( int n_sample_baseline=300 );

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
    PMT::Waveform_t fWaveform;

    double fBaseline=0.0; // Mean baseline value
    double fBaselineWidth=0.0; // Stdev baseline

};

#endif //__PMT_H
