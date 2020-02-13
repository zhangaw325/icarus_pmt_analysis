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

    PMT(int run, int subrun, int event, int board, int channel, int n_samples,
                                                Waveform::Rawdigits_t waveform);
    ~PMT();

    // Getters
    int getRun(){ return m_run; };
    int getEvent(){ return m_event; };
    int getBoard(){ return m_board; };
    int getChannel(){ return m_channel; }
    Waveform::Rawdigits_t getRawWaveform(){ return fRawWaveform; };
    Waveform::Waveform_t getWaveform(){ return fWaveform.waveform; };
    void getBaselineParams( double &mean, double &stdev );

    // Return quantites related to the pulse
    double getAmplitude(){ return pulse.amplitude; };
    double getIntegral(){ return pulse.integral; };
    double getWidth(){ return pulse.width; }
    double getStartTime(){ return pulse.start_time; }

    TH1D* getWaveformHist();

    // Helper to identify a PMT in list
    bool find(int run, int event, int board, int channel );

    // Helper that returns true if the signal on a given PMT is above a set threshold
    bool hasSignal( double threshold );
    bool hasPulse( double threshold );

  private:

    int m_run;
    int m_subrun;
    int m_event;
    int m_board;
    int m_channel;
    int m_nsamples;

    double m_sampling_period;

    vector<unsigned short> fRawWaveform;
    Waveform fWaveform;

    Pulse pulse;

};

#endif //__PMT_H
