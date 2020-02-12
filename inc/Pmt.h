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

  typedef vector<unsigned short> Waveform_t;

  public:
    PMT(int run, int event, int board, int channel, int n_samples);
    ~PMT();

    // Setters
    void setWaveform( PMT::Waveform_t waveform  ){ fWaveform = waveform; }

    // Getters
    int getRun(){ return fRun; };
    int getEvent(){ return fEvent; };
    int getBoard(){ return fBoard; };
    int getChannel(){ return fChannel; }

    TH1D* getWaveformHist();

    // Helper to identify a PMT in list
    bool find(int run, int event, int board, int channel ){
      return (run==fRun && event==fEvent && board==fBoard && channel==fChannel);
    };

  private:
    int fRun;
    int fEvent;
    int fBoard;
    int fChannel;
    int fNSamples;

    double m_sampling_period;

    PMT::Waveform_t fWaveform;

};

#endif //__PMT_H
