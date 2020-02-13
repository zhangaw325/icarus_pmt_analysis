////////////////////////////////////////////////////////////////////////////////
// Waveform class definition
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef  __WAVEFORM_H
#define __WAVEFORM_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "TH1D.h"

using namespace std;

struct Pulse
{
    double start_time;
    double width;
    double amplitude;
    double integral;
};

class Waveform
{

    public:

      typedef vector<unsigned short> Rawdigits_t;
      typedef vector<double> Waveform_t;

      Waveform();
      Waveform( vector<unsigned short> rawWaveform );
      ~Waveform();

      bool hasSignal(double n_sigma);

      bool hasPulse(double n_sigma, Pulse & pulse);

      Waveform_t waveform;

      double baselineMean;
      double baselineWidth;

    private:

      int n_sample_baseline;

      void removeBaseline();

};


#endif //__WAVEFORM_H
