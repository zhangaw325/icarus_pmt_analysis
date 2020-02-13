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

class Waveform
{

  public:

    typedef vector<unsigned short> Rawdigits_t;
    typedef vector<double> Waveform_t;

    Waveform();
    Waveform( vector<unsigned short> waveform );
    ~Waveform();

    Waveform_t signal;

    double baselineMean;
    double baselineWidth;

  private:

    void removeBaseline();

    int n_sample_baseline;

};

#endif //__WAVEFORM_H
