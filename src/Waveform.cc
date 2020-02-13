////////////////////////////////////////////////////////////////////////////////
// Waveform methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"

//------------------------------------------------------------------------------

Waveform::Waveform()
  : n_sample_baseline(300)
  , baselineMean(0.0)
  , baselineWidth(0.0)
{ };

//------------------------------------------------------------------------------

Waveform::Waveform( Rawdigits_t waveform )
  : Waveform()
{
  for( auto w : waveform ){ signal.push_back( double(w) ); }
  removeBaseline();
};

//------------------------------------------------------------------------------

Waveform::~Waveform(){};

//------------------------------------------------------------------------------

void Waveform::removeBaseline()
{

  int n_samples=signal.size();

  // Calculate the baseline as the mean values on the last part of the spectrum
  for(int t=n_samples-1; t>n_samples-1-n_sample_baseline; t--)
  {
    baselineMean += signal.at(t);
  }
  baselineMean /= n_sample_baseline;

  // Calculate the stdev of the baseline
  for(int t=n_samples-1; t>n_samples-1-n_sample_baseline; t--)
  {
    baselineWidth += pow(signal.at(t)-baselineMean, 2);
  }
  baselineWidth = sqrt( baselineWidth / (n_sample_baseline-1) );

  // Subtract the baseline from the waveform
  std::transform(signal.begin(), signal.end(), signal.begin(),
                                   [&] (double x) { return x - baselineMean; });

};
