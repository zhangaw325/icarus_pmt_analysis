////////////////////////////////////////////////////////////////////////////////
// Waveform methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"

//------------------------------------------------------------------------------

Waveform::Waveform()
  : n_sample_baseline(200)
  , baselineMean(0.0)
  , baselineWidth(0.0)
{ };

//------------------------------------------------------------------------------

Waveform::Waveform( Rawdigits_t rawWaveform )
  : Waveform()
{
  for( auto w : rawWaveform ){ waveform.push_back( double(w) ); }
  removeBaseline();
};

//------------------------------------------------------------------------------

Waveform::~Waveform(){};

//------------------------------------------------------------------------------

void Waveform::removeBaseline()
{

  int n_samples=waveform.size();

  // Calculate the baseline as the mean values on the last part of the spectrum
  for(int t=n_samples-1; t>n_samples-1-n_sample_baseline; t--)
  {
    baselineMean += waveform.at(t);
  }
  baselineMean /= n_sample_baseline;

  // Calculate the stdev of the baseline
  for(int t=n_samples-1; t>n_samples-1-n_sample_baseline; t--)
  {
    baselineWidth += pow(waveform.at(t)-baselineMean, 2);
  }
  baselineWidth = sqrt( baselineWidth / (n_sample_baseline-1) );

  // Subtract the baseline from the signal waveform
  std::transform(waveform.begin(), waveform.end(), waveform.begin(),
                                 [ & ] (double x) { return x - baselineMean; });

};

//------------------------------------------------------------------------------

bool Waveform::hasSignal(double n_sigma)
{
  // Find if the waveform has n consecutive counts above a threshold expressed
  // in number of sigmas

  bool has_signal=false;
  int n_counts=5.0;
  int counts=0;

  for( double value : waveform )
  {
    if( abs(value) > n_sigma*baselineWidth ){ counts++; }
    else{ counts=0; } // Reset the counts

    if( counts > n_counts ){ has_signal=true; break; }
  }

  return has_signal;
};

//------------------------------------------------------------------------------

bool Waveform::hasPulse(double n_sigma, Pulse & pulse)
{
  // Define the Pulse region of the signal as the region with >n consecutive
  // counts above a threshold expressed in number of sigmas.

  bool has_pulse=false;
  int counts=0;
  int n_counts=5.0;

  int t_start=0, t_end=0;

  for(int t=0; t<int(waveform.size()); t++)
  {
    double value = waveform.at(t);

    if( abs(value) > n_sigma*baselineWidth )
    {
      counts++;
      if( counts==1 ){ t_start=t-2; }
    }
    else
    {
      // if it is a pulse region, then this is the ending clause
      if( counts>n_counts && t_start>0 )
      {
        t_end=t+3;
        has_pulse=true;
        break; // Force to have only one pulse
      }
      else
      {
        // reset everything
        counts=0; t_start=0; t_end=0;
      }
    }
  }

  // If no pulse is found, we end the games here..
  if( !has_pulse ) { return has_pulse; }

  // Here we define the characteristics of the pulse
  double amp = 0., charge=0;
  for( int t=t_start; t<t_end; t++ )
  {
    if( waveform.at(t) < amp ){ amp = waveform.at(t); }
    charge += abs(waveform.at(t));
  }

  pulse.start_time = t_start;
  pulse.width = (t_end-t_start);
  pulse.amplitude = amp;
  pulse.integral = charge;

  return has_pulse;
};

//------------------------------------------------------------------------------
