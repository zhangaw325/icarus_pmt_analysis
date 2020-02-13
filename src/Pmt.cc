////////////////////////////////////////////////////////////////////////////////
// PMT Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Pmt.h"

PMT::PMT(int run, int event, int board, int channel, int n_samples)
  : fRun(run)
  , fEvent(event)
  , fBoard(board)
  , fChannel(channel)
  , fNSamples(n_samples)
  , m_sampling_period(1.0)
{
  // init
};

//------------------------------------------------------------------------------

PMT::~PMT(){};

//------------------------------------------------------------------------------

void PMT::setWaveform( vector<unsigned short> waveform  )
{

  // fRawWaveform store the original waveform
  fRawWaveform = waveform;

  // Cast type to double for further tests
  for( auto w : waveform ){ fWaveform.push_back( double(w) ); }

  // Remove baseline automatically
  this->removeBaseline();

};

//------------------------------------------------------------------------------

bool PMT::find(int run, int event, int board, int channel )
{
  return (run==fRun && event==fEvent && board==fBoard && channel==fChannel);
};

//------------------------------------------------------------------------------

void PMT::removeBaseline( int n_sample_baseline )
{

  // Calculate the baseline as the mean values on the last part of the spectrum
  for(int t=fNSamples-1; t>fNSamples-1-n_sample_baseline; t--)
  {
    fBaseline += fWaveform.at(t);
  }
  fBaseline /= n_sample_baseline;

  // Calculate the stdev of the baseline
  for(int t=fNSamples-1; t>fNSamples-1-n_sample_baseline; t--)
  {
    fBaselineWidth += pow(fWaveform.at(t)-fBaseline, 2);
  }
  fBaselineWidth = sqrt( fBaselineWidth / (n_sample_baseline-1) );

  // Subtract the baseline from the waveform
  std::transform(fWaveform.begin(), fWaveform.end(), fWaveform.begin(),
                                      [&] (double x) { return x - fBaseline; });

};

//------------------------------------------------------------------------------

void PMT::getBaselineParams( double &mean, double &stdev )
{
  mean = fBaseline;
  stdev = fBaselineWidth;
};

//------------------------------------------------------------------------------

TH1D* PMT::getWaveformHist()
{
  TH1D *hist = new TH1D("hist", "hist", fNSamples, 0, fNSamples*m_sampling_period);
  for(int t=0; t<fNSamples; t++){ hist->Fill( t, double(fWaveform.at(t)) ); }
  return hist;
};
