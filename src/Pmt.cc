////////////////////////////////////////////////////////////////////////////////
// PMT Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Pmt.h"

PMT::PMT(int run, int event, int board, int channel, int n_samples,
                                                 Waveform::Rawdigits_t waveform)
  : fRun(run)
  , fEvent(event)
  , fBoard(board)
  , fChannel(channel)
  , fNSamples(n_samples)
  , m_sampling_period(1.0)
{
  fRawWaveform=waveform;
  fWaveform = Waveform(waveform);
};

//------------------------------------------------------------------------------

PMT::~PMT(){};

//------------------------------------------------------------------------------

bool PMT::find(int run, int event, int board, int channel )
{
  return (run==fRun && event==fEvent && board==fBoard && channel==fChannel);
};

//------------------------------------------------------------------------------

void PMT::getBaselineParams( double &mean, double &stdev )
{
  mean = fWaveform.baselineMean;
  stdev = fWaveform.baselineWidth;
};

//------------------------------------------------------------------------------

TH1D* PMT::getWaveformHist()
{
  TH1D *hist = new TH1D("hist", "hist", fNSamples, 0, fNSamples*m_sampling_period);
  for(int t=0; t<fNSamples; t++){ hist->Fill( t, fWaveform.signal.at(t) ); }
  return hist;
};
