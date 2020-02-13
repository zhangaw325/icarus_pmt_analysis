////////////////////////////////////////////////////////////////////////////////
// PMT Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Pmt.h"

PMT::PMT(int run, int subrun, int event, int board, int channel, int n_samples,
                                                 Waveform::Rawdigits_t waveform)
  : m_run(run)
  , m_subrun(subrun)
  , m_event(event)
  , m_board(board)
  , m_channel(channel)
  , m_nsamples(n_samples)
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
  return (run==m_run && event==m_event && board==m_board && channel==m_channel);
};

//------------------------------------------------------------------------------

bool PMT::hasSignal( double threshold )
{
  // Just checking for the pmts having a signal
  return fWaveform.hasSignal(threshold);
};

//------------------------------------------------------------------------------

bool PMT::hasPulse( double threshold )
{
  // If a pulse is found, then it fills also the object pulse
  return fWaveform.hasPulse(threshold, pulse);
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
  TH1D *hist = new TH1D("hist", "hist", m_nsamples, 0, m_nsamples*m_sampling_period);
  for(int t=0; t<m_nsamples; t++){ hist->Fill( t, fWaveform.waveform.at(t) ); }
  return hist;
};
