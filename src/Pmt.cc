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

TH1D* PMT::getWaveformHist()
{
  TH1D *hist = new TH1D("hist", "hist", fNSamples, 0, fNSamples*m_sampling_period);
  for(int t=0; t<fNSamples; t++){ hist->Fill( t, double(fWaveform.at(t)) ); }
  return hist;
};
