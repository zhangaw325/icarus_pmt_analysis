////////////////////////////////////////////////////////////////////////////////
// Event methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Event.h"

Event::Event(){ };

//------------------------------------------------------------------------------

Event::Event(int run, int subrun, int event, int nevents)
  : m_run(run)
  , m_subrun(subrun)
  , m_event(event)
  , m_nevents(nevents)
{ };

//------------------------------------------------------------------------------

Event::Event(int run, int subrun, int event, int nevents, Pmts_t pmt_vector)
  : m_run(run)
  , m_subrun(subrun)
  , m_event(event)
  , m_nevents(nevents)
{
  this->setPmtVector( pmt_vector );
};

//------------------------------------------------------------------------------

Event::~Event()
{
  this->clean();
};

//------------------------------------------------------------------------------

Event::Pmts_t Event::getValidPmts(double n_sigma)
{
  // Get the PMTs with a signal n_sigma above the baseline

  Event::Pmts_t valid_pmts;
  for(auto pmt : m_pmt_vector)
  {
    if( pmt->hasSignal(n_sigma) ){ valid_pmts.push_back(pmt); }
  }

  return valid_pmts;
};

//------------------------------------------------------------------------------

Event::Pmts_t Event::getPmtsBoard(int board)
{
  // Get the PMTs with a signal n_sigma above the baseline

  Event::Pmts_t valid_pmts;
  for(auto pmt : m_pmt_vector)
  {
    if( pmt->getBoard() == board ){ valid_pmts.push_back(pmt); }
  }

  return valid_pmts;
};

//------------------------------------------------------------------------------

Event::Pmts_t Event::getPmtsChannel(int channel)
{
  // Get the PMTs with a signal n_sigma above the baseline

  Event::Pmts_t valid_pmts;
  for(auto pmt : m_pmt_vector)
  {
    if( pmt->getChannel() == channel ){ valid_pmts.push_back(pmt); }
  }

  return valid_pmts;
};

//------------------------------------------------------------------------------

/*
void Event::printWaveform(int run, int event, vector<PMT*> pmts )
{

    for( auto pmt : pmts )
    {
      int run = this->m_run;
      int event = this->m_event;
      int channel=pmt->getChannel();
      int board=pmt->getBoard();

      auto it = find_if( pmts.begin(), pmts.end(), [ & ] ( PMT *pmt ){
                                return pmt->find(run, event, board, channel); } );

      if( it != pmts.end() )
      {

        printf("Plot RUN %d Event %d Board %d Channel %d \n", run, event, board, channel);

        PMT *pmt = *it;

        cout << pmt->getAmplitude() << " " << pmt->getIntegral() << " " << pmt->getStartTime() << " " << pmt->getWidth() << endl;

        char name[100];
        sprintf(name, "Run%d-Event%d-Board%d-Channel%d", run, event, board, channel);
        TCanvas *c = new TCanvas(name, name, 600, 400);

        TH1D *hist = pmt->getWaveformHist();

        hist->SetName(name);
        hist->SetTitle(name);
        hist->Draw("hist");

      }
      else
      {
        cout << "Impossible to find PMTs  with board " << board
                                          << " and channel " << channel << endl;
        return;
      }
    }
}
*/
