////////////////////////////////////////////////////////////////////////////////
// PMT Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////
#include "Utils.h"
#include "Pmt.h"

PMT::PMT( )
{ };

//------------------------------------------------------------------------------

PMT::PMT( int run, int board, int channel )
  : m_run(run)
  , m_board(board)
  , m_channel(channel)
{
  m_pmt_number=m_channel*geo::nchannels*m_board;
  this->initHist();
};

//------------------------------------------------------------------------------

PMT::PMT( int run, int pmt_number)
  : m_run(run)
  , m_pmt_number(pmt_number)
{
  m_board=0; m_channel=0;
  this->getBoardAndChannel(m_pmt_number, m_board, m_channel);
  this->initHist();
};

//------------------------------------------------------------------------------

PMT::~PMT()
{};

//------------------------------------------------------------------------------

void PMT::initHist()
{
  char hname[100];

  // Amplitude
  sprintf(hname, "hist_board%d_channel%d_amplitude", m_board,  m_channel);
  h_amplitude = new TH1D(hname, hname, 100, 0, 2000);

  sprintf(hname, "hist_board%d_channel%d_amplitude_low", m_board,  m_channel);
  h_amplitude_low = new TH1D(hname, hname, 100, 0, 150);
};

//------------------------------------------------------------------------------

void PMT::getBoardAndChannel( int pmt_number, int & board, int & channel )
{
  //Return the board and channel number given the absolute pmt number
  board = pmt_number/geo::nchannels;
  channel = pmt_number%geo::nchannels;
};

//------------------------------------------------------------------------------

void PMT::loadWaveform( Waveform *waveform )
{

  // Fetch the relevant quantities from the waveform
  double amplitude = abs( waveform->getAmplitude() );
  double integral = abs( waveform->getIntegral() );

  // Fill the histograms
  h_amplitude->Fill(amplitude);
  h_amplitude_low->Fill(amplitude);

  // Fill the arrays
  // m_amplitude_array.push_back(amplitude);

};

//------------------------------------------------------------------------------

void PMT::writeHist()
{

  // Write the histograms to file
  h_amplitude->Write();
  h_amplitude_low->Write();

};

//------------------------------------------------------------------------------

void PMT::clean()
{
  // m_amplitude_array.clear();
};

//------------------------------------------------------------------------------
