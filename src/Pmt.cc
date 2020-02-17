////////////////////////////////////////////////////////////////////////////////
// PMT Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Pmt.h"

PMT::PMT( )
{};

//------------------------------------------------------------------------------

PMT::PMT( int run, int board, int channel )
  : m_run(run)
  , m_board(board)
  , m_channel(channel)
{};

//------------------------------------------------------------------------------

PMT::~PMT()
{};

//------------------------------------------------------------------------------

void PMT::loadWaveform( Waveform *waveform )
{
  m_waveform_array.push_back( waveform );
};

//------------------------------------------------------------------------------

bool PMT::hasWaveform()
{
  return ( m_waveform_array.size() > 0 );
};

//------------------------------------------------------------------------------

vector<double> PMT::getAmplitudes()
{
  vector<double> m_amplitude_array;

  for( auto waveform : m_waveform_array )
  {
    if(waveform->hasPulse(1.0))
    {
      m_amplitude_array.push_back( waveform->getAmplitude() );
    }
  }

  return m_amplitude_array;

}

//------------------------------------------------------------------------------

void PMT::clean()
{
  m_waveform_array.clear();
};

//------------------------------------------------------------------------------
