////////////////////////////////////////////////////////////////////////////////
// Event Class definition
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef  __EVENT_H
#define __EVENT_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "Pmt.h"

using namespace std;

class Event
{
  public:

    typedef vector<PMT*> Pmts_t;

    Event();
    Event(int run, int subrun, int event, int nevents );
    Event(int run, int subrun, int event, int nevents, Pmts_t pmt_vector );
    ~Event();

    void setPmtVector( Pmts_t pmt_vector ){ m_pmt_vector = pmt_vector; };

    void setPmt( PMT* pmt ){ m_pmt_vector.push_back(pmt); };

    Pmts_t getPmtVector(){ return m_pmt_vector; };

    // Selective getters
    Pmts_t getPmtsBoard( int board );
    Pmts_t getPmtsChannel( int channel );
    Pmts_t getValidPmts( double n_sigma );

    void clean(){ m_pmt_vector.clear(); };

    int m_run;
    int m_subrun;
    int m_event;
    int m_nevents;

  private:

    Pmts_t m_pmt_vector;
};

#endif //__EVENT_H
