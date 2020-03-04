////////////////////////////////////////////////////////////////////////////////
// RUN Methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Run.h"

RUN::RUN()
{
  // Produce association between run-subrun and optical channel
  this->loadDatabase();
};

//------------------------------------------------------------------------------

RUN::RUN(int run, int subrun)
  : m_run(run)
  , m_subrun(subrun)
{
  // Produce association between run-subrun and optical channel
  this->loadDatabase();
};

//------------------------------------------------------------------------------

RUN::RUN(string filename)
  : m_run_file(filename)
{
  m_run=0; m_subrun=0;
  this->getRundata(m_run_file, m_run, m_subrun);

  // Produce association between run-subrun and optical channel
  this->loadDatabase();
};

//------------------------------------------------------------------------------

RUN::~RUN()
{ };

//------------------------------------------------------------------------------


void RUN::loadDatabase()
{
  // Produce association between run-subrun and optical channel

  m_optical_channel=18;
}


//------------------------------------------------------------------------------

void RUN::getRundata(string filename, int & run, int & subrun)
{
  // NB: It assumes filename with structure *runRUN_SUBRUN_*.root
  int suffix[3];
  int end = filename.find(".root");

  suffix[0] = filename.find("run")+3;
  suffix[1] = suffix[0]+filename.substr(suffix[0], end-suffix[0]).find("_");
  suffix[2] = suffix[1]+filename.substr(suffix[1]+1, end-suffix[1]).find("_");

  run = std::stoi(filename.substr(suffix[0], suffix[1]-suffix[0]));
  subrun = std::stoi(filename.substr(suffix[1]+1, suffix[2]-suffix[1]));

  return ;
};

//------------------------------------------------------------------------------
