////////////////////////////////////////////////////////////////////////////////
// RUN Class definition
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef  __RUN_H
#define __RUN_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "Utils.h"

using namespace std;

class RUN
{

  //typedef map<int, int> Runinfo_t;
  //typedef map<int, Runinfo_t> Opticalinfo_t;

  public:

    RUN();
    RUN(int run, int subrun);
    RUN(string filename);
    ~RUN();

    void loadFile(string filename){ m_run_file=filename; };
    void getRundata( string filename, int & run, int & subrun );
    int getRun(){ return m_run; }
    int getSubrun(){ return m_subrun; }
    int getOpticalChannel();

  private:

    void loadDatabase();

    string m_run_file;
    int m_run;
    int m_subrun;
    int m_optical_channel=-1;
    vector<vector<string> > m_dataList;

    //Opticalinfo_t m_optical_map;

};


#endif //__RUN_H
