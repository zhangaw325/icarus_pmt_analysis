#include <iostream>
#include <stdio.h>


namespace utils
{

////////////////////////////////////////////////////////////////////////////////

  void get_rundata( string filename, int & run, int & subrun )
  {
	  // NB: assumes file with structure data_dl1_run*_*_20200204T225822_dl3_decoded.root
    int suffix[3];
    int end = filename.find(".root");

    suffix[0] = filename.find("run")+3;
    suffix[1] = suffix[0]+filename.substr(suffix[0], end-suffix[0]).find("_");
    suffix[2] = suffix[1]+filename.substr(suffix[1]+1, end-suffix[1]).find("_");

    run = std::stoi(filename.substr(suffix[0], suffix[1]-suffix[0]));
    subrun = std::stoi(filename.substr(suffix[1]+1, suffix[2]-suffix[1]));

    return ;
  }

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////

} //end namespace
