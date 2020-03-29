////////////////////////////////////////////////////////////////////////////////
// This is the space where spurios functions can be defined
//
// maito: ascarpell@bnl.gov
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>

#include "Geometry.h"

using namespace std;

namespace utils
{

  //////////////////////////////////////////////////////////////////////////////

  class CSVReader
  {
	  string m_filename;
	  char m_delimeter;

    public:
	     CSVReader(string filename, char delm = ',')
			 : m_filename(filename), m_delimeter(delm)
	     { }

	     // Function to fetch data from a CSV File
	     vector<vector<string> > getData();
  };

  //////////////////////////////////////////////////////////////////////////////

  //void get_rundata( string filename, int & run, int & subrun )
  //{
    // NB: assumes file with structure run*_***.root
  //  int suffix[3];
  //  suffix[0] = filename.find("run")+3;
  //  suffix[1] = filename.find("_");
  //  suffix[2] = filename.find(".root");
  //
  //  run = stoi(filename.substr(suffix[0], suffix[1]-suffix[0]));
  //  subrun = stoi(filename.substr(suffix[1]+1, suffix[2]-suffix[1]-1));

  //  return;
 //}

  //////////////////////////////////////////////////////////////////////////////

} //end namespace

#endif //UTILS_H
