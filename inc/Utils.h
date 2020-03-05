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

  // Parses through csv file line by line and returns the data
  // in vector of vector of strings.

  vector<vector<string> > CSVReader::getData()
  {
	   ifstream file(m_filename);
	   vector<vector<string> > dataList;

	   string line = "";
	   // Iterate through each line and split the content using delimeter
	   while (getline(file, line))
	   {
       istringstream s_stream(line);
       string token;
		   vector<string> vec;
       while (getline(s_stream, token, m_delimeter))
       {
         //cout << token << endl;
         vec.push_back(token);
       }
       //cout << vec << endl;
       dataList.push_back(vec);
	   }
	   // Close the File
	   file.close();

	   return dataList;
   }



////////////////////////////////////////////////////////////////////////////////

} //end namespace

#endif //UTILS_H
