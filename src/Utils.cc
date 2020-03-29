////////////////////////////////////////////////////////////////////////////////
// Utils.h function definition
// maito: ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Utils.h"

vector<vector<string> > utils::CSVReader::getData()
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
 };

//------------------------------------------------------------------------------
