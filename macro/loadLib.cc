////////////////////////////////////////////////////////////////////////////////
// Basic library loader into ROOT
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

void loadLib(){

  //gStyle->SetPalette(kRainBow);
  //gStyle->SetOptFit(11111);

  // Target a build folder
  gSystem->SetBuildDir("obj",true);

  // For some reasons ROOT has not a dictionary for vector of vectors
  gInterpreter->GenerateDictionary("vector<vector<unsigned short>>", "vector");

  char *eigen_install = getenv("EIGEN_INSTALL");
  if(eigen_install == NULL)
  {
    cout << " Set the EIGEN_INSTALL variable! " << endl;
  }
  gInterpreter->AddIncludePath(eigen_install);

  // Point to the location of the other libraries
  // It assumes you use loadLib inside the macro folder
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  string currentdir( cwd );

  string includepath=currentdir+"/../inc/";
  gInterpreter->AddIncludePath(includepath.c_str());

  // Compile the custom classes
  gROOT->LoadMacro((currentdir+"/../src/Utils.cc+").c_str());
  gROOT->LoadMacro((currentdir+"/../src/Run.cc+").c_str());
  gROOT->LoadMacro((currentdir+"/../src/Waveform.cc+").c_str());
  gROOT->LoadMacro((currentdir+"/../src/Pmt.cc+").c_str());

  #define __INITIALIZED__

}//end loadLib
