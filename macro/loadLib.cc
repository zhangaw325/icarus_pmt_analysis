////////////////////////////////////////////////////////////////////////////////
// Basic library loader into ROOT
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

void loadLib(){

  char cwd[1024];
  getcwd(cwd, sizeof(cwd));

  //gStyle->SetPalette(kRainBow);
  //gStyle->SetOptFit(11111);

  std::string currentdir( cwd );
  std::string includepath=currentdir+"/../inc/";

  gSystem->SetBuildDir("obj",true);
  gInterpreter->AddIncludePath(includepath.c_str());
  gROOT->LoadMacro((currentdir+"/../src/Event.cc+").c_str());
  gROOT->LoadMacro((currentdir+"/../src/Pmt.cc+").c_str());
  gROOT->LoadMacro((currentdir+"/../src/Waveform.cc+").c_str());

  //gROOT->LoadMacro((currentdir+"/source/311Lib.cc+").c_str());
  //gROOT->LoadMacro((currentdir+"/source/311style.cc+").c_str());

  #define __INITIALIZED__

}//end loadLib
