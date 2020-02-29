////////////////////////////////////////////////////////////////////////////////
// Macro to analyze the PMT signal
// Process the pulse analyis for the filename given as argument. Edit the input
// filename with the file you prefer. Note: the output will always be in the format
// "pmt_run%d_00%d.root". you must change manually the variables run and subrun for now.
//
// Usage: root -l loadLib.cc load_pmt_tchain.cc
// Out of the box, it will produce a root file in the macro folder called "pmt_run1067_1.root"
// with the ADC distribution of all the pmts in
// after baseline subtraction.
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"
#include "Pmt.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"

#include <fstream>
#include <sstream>

using namespace std;



void load_pmt_tchain()
{

  // CHANGE THE FILENAME BELOW IF WANT TO USE ONE OTHER FILE FOR YOUR ANALYSIS
  string filename="/media/disk_a/ICARUS/PMT_gainData_20200205_decoded/data_dl1_run1067_1_20200204T225822_dl3_decoded.root";
  int run=1067; // << CHANGE ACCORDING TO THE FILE USED
  int subrun=1; // << CHANGE ACCORDING TO THE FILE USED

  const int nboards=12; // << DO NOT CHANGE!
  const int nchannels=16; // << DO NOT CHANGE!


  //****************************************************************************
  // Input
  //

  TChain *tchain =  new TChain("caenv1730dump/events");

  // HERE you filename is loaded inside the macro to be processed!
  tchain->Add(filename.c_str());

  std::vector<std::vector<uint16_t> > *data=0; //unsigned short
  tchain->SetBranchAddress("fWvfmsVec", &data);


  //****************************************************************************
  // Here we create and initialize a vector with the class PMT. If you want you
  // can define also other objects (like TH1D classes) in here. See macro
  // noise_ana.cc with an example
  //

  PMT* pmts[nboards][nchannels];
  for( int board=0; board<nboards; board++ )
  {
    for( int channel=0; channel<nchannels; channel++ )
    {
      pmts[board][channel] = new PMT( run, board, channel );
    }
  }


  //****************************************************************************
  // Now it is time to loop over the events
  //


  for(int e=0; e<tchain->GetEntries(); e++)
  {
    cout << "Processing event: " << e << endl;

    // WE TAKE THE EVENT
    tchain->GetEvent(e);

    for(int board=0; board<nboards; board++)
    {
      for(int channel=0; channel<nchannels; channel++)
      {

        // Create the WAVEFORM OBJECT
        Waveform *waveform = new Waveform(run, subrun ,e, board, channel);


        // ONCE THE DATA ARE LOADED IN THE WAVEFORM OBJECT, BASELINE IS
        // AUTOMATICALLY SUBSTRACTED. YOU CAN USE THE FUNCTION "getWaveform()"
        // TO EXTRACT THE WAVEFORM VECTOR AFTER BASELINE SUBTRACTION OR
        // "getRawWaveform()" BEFORE BASELINE SUBTRACTION.
        waveform->loadData((*data).at(channel+nchannels*board));

        // THIS PART IS SOMETHING THAT YOU PROBABLY HAVE TO MODIFY:
        // WHAT THE CODE DOES HERE IS CHECKING IF THE WAFEFORM SATISFY MY CRITERIA
        // TO FIND A PULSE. IF YES, THIS WAVEFORM IS SAVED INSIDE THE PMT OBJECT
        if(waveform->hasPulse(1.0))
        {
          pmts[board][channel]->loadWaveform(waveform);
        }

      } // end channel
    } // end boards

  } // end events



  //****************************************************************************
  // Output: create an output TFile and write the histogram in it
  //

  char ofilename[100]; sprintf(ofilename, "pmt_run%d_%d.root", run, subrun);
  TFile ofile(ofilename, "RECREATE"); ofile.cd();

  for( int board=0; board<nboards; board++ )
  {
    for( int channel=0; channel<nchannels; channel++ )
    {
      pmts[board][channel]->writeHist();
    }
  }

  ofile.Close();

  cout << "All done" << endl;

  //****************************************************************************
  //
  // For now you want to focus mostly on the waveform class and change the logic
  // behind the pulse finder.
  //
  //****************************************************************************

} //end macro
