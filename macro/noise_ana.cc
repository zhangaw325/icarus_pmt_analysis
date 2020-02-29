////////////////////////////////////////////////////////////////////////////////
// Macro to analyze the PMT noise
// Process the noise analyis for the filename given as argument. Edit the input
// filename with the file you prefer. Note: the output will always be in the format
// "noise_run%d_00%d.root". you must change manually the variables run and subrun for now.
//
// Usage: root -l loadLib.cc noise_ana.cc
// Out of the box, it will produce a root file in the macro folder called "noise_run1264_1.root"
// with the ADC distribution of all the pmts in data_dl1_run1264_1_20200227T235326-decoded.root
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

using namespace std;



void noise_ana()
{

  // CHANGE THE FILENAME BELOW IF WANT TO USE ONE OTHER FILE FOR YOUR ANALYSIS
  string filename="/media/disk_a/ICARUS/PMTDataMonitoring/data_dl1_run1264_1_20200227T235326-decoded.root";
  int run=1264; //<< TO BE CHANGE ACCORDING TO THE FILENAME
  int subrun=1; //<< TO BE CHANGE ACCORDING TO THE FILENAME

  const int nboards=1;    // << This type of files in /PMTDataMonitoring/ have only one board. CHANGE ONLY IF YOU ARE LOOKING AT ONE OTHER TYPE OF FILE!
  const int nchannels=16; // << DON'T CHANGE

  //****************************************************************************
  // Input
  //

  TChain *tchain =  new TChain("caenv1730dump/events");

  // HERE you filename is loaded inside the macro to be processed!
  tchain->Add(filename.c_str());

  std::vector<std::vector<uint16_t> > *data=0; //unsigned short
  tchain->SetBranchAddress("fWvfmsVec", &data);


  //****************************************************************************
  // Now we create and initialize an histogram per each board and channel
  // This is an example, you can create more histograms using this model

  TH1D *h_pmt_rms[nboards][nchannels];

  for(int board=0; board<nboards; board++)
  {
    for(int channel=0; channel<nchannels; channel++)
    {
      //DEFINE HERE THE HISTOGRAMS. EACH HISTOGRAM SHOULD HAVE A DIFFERENT NAME
      //TO AVOID CONFUSION WHEN YOU SAVE THEM TO FILE.
      char hname[100];
      sprintf(hname, "hnoise_run%d_00%d_board%d_channel%d", run, subrun, board, channel);
      h_pmt_rms[board][channel]= new TH1D(hname, hname, 20, -20, 20);
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

        //LOOP OVER THE RAW WAVEFORM, UNCOMMENT THE LOOP IF NEED TO DO IT
        //for( float entry : waveform->getWaveform() )
        //{
        //  //DO SOMETHING
        //}

        // Loop over the entries of the waveform after the baseline subtraction
        // and fill the entries of the histogram h_pmt_rms
        for( float entry : waveform->getWaveform() )
        {
          h_pmt_rms[board][channel]->Fill( entry );
        }

      } // channel
    } // boards
  } // event



  //****************************************************************************
  // Output: create an output TFile and write the histogram in it
  //

  char ofilename[100]; sprintf(ofilename, "noise_run%d_%d.root", run, subrun);
  TFile ofile(ofilename, "RECREATE"); ofile.cd();

  for(int board=0; board<nboards; board++)
  {
    for(int channel=0; channel<nchannels; channel++)
    {
      //Write the histogram to the output TFILE
      h_pmt_rms[board][channel]->Write();
    }
  }

  ofile.Close();

  cout << "All done" << endl;

  //****************************************************************************
  //
  // Now you can open the output file "noise_run%d_00%d.root" using the command:
  // root -l "noise_run%d_00%d.root". With the TBrowser you should be able to
  // visualize all the histograms. Alternatively you can create one other root
  // macro which import this TFILE and make TGraphs with the quantities you find
  // more relevant.
  //
  //****************************************************************************

} //end macro
