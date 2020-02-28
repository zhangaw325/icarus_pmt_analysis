////////////////////////////////////////////////////////////////////////////////
// Macro to analyze the PMT noise
// Process the noise analyis for the filename given as argument. Edit the function
// argument with the file you prefer. Note: the output will always be in the format
// "noise_run%d_00%d.root". you must change manually the variables run and subrun for now.
//
// Usage: root -l loadLib.cc noise_ana.cc
// Out of the box it will produce a root file in the macro folder called "noise_run0_000.root"
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

//------------------------------------------------------------------------------

PMT* locatePMT(int run, int board, int channel , vector<PMT*> pmts)
{
  int index = -1;
  auto find = [&]( PMT *pmt )
  {
    return (pmt->getRun()==run && pmt->getBoard()==board
                                                 && pmt->getChannel()==channel);
  };

  auto it = std::find_if( pmts.begin(), pmts.end(), find );
  return *it;
}

//------------------------------------------------------------------------------

void noise_ana(string filename="/media/disk_a/ICARUS/PMTDataMonitoring/data_dl1_run1264_1_20200227T235326-decoded.root")
{

  int run=0.0;
  int subrun=0.0;
  int nboards=11;
  int nchannels=16;

  //****************************************************************************
  // Data structures definitions

  vector<PMT*> pmts;

  for( int board=0; board<nboards; board++ )
  {
    for( int channel=0; channel<nchannels; channel++ )
    {
      PMT *pmt = new PMT( run, board, channel );
      pmts.push_back(pmt);
    }
  }

  //****************************************************************************
  // Input

  string treename="caenv1730dump/events";
  TChain *tchain =  new TChain(treename.c_str());

  tchain->Add(filename.c_str());

  // Set Branch address
  std::vector<std::vector<uint16_t> > *data=0; //unsigned short
  tchain->SetBranchAddress("fWvfmsVec", &data);

  //****************************************************************************
  // Event loop

  // Loop over the events
  for(int e=0; e<tchain->GetEntries(); e++)
  {
    cout << "Processing event: " << e << endl;

    tchain->GetEvent(e);

    const int n_channels = 16;
    const int n_samples = (*data)[0].size();
    const int n_boards = (*data).size()/n_channels;

    for(int board=0; board<n_boards; board++)
    {
      for(int channel=0; channel<n_channels; channel++)
      {

        // Create the PMT object
        Waveform *waveform = new Waveform(run, subrun ,e, board, channel);
        waveform->loadData((*data).at(channel+n_channels*board));

        PMT *pmt_temp = locatePMT(run, board, channel, pmts);
        pmt_temp->loadWaveform(waveform);

      } // end channel
    } // end boards

  } // end events

  //****************************************************************************
  // Now save the PMTs

  // Open TFile (remember to change the output filename if necessary)
  char ofilename[100]; sprintf(ofilename, "noise_run%d_00%d.root", run, subrun);
  TFile ofile(ofilename, "RECREATE");
  ofile.cd();

  for(auto pmt_it : pmts)
  {
    cout << pmt_it->getBoard() << " " << pmt_it->getChannel() << endl;
    TH1D* h_noise = pmt_it->getNoiseHist();
    h_noise->Write();
  }

  ofile.Close();

  cout << "All done" << endl;

} //end macro
