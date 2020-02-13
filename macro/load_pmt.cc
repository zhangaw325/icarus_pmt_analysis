////////////////////////////////////////////////////////////////////////////////
// Macro that interfaces the output of the CAEN decoder script and produces
// a file format compatible with this small analysis framework
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Pmt.h"
#include "Event.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

using namespace std;

void printWaveform(int run, int event, int board, int channel, vector<PMT*> pmts )
{
    char name[100];
    sprintf(name, "Run%d-Event%d-Board%d-Channel%d", run, event, board, channel);
    TCanvas *c = new TCanvas(name, name, 600, 400);

    auto it = find_if( pmts.begin(), pmts.end(), [ & ] ( PMT *pmt ){
                                return pmt->find(run, event, board, channel); } );

    if( it != pmts.end() )
    {
      PMT *pmt = *it;

      TH1D *hist = pmt->getWaveformHist();

      hist->SetName(name);
      hist->SetTitle(name);
      hist->Draw("hist");
    }
    else
    {
      cout << "Impossible to find PMTs  with board " << board
                                            << " and channel " << channel << endl;
      return;
    }
}

//==============================================================================

void load_pmt(string filename="../data/data_dl2_run1077_23_20200205T183344_dl4.root",
                                                  string treename="caenv1730dump/events")
{
  // Open TFile
  TFile* ifile = new TFile(filename.c_str(), "READ");
  cout << "Open TFile"+filename << endl;

  // Get the TTres
  TTree* tevents = (TTree*)ifile->Get(treename.c_str());
  int n_events = tevents->GetEntries();
  cout << "TTree"+treename << " has " << n_events << " events" << endl;

  // Set Branch address
  std::vector<std::vector<uint16_t> > *data=0; //unsigned short
  tevents->SetBranchAddress("fWvfmsVec", &data);


  vector<PMT*> pmtVector;

  // Loop over the events
  for(int event=1; event<n_events; event++)
  {
    cout << "Processing event: " << event << endl;

    tevents->GetEvent(event);

    const int n_channels = 16;
    const int n_samples = (*data)[0].size();
    const int n_boards = (*data).size()/n_channels;

    // Loop over the boards
    for(int board=0; board<n_boards; board++)
    {

      // Loop over the PMTs
      for(int channel=0; channel<n_channels; channel++)
      {

        // Fill the PMT object
        PMT *pmt = new PMT(1, event, board, channel, n_samples,
                                          (*data).at(channel+n_channels*board));
        pmtVector.push_back(pmt);

      } // end channel
    } // end boards

    // Fill the Event Object
    // ..

  } // end events

  // ..

  printWaveform(1, 1, 6, 1, pmtVector);

  cout << "All done ..." << endl;

} //end macro
