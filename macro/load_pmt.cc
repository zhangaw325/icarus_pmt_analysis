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

void load_pmt()
{

  // Input
  int run=1077;
  int subrun=23;

  string filename="../data/data_dl2_run1077_23_20200205T183344_dl4.root";
  string treename="caenv1730dump/events";


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

  vector<Event*> events;

  // Loop over the events
  for(int e=1; e<3; e++)
  {
    cout << "Processing event: " << e << endl;

    tevents->GetEvent(e);

    // Initalize the Event object
    Event *event = new Event(run, subrun, e, n_events);

    const int n_channels = 16;
    const int n_samples = (*data)[0].size();
    const int n_boards = (*data).size()/n_channels;

    for(int board=0; board<n_boards; board++)
    {
      for(int channel=0; channel<n_channels; channel++)
      {

        // Create the PMT object
        PMT *pmt = new PMT(run, subrun ,e, board, channel, n_samples,
                                          (*data).at(channel+n_channels*board));

        // Save the PMT in the event
        event->setPmt(pmt);

      } // end channel
    } // end boards

    // Save the event in a vector
    events.push_back(event);

  } // end events

  cout << "All done ..." << endl;

} //end macro
