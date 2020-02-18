////////////////////////////////////////////////////////////////////////////////
// Macro that interfaces the output of the CAEN decoder script and produces
// a file format compatible with this small analysis framework
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"
//#include "Event.h"
#include "Pmt.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"

#include <fstream>
#include <sstream>

using namespace std;

void addTChain( string list_name, TChain *tchain, int ncount)
{
  /* Read the list of filenames to add to the TChain */
  fstream input( list_name, ios::in );
  int count =0;
  if (input.is_open()){   //checking whether the file is open
      string tp;
      while(getline(input, tp) && count < ncount ){ //read data from file object and put it into string.
         tchain->Add( tp.c_str() );
         count++;
      }
      input.close(); //close the file object.
   }

}

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

void load_pmt_tchain()
{
  //****************************************************************************
  // General definitions ( ..import from DB )

  int run=1067;
  int subrun=11;
  int nboards=12;
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

  addTChain( "../production/new_decoded_west_m50.list", tchain, 1);

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

        if(waveform->hasPulse(1.0))
        {
          PMT *pmt_temp = locatePMT(run, board, channel, pmts);
          pmt_temp->loadWaveform(waveform);
        }

      } // end channel
    } // end boards

  } // end events

  //****************************************************************************
  // Now save the PMTs

  // Open TFile
  TFile ofile("../data/pmt_west_m50.root", "RECREATE");
  ofile.cd();

  for(auto pmt_it : pmts){ pmt_it->writeHist(); }

  ofile.Close();

  cout << "All done ..." << endl;

} //end macro
