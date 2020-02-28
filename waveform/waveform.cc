////////////////////////////////////////////////////////////////////////////////
//
// Analyze the single waveform and place the results in the PMT histograms
// Usage: ./build/waveform -i /path/to/decoded_file.root -o /path/to/pmt_file.root
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"
#include "Pmt.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <fstream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

PMT* locatePMT(int run, int board, int channel , vector<PMT*> pmts)
{
  // Locate a PMT in the PMT array given the channel and board number

  auto find = [&]( PMT *pmt )
  {
    return (pmt->getRun()==run && pmt->getBoard()==board
                                                 && pmt->getChannel()==channel);
  };

  auto it = std::find_if( pmts.begin(), pmts.end(), find );

  return *it;
}

//------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{

  //****************************************************************************
  // General definitions ( ..import from DB ? )

  int run=1067; // Useless
  int subrun=11; // Useless
  const int nboards=12;
  const int nchannels=16;

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
  // Handle the input and output

  // Parse the arguments
  string ifilename, ofilename;
  for ( int i=1; i<argc; i=i+2 )
  {
    if      ( string(argv[i]) == "-i" ) ifilename = argv[i+1];
    else if ( string(argv[i]) == "-o" ) ofilename = argv[i+1];
    else {
      cout << "Unknown option " << argv[i+1] << endl;
      return 1;
    }
  }

  // Open TFile
  TFile* ifile = new TFile(ifilename.c_str(), "READ");
  if( !ifile->IsOpen() )
  {
    cout << "File "+ ifilename + " not open!" << endl;
    return 1;
  }

  // Get the TTres
  TTree* tevents = (TTree*)ifile->Get("caenv1730dump/events");
  int nevents = tevents->GetEntries();

  // Set Branch address
  std::vector<std::vector<unsigned short> > *data=0; //unsigned short
  tevents->SetBranchAddress("fWvfmsVec", &data);

  //****************************************************************************
  // Event looper

  for(int event=0; event<nevents; event++)
  {

    if ( event % 100 == 0 )
    {
      cout<<"["<<int((float(event)/float(nevents))*100)<<"%] Completed"<<endl;
    }

    tevents->GetEvent(event);

    const int n_samples = (*data)[0].size();
    const int n_boards = (*data).size()/nchannels;

    for(int board=0; board<n_boards; board++)
    {
      for(int channel=0; channel<nchannels; channel++)
      {

        // Create the waveform object and load raw data in
        Waveform *waveform = new Waveform(run, subrun ,event, board, channel);
        waveform->loadData((*data).at(channel+nchannels*board));

        // Other operations on the single waveform
        // ...

        // Extract the signal: hasPulse method enables all the pulse quantities
        if(waveform->hasPulse(1.0))
        {
          PMT *pmt_temp = locatePMT(run, board, channel, pmts);
          pmt_temp->loadWaveform(waveform);
        }

      } // end channel
    } // end board

    // Other event operations (if applicable)
    // ...

  } // end events

  ifile->Close();

  //****************************************************************************
  // Now save the PMTs

  // Open TFile
  TFile *ofile = new TFile(ofilename.c_str(), "RECREATE");
  if( !ofile->IsOpen() )
  {
    cout << "File "+ ofilename + " not open!" << endl;
    return 1;
  }
  ofile->cd();

  // Now loop over the pmt array and fill the quantities in file
  for(auto pmtit : pmts)
  {

    // Other single pmt analysis (if applicable)
    // ...

    // Fill the file with the histogram only
    pmtit->writeHist();

  }

  ofile->Close();

  // Is there anything else to do?
  // ...

  cout << "All done ..." << endl;
  return 0;
}
