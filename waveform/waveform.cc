////////////////////////////////////////////////////////////////////////////////
//
// Analyze the single waveform and place the results in the PMT histograms
// Usage:
//      ./build/waveform -i /path/to/decoded_file.root -o /path/to/pmt_file.root
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Geometry.h"
#include "Run.h"
#include "Waveform.h"
#include "Pmt.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <fstream>
#include <sstream>

using namespace std;

int main( int argc, char* argv[] )
{

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

  // THE RUN CLASS SHOULD ASSIGN THE RIGTH RUN AND SUBRUN NUMBERS ACCORDING TO
  // THE FILENAME GIVEN
  RUN my_run(ifilename);
  int run=my_run.getRun();
  int subrun=my_run.getSubrun();

  // We get the optical channel and find if it is straighforwardly invalid
  int optical_channel = my_run.getOpticalChannel();
  if(optical_channel < 0)
  {
    cout<<"ERROR! Check if the optical channel is correctly identified!"<< endl;
  }

  // Print out what the macro is doing so far
  cout << "Processing run: " << run << " and subrun: " << subrun << endl;
  cout << "Optical channel is: " << optical_channel << endl;

  // Initialize the PMT
  PMT* pmts[geo::nboards][geo::nchannels];
  for( int board=0; board<geo::nboards; board++ )
  {
    for( int channel=0; channel<geo::nchannels; channel++ )
    {
      pmts[board][channel] = new PMT( run, board, channel );
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
  // EVENT LOOPER

  for(int event=0; event<nevents; event++)
  {

    if ( event % 100 == 0 )
    {
      cout<<"["<<int((float(event)/float(nevents))*100)<<"%] Completed"<<endl;
    }

    tevents->GetEvent(event);

    for(int board=0; board<geo::nboards; board++)
    {
      for(int channel=0; channel<geo::nchannels; channel++)
      {

        // Apply a selection using the optical channel
        if( pmts[board][channel]->isIlluminated(optical_channel) ){ continue; }

        // Create the waveform object and load raw data in
        Waveform *waveform = new Waveform(run, subrun ,event, board, channel);
        waveform->loadData((*data).at(channel+geo::nchannels*board));

        if(!waveform->isValidWaveform()){ continue; }

        // Extract the signal: hasPulse method enables all the pulse quantities
        if(waveform->hasPulse(1.0))
        {
          pmts[board][channel]->loadWaveform(waveform);
        }

        // SOMETHING ELSE TO DO IN THIS LOOP?

      }
    }
  } // end events

  ifile->Close();

  //****************************************************************************
  // NOW SAVE THE INFO STORED IN THE PMT CLASS

  // Open TFile
  TFile *ofile = new TFile(ofilename.c_str(), "RECREATE");
  if( !ofile->IsOpen() )
  {
    cout << "File "+ ofilename + " not open!" << endl;
    return 1;
  }
  ofile->cd();

  for( int board=0; board<geo::nboards; board++ )
  {
    for( int channel=0; channel<geo::nchannels; channel++ )
    {
      pmts[board][channel]->writeHist();
    }
  }

  ofile->Close();

  cout << "All done!" << endl;
  return 0;
}
