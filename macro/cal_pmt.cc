////////////////////////////////////////////////////////////////////////////////
//
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"
//#include "Event.h"
#include "Pmt.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

const int nboards=12;
const int nchannels=16;

TH1D *h_amplitude[nboards][nchannels];

void cal_pmt()
{
  //****************************************************************************
  TCanvas c("amplitude", "amplitude", 600, 500); c.Print("amplitude.pdf[");

  //****************************************************************************
  // Read the file in

  TFile* ifile = new TFile("test_pmt.root", "READ");
  TTree *pmt_tree = (TTree*)ifile->Get("pmts");
  PMT *pmt = new PMT(); pmt_tree->SetBranchAddress( "pmt", &pmt );

  //****************************************************************************
  // Loop over pmts

  for(int index=0; index<pmt_tree->GetEntries(); index++ )
  {
    pmt_tree->GetEntry(index);

    //Initialize the histograms
    char hname[100];

    int run = pmt->getRun();
    int board = pmt->getBoard();
    int channel = pmt->getChannel();

    vector<Waveform*> waveforms = pmt->getWaveforms();
    if( waveforms.size() < 10 ){ continue; }

    printf("Run %d Board%d Channel%d \n", run, board, channel);

    sprintf(hname, "Run%d_Board%d_Channel%d_amplitude", run, board, channel);
    h_amplitude[board][channel] = new TH1D(hname, hname, 100, -2000, 0);

    for( auto waveform : waveforms )
    {
      if(waveform->hasPulse(1.0))
      {
        h_amplitude[board][channel]->Fill( waveform->getAmplitude() );
      }
    }

    //
    h_amplitude[board][channel]->Draw();
    c.Print("amplitude.pdf");

  } // end pmts

  c.Print("amplitude.pdf]");

} //end macro
