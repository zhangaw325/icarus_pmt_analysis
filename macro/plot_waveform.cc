////////////////////////////////////////////////////////////////////////////////
//
// Quick macro to visualize a waveform given a filename, an event, a board and a
// channel number. Just edit the argument of plot_waveform() as you prefer
// Usage: root -l loadLib.cc plot_waveform.cc
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#include "Waveform.h"
#include "Pmt.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

using namespace std;

void plot_waveform( string filename="/media/disk_a/ICARUS/PMTDataMonitoring/data_dl1_run1264_1_20200227T235326-decoded.root", int event=1, int my_board = 0, int my_channel = 1)
{

  // they can stay set at 0 both of them for now
  int run=0;
  int subrun=0;
  int nchannels=16;

  // Open TFile
  TFile* ifile = new TFile(filename.c_str(), "READ");
  cout << "Open TFile"+filename << endl;

  // Get the TTres
  TTree* tevents = (TTree*)ifile->Get("caenv1730dump/events");
  int nevents = tevents->GetEntries();

  // Set Branch address
  std::vector<std::vector<uint16_t> > *data=0; //unsigned short
  tevents->SetBranchAddress("fWvfmsVec", &data);
  tevents->GetEvent(event);

   // Create the PMT object
  Waveform *waveform = new Waveform(run, subrun ,event, my_board, my_channel);
  waveform->loadData((*data).at(my_channel+nchannels*my_board));

  char cname[100];
  sprintf(cname, "Waveform-Run%d-Subrun%d-Event%d-Board%d-Channel%d", run, subrun ,event, my_board, my_channel);
  printf("%s \n", cname);
  TH1D* hist = waveform->getWaveformHist();
  TCanvas *c = new TCanvas(cname, cname, 1200, 400); c->cd();
  hist->GetXaxis()->SetRangeUser(1500, 3500);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleOffset(0.33);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetRangeUser(-20, 20);
  hist->Draw("hist");
  cout << waveform->getBaselineMean() << " " << waveform->getBaselineWidth() << endl;

} //end macro
