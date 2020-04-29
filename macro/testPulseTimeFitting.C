void testPulseTimeFitting(){
  TFile* f = new TFile("Waveforms_run207.root","read");
  TH1D* h = (TH1D*)f->Get("Run207-Subrun2-Event0-Board0-Channel1_hist");

  TCanvas* c = new TCanvas();
  h->Draw();
}
