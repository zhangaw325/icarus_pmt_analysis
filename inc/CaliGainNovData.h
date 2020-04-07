#include "Pmt.h"
#include "TH1F.h"

class CaliGainNovData{

  public:

    CaliGainNovData(int board, int channel):m_board(board),m_ch(channel){
      chargehist.resize(0);
      voltages.resize(0);
    }
    ~CaliGainNovData(){ for(int i=0; i<chargehist.size(); i++) delete chargehist[i]; }
    void loadOneHist(TH1D* hist, double hv){
        chargehist.push_back(hist);
        voltages.push_back(hv);
    }
    void plotHist(){
      cout<<chargehist.size()<<endl;
      
      char text[100]; 
      TCanvas* cTest = new TCanvas();
      gStyle->SetOptStat(0);
      for(int i=0; i<chargehist.size(); i++){
        if(! chargehist[i]) continue;
        if(i==0) chargehist[i]->Draw();
        else chargehist[i]->Draw("same");
        chargehist[i]->SetLineColor(i+1);
        cout<<chargehist[i]->GetNbinsX()<<"\t"<<voltages[i]<<endl;
      }
      TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
      leg->SetHeader("HV values");
      for(int i=0; i<chargehist.size(); i++){
        sprintf(text, "%.1f", voltages[i]);
        leg->AddEntry(chargehist[i],text,"lp");
      }
      leg->Draw();
      sprintf(text,"chargeDist_bd%d_ch%d.png",m_board, m_ch);
      cTest->SaveAs(text);
      cTest->Close();
      
    }

  private:

    vector<TH1D*> chargehist;
    vector<double> voltages;
    int m_board;
    int m_ch;
};
