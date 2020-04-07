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
      /*
      TCanvas* cTest = new TCanvas();
      for(int i=0; i<chargehist.size(); i++){
        if(! chargehist[i]) continue;
        if(i==0) chargehist[i]->Draw();
        else chargehist[i]->Draw("same");
        chargehist[i]->SetLineColor(i+1);
      }
      cTest->SaveAs("testplot.png");
      */
    }

  private:

    vector<TH1D*> chargehist;
    vector<double> voltages;
    int m_board;
    int m_ch;
};
