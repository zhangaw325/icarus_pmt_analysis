#include "Pmt.h"
#include "TH1F.h"

class CaliGainNovData{

  public:

    CaliGainNovData(int board, int channel);

    ~CaliGainNovData();

    void loadOneHist(TH1D* hist, double hv);
    void plotHist();

  private:

    vector<TH1D*> chargehist;
    vector<double> voltages;
    int m_board;
    int m_ch;
};
