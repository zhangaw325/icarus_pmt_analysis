double idealresponse(double* x, double * par){
      //function implementation using class data members
      double mu = par[0];
      double q = par[1];
      double sigma = par[2];
      double amplitude = par[3];
      double sum=0;
      for(int n=1; n<100; n++){
        sum += (TMath::Power(mu,n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n)
                *TMath::Exp(-1.0*(x[0]-q*n)*(x[0]-q*n)
                  /(2.0*n*sigma*sigma))/(sigma*TMath::Sqrt(2.0*TMath::Pi()*n)));
      }
      return amplitude * sum;
}

void test(){
  TFile* f = new TFile("gain_curves_result_run217.root");
  //TGraphErrors* g = (TGraphErrors*)f->Get("gain_bd0_ch0");
  //g->Draw("AP");
/*
  TF1* fGain = new TF1("fGain","[0]*TMath::Power(x,[1])",1,2000);
  fGain->SetParNames("a","b");
  g->Fit("fGain","QR","",1150, 1270);
  for(int i=0; i<5;i++){
    double a=fGain->GetParameter(0);
    double b=fGain->GetParameter(1);
    fGain->SetParameters(a,b);
    g->Fit("fGain","QR","",1150, 1270);
  }
*/

  TCanvas* c = (TCanvas*)f->Get("cChargeDis_bd0_ch7");
  TH1D* h = (TH1D*)c->GetPrimitive("hLowQ_run213_board0_ch7");
  TCanvas* c1 = new TCanvas();
  h->Draw();
  h->GetYaxis()->SetRangeUser(0,600);
  h->GetXaxis()->SetRangeUser(-1,15);
  TF1* func = (TF1*)c->GetPrimitive("idealfunc0");
  func->Draw("same");

  TF1* f2 = new TF1("f2",idealresponse,0,100,4);
  f2->SetParameters(0.6, 3, 0.6, 800);
  h->Fit("f2","RQ","",0.8,15);
  double par[4];
  for(int k=0; k<3;k++){
    for(int i=0;i<4;i++){
      par[i]=f2->GetParameter(i);
    }
    f2->SetParameters(par[0],par[1],par[2],par[3]);
    h->Fit("f2","RQ","",0.8,15);
  }
}
