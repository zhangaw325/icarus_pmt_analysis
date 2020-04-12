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

double responsewithpedestal(double* x, double * par){
      //function implementation using class data members
      double mu = par[0];
      double q = par[1];
      double sigma = par[2];
      double amplitude = par[3];
      double sum=0;
      for(int n=1; n<100; n++){
        sum += amplitude* (TMath::Power(mu,n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n)
                *TMath::Exp(-1.0*(x[0]-q*n)*(x[0]-q*n)
                  /(2.0*n*sigma*sigma))/(sigma*TMath::Sqrt(2.0*TMath::Pi()*n))) ;
                       
      }

     sum = sum + par[6]/TMath::Sqrt(2.0*TMath::Pi()*par[5]*par[5])
                           *TMath::Exp(-1.0*(x[0]-par[4])*(x[0]-par[4])/(2.0*par[5]*par[5])); // added pedestal

      return sum ;
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
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  h->Draw();
  h->GetYaxis()->SetRangeUser(0,600);
  h->GetXaxis()->SetRangeUser(-1,30);
  TF1* func = (TF1*)c->GetPrimitive("idealfunc0");
  func->Draw("same");

  TF1* f2 = new TF1("f2",responsewithpedestal,-1,100,7);
  //f2->SetParameters(0.8, 0.4, 0.14, 400, 0, 0.14, 700);
  //f2->Draw();
  f2->SetParameters(0.6, 3, 0.6, 50, 0, 0.1, 700);

  h->Fit("f2","RQ","",-0.5,30);
  double par[7];
  for(int k=0; k<3;k++){
    for(int i=0;i<7;i++){
      par[i]=f2->GetParameter(i);
    }
    f2->SetParameters(par[0],par[1],par[2],par[3], par[4], par[5], par[6]);
    h->Fit("f2","RQ","",-1,30);
  }


  
}
