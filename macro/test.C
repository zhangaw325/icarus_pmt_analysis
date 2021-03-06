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

double responsewithpedestal_expGaus(double* x, double * par){
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

     // add the exponentially modified guassian function:
     // https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
     // par[6]: amplitude
     // par[4]: mean
     // par[5]: sigma
     // par[7]: exponent time (I use the alternative form)

     sum = sum + par[6]*par[5]/par[7]*TMath::Sqrt(TMath::Pi()/2.0)*TMath::Exp(0.5*(par[5]*par[5])/(par[7]*par[7])-(x[0]-par[4])/par[7])*TMath::Erfc(TMath::Sqrt(1./2.)*(par[5]/par[7]-(x[0]-par[4])/par[5]));

      return sum ;
}

double myGausPDF(double x, double mean, double sigma){
  return 1.0/TMath::Sqrt(2.0*TMath::Pi()*sigma*sigma)*TMath::Exp(-1.0/2.0*TMath::Power( (x-mean)/sigma ,2));
}

double myExpGausPDF(double x, double mean, double sigma, double c0){
  return sigma/c0 * TMath::Sqrt(TMath::Pi()/2.0)
         * TMath::Exp( 1./2.*(sigma/c0)*(sigma/c0)- (x-mean)/c0 )
         * TMath::Erfc( 1./TMath::Sqrt(2.) * (sigma/c0 - (x-mean)/sigma) );
  //return c0*0.5*TMath::Exp(0.5*c0*(2.0*mean+c0*sigma*sigma-2.0*x))
  //       * TMath::Erfc( (mean+c0*sigma*sigma-x)/sigma/1.414);
}

double myChargeFunc(double * x, double * par){
  // following Milind's paper, arXiv: 1909.05373
  double mu = par[0];  // mean npe
  double q = par[1];   // charge of spe
  double sigma = par[2]; // width of spe
  double q0 = par[3];   // mean of backgroud
  double sigma0 = par[4]; // width of backgroud
  double c0 = par[5];    // exp. component in the backgroud
  double w = par[6]; // the weight: probability of dark pulse 
  double a = par[7]; // final scale factor or amplitude

  double pdf1 = TMath::Exp(-1.0*mu) * ( (1.0-w)*myGausPDF(x[0],q0,sigma0)
                                     +w*myExpGausPDF(x[0],q0,sigma0,c0) );

  if(TMath::IsNaN(q)) cout<<"q is nan"<<endl;
  if(TMath::IsNaN(sigma)) cout<<"sigma is nan"<<endl;
  if(TMath::IsNaN(q0)) cout<<"q0 is nan"<<endl;
  if(TMath::IsNaN(sigma0)) cout<<"sigma0 is nan"<<endl;

  double pdf2 = mu*TMath::Exp(-1.0*mu) * ( (1.0-w)*myGausPDF(x[0],q+q0,TMath::Sqrt(sigma*sigma+sigma0*sigma0))
                                            + w*myExpGausPDF(x[0], q+q0, TMath::Sqrt(sigma*sigma+sigma0*sigma0), c0) );
  double pdf3 = 0;
  for(int k=2; k<20; k++){
    pdf3 += TMath::Power(mu, k)*TMath::Exp(-1.0*mu)/TMath::Factorial(k)
           * ( (1.0-w)*myGausPDF(x[0],q*k+q0, TMath::Sqrt(sigma*sigma*k+sigma0*sigma0))
               + w*myExpGausPDF(x[0],q*k+q0, TMath::Sqrt(sigma*sigma*k+sigma0*sigma0), c0) ) ;
  }

  return a*(pdf1+pdf2+pdf3);
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
  TH1D* h = (TH1D*)c->GetPrimitive("hLowQ_run217_board0_ch7");
  TCanvas* c1 = new TCanvas();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  h->Draw();
  //h->GetYaxis()->SetRangeUser(0,600);
  h->GetXaxis()->SetRangeUser(-1,100);
  //TF1* func = (TF1*)c->GetPrimitive("idealfunc0");
  //func->Draw("same");

  //TF1* f2 = new TF1("f2",responsewithpedestal,-1,100,7);
  //f2->SetParameters(0.8, 0.4, 0.14, 400, 0, 0.14, 700);
  //f2->Draw();
  //f2->SetParameters(0.6, 3, 0.6, 50, 0, 0.1, 700);

  //TF1* f2 = new TF1("f2",responsewithpedestal_expGaus,-1,100,8);
  //f2->SetParameters(0.6, 1, 0.6, 1600, 0, 0.1, 2400, 0.1);

  TF1* f2 = new TF1("f2",myChargeFunc,-1,100,8);
  f2->SetParameters(0.5, 1, 0.6, 0, 0.1, 10, 0.3, 1600);
  //f2->FixParameter(3,0);
  //f2->FixParameter(4,0);
  //TCanvas* c2 = new TCanvas();
  //f2->Draw();
 

  h->Fit("f2","R","",-0.5,100);
  double par[8];
  for(int k=0; k<3;k++){
    for(int i=0;i<8;i++){
      //if(i==3 || i==4) continue;
      par[i]=f2->GetParameter(i);
      f2->SetParameter(i, par[i]);
    }
    //f2->SetParameters(par[0],par[1],par[2],par[3], par[4], par[5], par[6]);
    h->Fit("f2","R","",-1,100);
  }


  
}
