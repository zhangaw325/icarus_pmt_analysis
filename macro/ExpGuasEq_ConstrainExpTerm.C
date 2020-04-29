double ResponseWithPedestal_expGaus_constrainNpe(double npe,
                                      double Q1, double w1, double a1,            
                                      double Q2, double w2, double a2,
                                      double Q3, double w3, double a3,
                                      double pedm1, double pedw1, double pedamp1, double pedtau1,
                                      double pedm2, double pedw2, double pedamp2, double pedtau2,
                                      double pedm3, double pedw3, double pedamp3, double pedtau3)
{
  //cout<<"idealResponse is called ..."<<endl;
  const int npoints = (int) ds.ghist.size(); //3;//this->getNbOfDataPoints();
  //cout<<"check "<<global_hist.size()<<endl;
  if(npoints==0) {cout<<"fitting no data in CaliGainNovData:myfcn"<<endl; exit(0);}

  // I have these parameters:
  //double npe = mu; // mean npe
  double q[npoints]; q[0]=Q1; q[1]=Q2; q[2]=Q3;
  double sigma[npoints]; sigma[0]=w1; sigma[1]=w2; sigma[2]=w3;
  double amp[npoints]; amp[0]=a1; amp[1]=a2; amp[2]=a3;
  double q0[npoints]; q0[0]=pedm1; q0[1]=pedm2; q0[2]=pedm3; 
  double w0[npoints]; w0[0]=pedw1; w0[1]=pedw2; w0[2]=pedw3;
  double a0[npoints]; a0[0]=pedamp1; a0[1]=pedamp2; a0[2]=pedamp3;
  double tau0[npoints]; tau0[0]=pedtau1; tau0[1]=pedtau2; tau0[2]=pedtau3;

  // compute the chi2
  int ndf=0;
  double chi2[3]; 
  for(int i=0; i<npoints; i++){
    chi2[i]=0.0;
    int nbins = (int)ds.ghist[i]->GetNbinsX();//200;

    for(int j=0; j<nbins-1; j++){
      double x = ds.ghist[i]->GetBinCenter(j+1);
      double y = ds.ghist[i]->GetBinContent(j+1);
      double prediction = 0;

      // ideal response term
      double p2 = 0;
      for(int k=1;k<100;k++){
          p2 +=    (TMath::Power(npe,k)*TMath::Exp(-1.0*npe)/TMath::Factorial(k)
            *TMath::Exp(-1.0*(x-q[i]*k)*(x-q[i]*k) 
              /(2.0*k*sigma[i]*sigma[i]))/(sigma[i]*TMath::Sqrt(2.0*TMath::Pi()*k)) ) ;
      }

      prediction = (p2*amp[i] 
                         + a0[i]*w0[i]/tau0[i]*TMath::Sqrt(TMath::Pi()/2.0)
                           *TMath::Exp(0.5*(w0[i]*w0[i])/(tau0[i]*tau0[i])-(x-q0[i])/tau0[i])
                           *TMath::Erfc(TMath::Sqrt(1./2.)*(w0[i]/tau0[i]-(x-q0[i])/w0[i])) ) *TMath::Exp(-1.0*npe);

      if(y!=0){
        chi2[i] += TMath::Power( y-prediction , 2)/y; // the chi2 function
        ndf++;
      }
    }
  }

  ds.ndf=ndf;
  return chi2[0]+chi2[1]+chi2[2];

}


void myfcnWithPedestal_expGaus_constrainNpe(int &npar, double * deriv, double &f, double * par, int iflag){
  f = ResponseWithPedestal_expGaus_constrainNpe(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], par[10], par[11], par[12], par[13], par[14], par[15], par[16], par[17], par[18], par[19], par[20], par[21]);
}

double myChargeFunc_expGaus_constrainNpe (double* x, double * par){
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

      return sum*TMath::Exp(-1.0*mu) ;
}

void fitWithNormalExpGausBkg_constrainNpe(){

    // define number of parameters and their initial values

    int bin1[3];
    int bin2[3];
    for(int i=0;i<3;i++){
      bin1[i] = ds.ghist[i]->GetXaxis()->FindBin(-0.2);
      bin2[i] = ds.ghist[i]->GetXaxis()->FindBin(0.2);
    }
    //int bin3 = lowchargehist[0]->GetXaxis()->FindBin(0.5);
    //int bin4 = lowchargehist[0]->GetNbinsX();


    //int n = getNbOfDataPoints() ;
    const int n_pars = 22;//ds.npars_lowcharge;//10;//1 + 3*n ;
    double vstart[n_pars]={20, 
                       1.0, 0.5, ds.ghist[0]->Integral()/400.0*40,
                       1.0, 0.5, ds.ghist[1]->Integral()/400.0*40,
                       1.0, 0.5, ds.ghist[2]->Integral()/400.0*40,
                       0.0, 0.1, ds.ghist[0]->Integral(bin1[0],bin2[0])/10.0*4, 0.1,
                       0.0, 0.1, ds.ghist[1]->Integral(bin1[1],bin2[1])/10.0*4, 0.1,
                       0.0, 0.1, ds.ghist[2]->Integral(bin1[2],bin2[2])/10.0*4, 0.1};
    double step[n_pars]={0.02,
                          0.05, 0.1, 5, 
                          0.05, 0.1, 5, 
                          0.05, 0.1, 5, 
                          0.005, 0.01, 10, 0.01,
                          0.005, 0.01, 10, 0.01,
                          0.005, 0.01, 10, 0.01};
    double minVal[n_pars]={0,
                           0.01, 0.001, 0, 
                           0.01, 0.001, 0, 
                           0.01, 0.001, 0,
                           -0.5, 0.001, 0, 0,
                           -0.5, 0.001, 0, 0,
                           -0.5, 0.001, 0, 0};
    double maxVal[n_pars]={10, 10, 10,ds.ghist[0]->Integral()/4.0, 
                               10, 10,ds.ghist[1]->Integral()/4.0, 
                               10, 10,ds.ghist[2]->Integral()/4.0,
                               0.5, 0.5, ds.ghist[0]->Integral(bin1[0],bin2[0])*10.0, 0.5,
                               0.5, 0.5, ds.ghist[1]->Integral(bin1[1],bin2[1])*10.0, 0.5,
                               0.5, 0.5, ds.ghist[2]->Integral(bin1[2],bin2[2])*10.0, 0.5 };
    string parName[n_pars]={"npe",
                            "q1", "w1", "a1",
                            "q2", "w2", "a2",
                            "q3", "w3", "a3",
                            "pedm1", "pedw1", "peda1", "pedt1",
                            "pedm2", "pedw2", "peda2", "pedt2",
                            "pedm3", "pedw3", "peda3", "pedt3"
                           };

   double vresult[n_pars];
   double verrors[n_pars];

    // reset some of the initial parameters
    for(int i=0; i< (int)ds.ghist.size(); i++){
      //step[i*3+1] = 0.05;
      vstart[0] = ds.ghist[i]->GetMean()/1.6;
      //if(vstart[0]<1) vstart[0]=1;
    }

    // set the fitter: use TMinuit
    TMinuit* minimizer = new TMinuit(n_pars);
    // the fit function is set
    minimizer->SetFCN(myfcnWithPedestal_expGaus_constrainNpe);
    //  
    for(int i=0; i<n_pars; i++){
      minimizer->DefineParameter(i, parName[i].c_str(), vstart[i], step[i], 
                             minVal[i], maxVal[i]);
    }
    // do the fit
    minimizer->Command("SET PRINT -1");//
    //minimizer->SetPrintLevel(0);
    Int_t ierrs=0;
    minimizer->mnexcm("SET NOWarnings", 0, 0, ierrs);
    int fitstatus= minimizer->Migrad();

    // get parameters
    for(int i=0; i<n_pars; i++){
      minimizer->GetParameter(i, vresult[i], verrors[i]);
    }
    cout<<"check fit status "<<fitstatus<<endl;
    if(fitstatus==4){
      for(int i=0; i<n_pars; i++){
        minimizer->DefineParameter(i, parName[i].c_str(), vresult[i], step[i]/2.0, 
                               minVal[i], maxVal[i]);
      }
      fitstatus = minimizer->Migrad();
      cout<<"status of the second fit "<<fitstatus<<endl;
      for(int i=0; i<n_pars; i++){
        minimizer->GetParameter(i, vresult[i], verrors[i]);
      }
    }  

    // here we get the result from the fit
    double amin,edm,errdef;
    int nvpar,nparx,icstat;
    minimizer->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    cout<<amin<<"\t"<<ds.ndf-22<<"\t"<<edm<<"\t"<<errdef<<"\t"<<nvpar<<"\t"<<nparx<<"\t"<<icstat<<endl;  
    cout<<"---fit results "<<endl;
    for(int i=0; i<n_pars; i++){
      cout<<parName[i]<<"\t"<<vresult[i]<<"\t"<<verrors[i]<<endl;
    }

    // plot

   TCanvas* c = new TCanvas();
   ds.ghist[0]->Draw(); ds.ghist[0]->SetLineColor(1);
   ds.ghist[0]->GetXaxis()->SetRangeUser(-1,15);
   ds.ghist[1]->Draw("same"); ds.ghist[1]->SetLineColor(2);
   ds.ghist[2]->Draw("same"); ds.ghist[2]->SetLineColor(3);

   TF1* func[3];
   char text[50];
   for(int i=0; i<3; i++){
      sprintf(text, "func%d", i);
      func[i] = new TF1(text,myChargeFunc_expGaus_constrainNpe,-1,50,8);
      func[i]->SetParameters(vresult[0], 
                        vresult[i*3+1], vresult[i*3+2], vresult[i*3+3], vresult[i*4+10], vresult[i*4+11], vresult[i*4+12], vresult[i*4+13]);
      func[i]->Draw("same"); func[i]->SetLineColor(i+1); func[i]->SetLineStyle(2);
   }    

}

