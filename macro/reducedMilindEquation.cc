double mygainResponseReduced(double p0,
                                       double p1, double p2, double p3, double p4, double p5, double p6, double p7,
                                       double p8, double p9, double p10, double p11, double p12, double p13, double p14,
                                       double p15, double p16, double p17, double p18, double p19, double p20, double p21){
  // following Milind's paper, arXiv: 1909.05373
  const int nVs = 3; // number of voltage points
  int n=7; // for each charge distribution, number of pars (not counting npe)
  double mu = p0;//par[0];  // mean npe
  double q[nVs]; //= par[1];   // charge of spe
  double sigma[nVs];// = par[2]; // width of spe
  double q0[nVs]; //= par[3];   // mean of backgroud
  double sigma0[nVs]; //= par[4]; // width of backgroud
  double c0[nVs]; //= par[5];    // exp. component in the backgroud
  double w[nVs]; //= par[6]; // the weight: probability of dark pulse 
  double a[nVs]; //= par[7]; // final scale factor or amplitude
  /*
  for(int i=0; i<nVs;i++){
    q[i] = 0;//par[i*n+1];
    sigma[i] = 0;//par[i*n+2];
    q0[i] = 0;//par[i*n+3];
    sigma0[i] = 0;//par[i*n+4];
    c0[i] = 0;//par[i*n+5];
    w[i] = 0;//par[i*n+6];
    a[i] = 0;//par[i*n+7];
  }
  */

  q[0] = p1;     q[1] = p8;     q[2] = p15;
  sigma[0]=p2;   sigma[1]=p9;   sigma[2]=p16;
  q0[0]=p3;      q0[1]=p10;     q0[2]=p17;
  sigma0[0]=p4;  sigma0[1]=p11; sigma0[2]=p18;
  c0[0]=p5;      c0[1]=p12;     c0[2]=p19;
  w[0]=p6;       w[1]=p13;      w[2]=p20;
  a[0]=p7;       a[1]=p14;      a[2]=p21;

  double pdf1[nVs], pdf2[nVs], pdf3[nVs];
  double chi2[nVs];
  double prediction = 0;

  int ndf = 0;
  for(int i=0; i<nVs; i++){
      pdf1[i]=0; pdf2[i]=0; pdf3[i]=0;
      chi2[i] =0;
      int nbins = (int)ds.ghist[i]->GetNbinsX();

      prediction = 0;
      for(int j=0; j<nbins; j++){
          double x = ds.ghist[i]->GetBinCenter(j+1);
          double y = ds.ghist[i]->GetBinContent(j+1);

          //if(x<-1 || x>20) continue;

          pdf1[i] = TMath::Exp(-1.0*mu) * ( (1.0-w[i])*myGausPDF(x,q0[i],sigma0[i])
                                      +w[i]*myExpGausPDF(x,q0[i],sigma0[i],c0[i]) );
          // test nan value
          if(TMath::IsNaN(q0[i]) || TMath::IsNaN(sigma0[i]) )
            cout<<"weried place 1 for nan value. "<<x<<"\t"<<q0[i] <<"\t"<<sigma0[i]<<"\t"<<c0[i]<<"\t"<<c0[i]<<endl;

          pdf2[i] = mu*TMath::Exp(-1.0*mu) * ( (1.0-w[i])*myGausPDF(x,q[i]+q0[i],TMath::Sqrt(sigma[i]*sigma[i]+sigma0[i]*sigma0[i]))
                           + w[i]*myExpGausPDF(x, q[i]+q0[i], TMath::Sqrt(sigma[i]*sigma[i]+sigma0[i]*sigma0[i]), c0[i]) );
          // test nan value
          if( TMath::IsNaN(q[i]+q0[i]) || TMath::IsNaN(TMath::Sqrt(sigma[i]*sigma[i]+sigma0[i]*sigma0[i])) ){
             cout<<"wiried place 2 for nan: "<<x<<"\t"<<mu<<"\t"<<q[i]<<"\t"<<q0[i]<<"\t"
                                             <<sigma[i]<<"\t"<<sigma0[i]<<"\t"<<c0[i]<<"\t"<<w[i]<<endl;
          }

          for(int k=2; k<=100; k++){
              pdf3[i] = TMath::Power(mu, k)*TMath::Exp(-1.0*mu)/TMath::Factorial(k)
                      * ( myGausPDF(x,q[i]*k+q0[i], TMath::Sqrt(sigma[i]*sigma[i]*k+sigma0[i]*sigma0[i])) ) ;

              // test nan value
             if( TMath::IsNaN(q[i]*k+q0[i]) || TMath::IsNaN(TMath::Sqrt(sigma[i]*sigma[i]*k+sigma0[i]*sigma0[i])) ){
               cout<<"wiried place 3 for nan: "<<k<<"\t"<<x<<"\t"<<mu<<"\t"<<"\t"<<q[i]<<"\t"<<q0[i]<<"\t"<<sigma[i]<<"\t"<<sigma0[i]<<"\t"<<c0[i]<<"\t"<<w[i]<<endl;
             }

          }

          prediction = a[i]*(pdf1[i]+pdf2[i]+pdf3[i]);

          // apply the poisson likelihood
          //if(prediction!=0){
            if(y!=0){
              chi2[i] += TMath::Power( y-prediction, 2)/y;
              ndf++;
            }
            //else chi2[i] += 2.0*prediction;
            //cout<<"chi2 "<<i<<chi2[i]<<"\t"<<prediction<<"\t"<<y<<endl;
          //}
          //else{
            // just for test
            //cout<<"in fit p=0 : "<<"bin "<<j<<" , x= "<<x<<" , y= "<<y<<endl;
          //}
          //cout<<"thoise pdfs "<<i<<"\t"<<pdf1[i]<<"\t"<<pdf2[i]<<"\t"<<pdf3[i]<<endl;
      }
  }
  ds.ndf=ndf;
  //cout<<"return value "<<chi2[0]+chi2[1]+chi2[2]<<endl;
  return chi2[0]+chi2[1]+chi2[2];

}

void mygainfcnReduced(int &npar, double * deriv, double &f, double * par, int iflag){
  //cout<<"mygainfcn called with "<<npar<<" parameters : " <<endl;
  //for(int i=0; i<22;i++) cout<<par[i]<<"\t";
  //cout<<endl;
  f = mygainResponseReduced(par[0],
                             par[1],par[2],par[3],par[4],par[5],par[6],par[7],
                             par[8],par[9],par[10],par[11],par[12],par[13],par[14],
                             par[15],par[16],par[17],par[18],par[19],par[20],par[21]);
  //cout<<"f="<<f<<endl;
}

double myChargeFuncReduced(double * x, double * par){
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

  if(TMath::IsNaN(q)) cout<<"q is nan in myChargeFuncReduced"<<endl;
  if(TMath::IsNaN(sigma)) cout<<"sigma is nan in myChargeFuncReduced"<<endl;
  if(TMath::IsNaN(q0)) cout<<"q0 is nan in myChargeFuncReduced"<<endl;
  if(TMath::IsNaN(sigma0)) cout<<"sigma0 is nan in myChargeFuncReduced"<<endl;

  double pdf2 = mu*TMath::Exp(-1.0*mu) * ( (1.0-w)*myGausPDF(x[0],q+q0,TMath::Sqrt(sigma*sigma+sigma0*sigma0))
                                            + w*myExpGausPDF(x[0], q+q0, TMath::Sqrt(sigma*sigma+sigma0*sigma0), c0) );
  double pdf3 = 0;
  for(int k=2; k<100; k++){
    pdf3 += TMath::Power(mu, k)*TMath::Exp(-1.0*mu)/TMath::Factorial(k)
           * ( myGausPDF(x[0],q*k+q0, TMath::Sqrt(sigma*sigma*k+sigma0*sigma0)) ) ;
  }

  return a*(pdf1+pdf2+pdf3);
}

void   fitWithMilindEqReduced(){

    const int nVs = 3;
    const int npars=22;
    double vstart[npars];
    double step[npars];
    double minVal[npars]={0.01, 
                           0.01, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0,
                           0.01, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0,
                           0.01, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0};
    double maxVal[npars]={50, 
                           10, 1, 0.5, 1, 1e5, 1, 1e7,
                           10, 1, 0.5, 1, 1e5, 1, 1e7,
                           10, 1, 0.5, 1, 1e5, 1, 1e7};
    string parName[npars]={"npe",
          "q11", "#sigma_{11}", "q01", "#sigma_{01}", "c01", "w1","a1", 
          "q12", "#sigma_{12}", "q02", "#sigma_{02}", "c02", "w2","a2", 
          "q13", "#sigma_{13}", "q03", "#sigma_{03}", "c03", "w3","a3"};

   double vresult[npars];
   double verrors[npars];

    // initial values to the parameters
    vstart[0] = ds.ghist[2]->GetMean()/1.6;
    if(vstart[0]<0.5) vstart[0]=0.5;
    step[0] = 0.1;
    for(int i=0; i<nVs; i++){
        step[i*7+1] = 0.1; // step for q1
        step[i*7+2] = 0.05; // step for sigma1
        step[i*7+3] = 0.02; // step for q0
        step[i*7+4] = 0.02; // step for sigma0
        step[i*7+5] = 0.5; // step for c0;
        step[i*7+6] = 0.1; // step for w;
        step[i*7+7] = 100; // step for a

        vstart[i*7+1] = 1; // start for q1
        vstart[i*7+2] = 0.5; // start for sigma1
        vstart[i*7+3] = 0; // start for q0
        vstart[i*7+4] = 0.1; // start for sigma0
        vstart[i*7+5] = 1; // start for c0
        vstart[i*7+6] = 0.3; // start for w
        vstart[i*7+7] = ds.ghist[i]->Integral()*TMath::Exp(-1.0*vstart[0]); // start for a
    }

    TF1* hfit = new TF1("myfit", myChargeFuncReduced,-1,50,8);
    for(int i=0 ; i<3; i++){
      hfit->SetParameters(vstart[0], vstart[i*7+1], vstart[i*7+2], vstart[i*7+3], vstart[i*7+4], vstart[i*7+5], vstart[i*7+6], vstart[i*7+7] );
      ds.ghist[i]->Fit("myfit","RQ", "", -1, 20);
      cout<<"fit individually, I get these parameters"<<endl;
      cout<<"chi2/ndf "<< hfit->GetChisquare()<<"\t"<<hfit->GetNDF()<<endl;
      double thispar[8];
      for(int j=0; j<8; j++){
        thispar[j] = hfit->GetParameter(j);
        cout<<thispar[j]<<"\t"<<hfit->GetParError(j)<<"\t";
      }
      cout<<endl;
      vstart[0] = thispar[0];
      for(int j=0; j<7; j++){
        vstart[i*7+j+1] = thispar[j+1];
      }
    }

    // set the fitter: use TMinuit
    //gSystem->Load("libMinuit");
    TMinuit* minimizer = new TMinuit(npars); // was ds.npars here
    // the fit function is set
    //cout<<"Set the minimizer with "<<npar<<" parameters "<<endl;
    minimizer->SetFCN(mygainfcnReduced);
    //  
    for(int i=0; i<npars; i++){
      minimizer->DefineParameter(i, parName[i].c_str(), vstart[i], step[i], 
                             minVal[i], maxVal[i]);
    }
    // do the fit
    minimizer->Command("SET PRINT -1");//
    Int_t ierrs=0;
    minimizer->mnexcm("SET NOWarnings", 0, 0, ierrs);
    int fitstatus= minimizer->Migrad();

    // get parameters
    for(int i=0; i<npars; i++){
      minimizer->GetParameter(i, vresult[i], verrors[i]);
    }
    cout<<"check fit status "<<fitstatus<<endl;
    // if fit is unsuccessful, try again with updated initial values
    if(fitstatus==4){
      for(int i=0; i<npars; i++){
        minimizer->DefineParameter(i, parName[i].c_str(), vresult[i], step[i]/2.0, 
                               minVal[i], maxVal[i]);
      }
      fitstatus = minimizer->Migrad();
      cout<<"status of the second fit "<<fitstatus<<endl;
      for(int i=0; i<npars; i++){
        minimizer->GetParameter(i, vresult[i], verrors[i]);
      }
    }  

    // here we get the result from the fit
    double amin,edm,errdef;
    int nvpar,nparx,icstat;
    minimizer->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    cout<<amin<<"\t"<<ds.ndf-22<<"\t"<<edm<<"\t"<<errdef<<"\t"<<nvpar<<"\t"<<nparx<<"\t"<<icstat<<endl;  
    cout<<"---fit results "<<endl;
    for(int i=0; i<npars; i++){
      cout<<parName[i]<<"\t"<<vresult[i]<<"\t"<<verrors[i]<<endl;
    }

    //m_fitstatus = fitstatus;
    //minChi2 = amin;
  
   TCanvas* c = new TCanvas();
   ds.ghist[0]->Draw(); ds.ghist[0]->SetLineColor(1);
   ds.ghist[0]->GetXaxis()->SetRangeUser(-1,15);
   ds.ghist[1]->Draw("same"); ds.ghist[1]->SetLineColor(2);
   ds.ghist[2]->Draw("same"); ds.ghist[2]->SetLineColor(3);

   TF1* func[3];
   char text[50];
   for(int i=0; i<3; i++){
      sprintf(text, "func%d", i);
      func[i] = new TF1(text,myChargeFuncReduced,-1,50,8);
      func[i]->SetParameters(vresult[0], 
                        vresult[i*7+1], vresult[i*7+2], vresult[i*7+3], vresult[i*7+4], vresult[i*7+5], vresult[i*7+6], vresult[i*7+7]);
      func[i]->Draw("same"); func[i]->SetLineColor(i+1); func[i]->SetLineStyle(2);
   }
  
}
