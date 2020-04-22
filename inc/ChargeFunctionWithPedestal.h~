#ifndef  __ChargeFunctionWithPedestal_H
#define __ChargeFunctionWithPedestal_H

#include "TMath.h"

class ChargeFunctionWithPedestal {
  public:
    // use constructor to customize the function object
    double operator() (double* x, double * par){
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
};

class ChargeFunctionWithPedestal_expGaus { // expoentially-modified Gaussian
  public:
    // use constructor to customize the function object
    double operator() (double* x, double * par){
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
};

class myChargeFunctionWithRealPedestal { // expoentially-modified Gaussian
  public:
    // use constructor to customize the function object
    double operator() (double* x, double * par){
       // following Milind's paper, arXiv: 1909.05373
      double mu = par[0];  // mean npe
      double q = par[1];   // charge of spe
      double sigma = par[2]; // width of spe
      double q0 = par[3];   // mean of backgroud
      double sigma0 = par[4]; // width of backgroud
      double c0 = par[5];    // exp. component in the backgroud
      double w = par[6]; // the weight: probability of dark pulse 
      double a = par[7]; // final scale factor or amplitude

      double pdf1 =0;
      pdf1 = TMath::Exp(-1.0*mu) * ( (1.0-w)*myGausPDF(x[0],q0,sigma0)
                                          +w*myExpGausPDF(x[0],q0,sigma0,c0) );

      double pdf2 = mu*TMath::Exp(-1.0*mu) * ( (1.0-w)*myGausPDF(x[0],q+q0,TMath::Sqrt    (sigma*sigma+sigma0*sigma0))
                                            + w*myExpGausPDF(x[0], q+q0, TMath::Sqrt(sigma*sigma+sigma0*sigma0), c0) );

      double pdf3 = 0;
      for(int k=2; k<=20; k++){
        pdf3 += TMath::Power(mu, k)*TMath::Exp(-1.0*mu)/TMath::Factorial(k)
           * ( (1.0-w)*myGausPDF(x[0],q*k+q0, TMath::Sqrt(sigma*sigma*k+sigma0*sigma0))
               + w*myExpGausPDF(x[0],q*k+q0, TMath::Sqrt(sigma*sigma*k+sigma0*sigma0), c0)) ;
      }

      return a*(pdf1+pdf2+pdf3);
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
};

#endif
