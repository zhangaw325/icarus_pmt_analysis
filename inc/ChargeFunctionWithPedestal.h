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

#endif
