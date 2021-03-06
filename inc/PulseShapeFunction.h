#ifndef  __PulseShapeFunction_H
#define __PulseShapeFunction_H

#include "TMath.h"

class PulseShapeFunction_ExpGaus {
  public:
    // use constructor to customize the function object
    double operator() (double* x, double * par){
      double t0 = par[0];
      //double mu = par[1];
      double w = par[1];
      double c = par[2];
      double a = par[3];

      double t = x[0];

      return a*c/2.0*TMath::Exp(c*c*w*w/2.0)*TMath::Exp(-1.0*c*(t-t0))
              * TMath::Erfc( 1.0/1.414* (c*w-(t-t0)/w) );
    }
};

class PulseShapeFunction_LogNormal {
  public:
    // use constructor to customize the function object
    double operator() (double* x, double * par){
      double t0 = par[0];
      double mu = par[1];
      double w = par[2];
      double a = par[3];

      double t = x[0];
      if(t<=mu) return par[4]; //baseline
      else
        return par[4]+a/TMath::Sqrt(2.0*TMath::Pi()) / (w*(t-t0))
               * TMath::Exp( 1.0/2.0 * TMath::Power( (TMath::Log10(t-t0)-mu)/w ,2) );
    }
};

class PulseShapeFunction_TwoExpoentials{
  public:
    double operator() (double* x, double* par){
      double t=x[0];
      double t0 = par[0];
      double tau = par[1];
      double tauS = par[2];
      double a = par[3];
      double a_off = par[4];

      double value = 0;
      if(t<t0) value = a_off;
      else{
        value = a/(tau-tauS)*( TMath::Exp(-1.0*(t-t0)/tauS) - TMath::Exp(-1.0*(t-t0)/tau) ) + a_off;
      }
      return -1.0*value; // just inverted because the waveform is inverted
    }
};

#endif
