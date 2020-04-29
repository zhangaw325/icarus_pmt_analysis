double myGausPDF(double x, double mean, double sigma){
  // test output
  double value = 1.0/TMath::Sqrt(2.0*TMath::Pi()*sigma*sigma)*TMath::Exp(-1.0/2.0*TMath::Power( (x-mean)/sigma ,2));

  //if( TMath::IsNaN(value) ){
  //  cout<<"nan in myGausPDF"<<endl;
  //  cout<<x<<"\t"<<mean<<"\t"<<sigma<<endl;
  //}

  return value;
}

double myExpGausPDF(double x, double mean, double sigma, double c0){
  double value = sigma/c0 * TMath::Sqrt(TMath::Pi()/2.0)
         * TMath::Exp( 1./2.*(sigma/c0)*(sigma/c0)- (x-mean)/c0 )
         * TMath::Erfc( 1./TMath::Sqrt(2.) * (sigma/c0 - (x-mean)/sigma) );
  // test output
  //if( TMath::IsNaN(value) ){
  //  cout<<"nan in myExpGausPDF"<<endl;
  //  cout<<x<<"\t"<<mean<<"\t"<<sigma<<"\t"<<c0<<endl;
  //}
  return value;
  //return c0*0.5*TMath::Exp(0.5*c0*(2.0*mean+c0*sigma*sigma-2.0*x))
  //       * TMath::Erfc( (mean+c0*sigma*sigma-x)/sigma/1.414);
}

