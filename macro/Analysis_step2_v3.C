#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <bits/stdc++.h> 

#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
//#include <TSubString.h>
#include <TDirectory.h>

//using namespacce std;

double pulsefunc(double* x, double * par){
    double v=0;
    if( x[0]<par[3] ) v=0+par[4];
    else v = par[0]/(par[2]-par[1])*( TMath::Exp(-(x[0]-par[3])/par[1]) - TMath::Exp(-(x[0]-par[3])/par[2]) ) + par[4];
    return v;
}

double pulsefunc2(double* x, double * par){
    double v=0;
    if( x[0]<par[2] ) v=0+par[3];
    else v = par[0]/(par[1]*par[1])*(x[0]-par[2])*( TMath::Exp(-(x[0]-par[2])/par[1]) ) + par[3]; // - TMath::Exp(-(x[0]-par[3])/par[2]) ) + par[4];
    return v;
}

//single pe. waveform model: Soren JETTER et al: PMT waveform modeling at the Daya Bay experiment, Chinese Physics C 2012(36):733 
double pulsefunc1(double* x, double* par){
    double v=0;
    if( x[0]<par[3] ) v=0;
    else v=par[0]*TMath::Exp( -0.5*TMath::Power( TMath::Log( (x[0]-par[3]) /par[1])/par[2] ,2) );
    return v;
}

int Analysis_step2_v3(){ //Analysis_step2

    //std::string filename = "run207_01.root";
    TString filename("processed_data_result.root");
    TFile* f = new TFile(filename.Data(),"read");
    TTree* t = (TTree*)f->Get("caenv1730dump/events");

    // now prepare a output ROOT file
    //std::string outputfilename = filename.substr(0,filename.length()-5) + "_out.root";
    TString outputfilename=filename(0,filename.Length()-5);
    std::cout<<"*test "<<outputfilename<<std::endl;
    outputfilename += "_out_v3.root";
    std::cout<<"write a root file "<<outputfilename<<std::endl;
    TString outputtextfile = "all_pulse_amp_result.txt"; //filename(0,filename.Length()-5) + "_out_v3.txt";
    TFile* frootout = new TFile(outputfilename.Data(), "recreate");
    std::fstream fout(outputtextfile.Data(),ios::out | ios::app);
    fout<<"run file " << filename.Data() <<endl;
    //exit(0);

    std::vector<std::vector<uint16_t> > *data=0; //unsigned short
    t->SetBranchAddress("fWvfmsVec", &data);

    t->GetEntry(3);
    
    const int n = 16;// (*data).size();
    const int n1 = (*data)[0].size();
    int Nevts = t->GetEntries();

    const int n0 = (*data).size()/n;
    std::cout<<"number of fragements (boards): "<<n0<<std::endl;
    std::cout<<"number of channels per board : "<<n<<std::endl;
    std::cout<<"number of samples "<<n1<<std::endl;
    std::cout<<"number of triggers "<< Nevts<<std::endl;

    TH1F* h[n];// = new TH1F("waveform","",n1+1,0,n1+1);
    char name[100];
    /*
    for(int i=0;i<n;i++){
        sprintf(name,"ch%d",i);
        h[i] = new TH1F(name, name, n1, 0, n1);
        sprintf(name,"ch%d_baseline",i);
        //if(i==0){
        hBaseline[i] = new TH1F(name, name, 17000, 0, 17000);
        hBaseline[i]->SetXTitle("ADC counts");
        hBaseline[i]->SetYTitle("Counts");
        //}
    }
    */

    // preparing histograms
    TH1F* hBaseMean[n0][n];
    TH1F* hBaseWidth[n0][n];
    TH1F* hPulseAmp[n0][n];
    TH1F* hPulsePeakTimeBin[n0][n];
    TH1F* hPulseStartTimeBin[n0][n];
    TH1F* hPulseStartTimeBinFromFit[n0][n];
    TH1F* hPulseWidth[n0][n];
    TH1F* hPulseCharge[n0][n];
    TH1F* hPulseRiseTime[n0][n];
    TH1F* hPulseFallTime[n0][n];
    TH1F* hTimeDiffWrsC0[n0][n-1]; // pulse time difference with respect to the first ch
    TH2F* hPulseTimeCorr[n0][n-1];//look at pulse time colleration
    TH2F* hTimeDiffVsChargeDiff[n0][n-1]; // correlation of time diff and charge diff

    TH2F* hPulseChargeVsRiseTime[n0][n];
    TH2F* hPulseChargeVsFallTime[n0][n];

    TDirectory* dirWaveform[n0];
    TDirectory* dirPulseAmp[n0];

    for(int k=0; k<n0; k++){
        sprintf(name, "waveforms_board%d",k);
        dirWaveform[k] = frootout->mkdir(name);
        sprintf(name, "pulseAmp_board%d",k);
        dirPulseAmp[k] = frootout->mkdir(name);
        for(int i=0; i<n; i++){ // I don't need the last channel
            sprintf(name,"hBaselineMean_board%d_ch%d",k,i);
            hBaseMean[k][i] = new TH1F(name,name, 400,14500,15100);
            hBaseMean[k][i]->SetXTitle("Waveform baseline mean (ADC)");
            hBaseMean[k][i]->SetYTitle("Counts");
    
            sprintf(name,"hBaselineWidth_board%d_ch%d",k,i);
            hBaseWidth[k][i] = new TH1F(name,name, 100,0,20);
            hBaseWidth[k][i]->SetXTitle("Waveform baseline width (ADC)");
            hBaseWidth[k][i]->SetYTitle("Counts");
    
            sprintf(name,"hPulseAmp_board%d_ch%d",k,i);
            hPulseAmp[k][i] = new TH1F(name,name, 800,0,1600);
            //else hPulseAmp[k][i] = new TH1F(name,name, 1800,0,18000);
            hPulseAmp[k][i]->SetXTitle("Pulse amplitude (mV)");
            hPulseAmp[k][i]->SetYTitle("Counts");
    
            sprintf(name,"hPulsePeakTimeBin_board%d_ch%d",k,i);
            hPulsePeakTimeBin[k][i] = new TH1F(name,name, 500,400,500);
            hPulsePeakTimeBin[k][i]->SetXTitle("Time bin at pulse peak (ns)");
            hPulsePeakTimeBin[k][i]->SetYTitle("Counts");
    
            sprintf(name,"hPulseStartTime_board%d_ch%d",k,i);
            hPulseStartTimeBin[k][i] = new TH1F(name,name, 500,400,500);
            hPulseStartTimeBin[k][i]->SetXTitle("Time bin at pulse start (ns)");
            hPulseStartTimeBin[k][i]->SetYTitle("Counts");
    
            sprintf(name,"hPulseSTFromFit_board%d_ch%d",k,i);
            hPulseStartTimeBinFromFit[k][i] = new TH1F(name,name, 500,400,500);
            hPulseStartTimeBinFromFit[k][i]->SetXTitle("Time bin at pulse start from fitting the pulse (ns)");
            hPulseStartTimeBinFromFit[k][i]->SetYTitle("Counts");
    
            sprintf(name,"hPulseWidth_board%d_ch%d",k,i);
            hPulseWidth[k][i] = new TH1F(name,name, 300,0,300);
            hPulseWidth[k][i]->SetXTitle("Pulse width (ns)");
            hPulseWidth[k][i]->SetYTitle("Counts");
    
            sprintf(name,"hPulseCharge_board%d_ch%d", k,i);
            hPulseCharge[k][i] = new TH1F(name,name, 200,0,1000);
            hPulseCharge[k][i]->SetXTitle("Pulse charge (pC)");
            hPulseCharge[k][i]->SetYTitle("Counts");

            sprintf(name,"hPulseRiseTime_board%d_ch%d",k,i);
            hPulseRiseTime[k][i] = new TH1F(name,name, 100,0,10);
            hPulseRiseTime[k][i]->SetXTitle("Pulse rise time from fit (ns)");
            hPulseRiseTime[k][i]->SetYTitle("Counts");

            sprintf(name,"hPulseChargeVsRiseTime_board%d_ch%d",k,i);
            hPulseChargeVsRiseTime[k][i] = new TH2F(name,name, 100,0,10, 200, 0, 100);
            hPulseChargeVsRiseTime[k][i]->SetXTitle("Pulse rise time from fit (ns)");
            hPulseChargeVsRiseTime[k][i]->SetYTitle("Pulse charge (pc)");

            sprintf(name,"hPulseFallTime_board%d_ch%d",k,i);
            hPulseFallTime[k][i] = new TH1F(name,name, 200,0,20);
            hPulseFallTime[k][i]->SetXTitle("Pulse fall time from fit (ns)");
            hPulseFallTime[k][i]->SetYTitle("Counts");

            sprintf(name,"hPulseChargeVsFallTime_board%d_ch%d",k,i);
            hPulseChargeVsFallTime[k][i] = new TH2F(name,name, 200,0,20, 200, 0, 100);
            hPulseChargeVsFallTime[k][i]->SetXTitle("Pulse fall time from fit (ns)");
            hPulseChargeVsFallTime[k][i]->SetYTitle("Pulse charge (pc)");

            if(i>0){
                sprintf(name,"hPTdiff_bd%d_ch%d_toC0",k,i);
                hTimeDiffWrsC0[k][i-1] = new TH1F(name,name,200,-10,10);
                sprintf(name,"Pulse time diff b/t ch %d and ch0 in board %d", i,k);
                hTimeDiffWrsC0[k][i-1]->SetTitle(name);
                hTimeDiffWrsC0[k][i-1]->SetXTitle("Pulse time different (ns)");
                hTimeDiffWrsC0[k][i-1]->SetYTitle("Counts");

                sprintf(name,"hPTdiffVsQdiff_bd%d_ch%d_toC0",k,i);
                hTimeDiffVsChargeDiff[k][i-1] = new TH2F(name,name,200,-10,10, 40,-20,20);
                sprintf(name,"Pulse charge diff vs time diff b/t ch %d and ch0 in board %d", i,k);
                hTimeDiffVsChargeDiff[k][i-1]->SetTitle(name);
                hTimeDiffVsChargeDiff[k][i-1]->SetXTitle("Pulse time different (ns)");
                hTimeDiffVsChargeDiff[k][i-1]->SetYTitle("Charge difference (pC)");

                sprintf(name,"hPulseTimeCorr_bd%d_ch%d_toC0",k,i);
                hPulseTimeCorr[k][i-1] = new TH2F(name,name,500,400,500,500,400,500);
                sprintf(name,"Pulse time in Ch %d",i);
                hPulseTimeCorr[k][i-1]->SetXTitle("Pulse time in Ch0 (ns)");
                hPulseTimeCorr[k][i-1]->SetYTitle(name);
            }
        }
    } // end preparing histograms

    // temperature 

    // define the key variables
    double baseline=0, baselinew=0;
    double pulseamp=0, pulsepeakt=0;
    double pulseST=0, pulseST1=0, pulsewidth=0;
    double pulsecharge=0;
    double Vconverter = 0.122;// mV/ADC
    double Tconverter = 2.0; // ns per time bin;
    double Qconverter = 0.02; // pC,  = 1 mV* 1 ns / 50 Ohm
    std::vector<double> onewaveform; // its size should always be n1, the number of samples
    onewaveform.resize(n1);
    std::vector<double> subwaveform; // for baseline calcultion
    int N_base = 200;
    int Nsigma = 5;
    subwaveform.resize(N_base); // use the last two 200 samples in this analysis;
    // temp variables
    double asampleV=0, asampleV1=0;
    double pulsetimefitc0=0, pulsetimefit=0;
    double pulsechargec0=0; 
    TH1F* hPulse = new TH1F("pulse","",n1, 0, n1*Tconverter);
    TF1* myfitfunc = new TF1("myfitfunc",pulsefunc,400,1000,5); // function to fit the pulse
    myfitfunc->SetParNames("Amp","#tau_{s}","#tau", "t0", "Offset");

    TF1* myfitfunc2 = new TF1("myfitfunc2",pulsefunc2,400,1000,4); // function to fit the pulse
    myfitfunc2->SetParNames("Amp","#tau", "t0", "Offset");

    TF1* myfitfunc1 = new TF1("myfitfunc1",pulsefunc1,400,1000,4); // another function, using Daya Bay experience
    myfitfunc1->SetParNames("Amplitude", "shapeT", "widthT", "pulseT");

    for(int m=0; m<Nevts; m++){ // loop through triggers Nevts
        t->GetEntry(m);
        if(m % 100 ==0 ) std::cout<<"evt "<<m<<std::endl;
        for(int k=0; k<n0; k++){ // loop through digitizer boards
            //std::cout<<"    board "<<k<<std::endl;
            pulsetimefitc0=0; pulsechargec0=0;
            for(int i=0; i<n; i++){ // loop through channels
                //std::cout<<"           ch "<<i<<endl;
                for(int j=0; j<n1; j++){ // loop through waveform samples
                   //std::cout<<j<<"\t"<< 1.0* ( (*data)[i+n*k][j] ) <<std::endl;
                   // get one waveform
                   onewaveform[j] = 1.0* ( (*data)[i+n*k][j] );
                   if(j>= (n1-N_base) ){
                       subwaveform[j-n1+N_base] = 1.0* ( (*data)[i+n*k][j] );
                       baseline += subwaveform[j-n1+N_base] ;
                   }
                } // end waveform sample loop
                //**** baseline ****//
                baseline /= N_base; // this is the average of baseline
                //std::cout<<baseline<<std::endl;
                //now calculate baseline width (standard deviation)
                for(std::vector<double>::iterator it = subwaveform.begin(); it!=subwaveform.end(); ++it){
                    baselinew += TMath::Power( *it - baseline, 2);
                }
                baselinew = TMath::Sqrt( baselinew/(N_base-1) );
                //std::cout<<baselinew<<std::endl;
                hBaseMean[k][i]->Fill(baseline);
                hBaseWidth[k][i]->Fill(baselinew);

                //**** pulse amplitude ****//
                pulseamp = onewaveform[0];
                for(int j=1; j<n1; j++){
                    if(onewaveform[j]<pulseamp){ pulseamp = onewaveform[j] ; pulsepeakt = j+1;} 
                }
                hPulseAmp[k][i]->Fill( (baseline - pulseamp)*Vconverter );
                hPulsePeakTimeBin[k][i]->Fill( pulsepeakt*Tconverter );
               // std::cout<< baseline - pulseamp <<std::endl;

                if( (pulseamp >= baseline - Nsigma*baselinew) ){
                    std::cout<<"no pulse in evt "<<m<<" board "<< k << " ch "<<i<<std::endl;
                    baseline=0; baselinew=0; pulseST=0; pulsewidth=0; pulsecharge=0; pulseamp=0; pulsepeakt=0;
                    pulsetimefitc0=0; pulsechargec0=0;
                    continue;
                }

                //**** pulse start time bin ****//
                for(int j=pulsepeakt; j>=0; j--){
                    if(onewaveform[j]< (baseline - baselinew*Nsigma)) continue;
                    else{
                      pulseST = j-2;
                      break;
                    }
                }
                //**** pulse end time bin and pulse width ****//
                for(int j=pulsepeakt; j<n1; j++){
                    if(onewaveform[j]< (baseline - baselinew*Nsigma)) continue;
                    else{
                      pulsewidth = j+3 - pulseST;
                      break;
                    }
                }
                hPulseStartTimeBin[k][i]->Fill(pulseST *  Tconverter );
                hPulseWidth[k][i]->Fill(pulsewidth * Tconverter );

                //**** pulse charge ****//
                for(int j=pulseST; j<pulsewidth+pulseST; j++){
                    pulsecharge += baseline - onewaveform[j] ;
                }
                pulsecharge *= Qconverter* Vconverter * Tconverter ;
                hPulseCharge[k][i]->Fill(pulsecharge);
                //std::cout<<pulsecharge<<std::endl;

                // do the priliminary pulse fit
                sprintf(name,"hPulse_evt%d_board%d_ch%d",m,k,i);
                hPulse->SetName(name); hPulse->SetTitle(name);
                for(int j=0; j<n1; j++){
                    hPulse->SetBinContent(j+1, 1.0*(onewaveform[j]) ); // -baseline
                }
                myfitfunc->SetParameters(5000, 3.6, 8.2, pulseST*Tconverter, baseline);
                //myfitfunc->SetParLimits(1,2,6);
                //myfitfunc->SetParLimits(2,6,15);
                //hPulse->Fit("myfitfunc","QR","",pulsepeakt*Tconverter-30,pulsepeakt*Tconverter+20);
                hPulse->Fit("myfitfunc","QR","",180,400);
                hPulse->GetXaxis()->SetRangeUser(180,400);

                //myfitfunc2->SetParameters(5000, 5, pulseST*Tconverter, 0);
                //hPulse->Fit("myfitfunc2","QR","",pulsepeakt*Tconverter-30,pulsepeakt*Tconverter+20);
                
                //myfitfunc1->SetParameters(-300, 6, 1, pulseST*Tconverter);
                //hPulse->Fit("myfitfunc1","QR","",pulsepeakt*Tconverter-30,pulsepeakt*Tconverter+20);
                //hPulse->Fit("gaus","QR","",pulsepeakt*Tconverter-15,pulsepeakt*Tconverter+20);

                pulsetimefit = myfitfunc->GetParameter(3); //pulseST*Tconverter;//
                //pulsetimefit = hPulse->GetFunction("gaus")->GetParameter(1); //pulseST*Tconverter;//
                //pulsetimefit = myfitfunc2->GetParameter(2);
                //pulsetimefit = myfitfunc1->GetParameter(3);
                hPulseStartTimeBinFromFit[k][i]->Fill( pulsetimefit );

                hPulseRiseTime[k][i]->Fill( myfitfunc->GetParameter(1) );
                hPulseFallTime[k][i]->Fill( myfitfunc->GetParameter(2) );

                //fout<<m<<"\t"<<i<<"\t"<<myfitfunc->GetParameter(1)<<std::endl;

                //hPulseChargeVsRiseTime[k][i]->Fill( myfitfunc->GetParameter(1), pulsecharge );
                //hPulseChargeVsFallTime[k][i]->Fill( myfitfunc->GetParameter(2), pulsecharge );

                // here I do another linear fit to the rise of the pulse
                //double tstart = pulsetimefit;
                //double tend = pulsepeakt*Tconverter;
                //if(tend>500) tend=448;
                //std::cout<<"Before fit, check range "<<pulsepeakt*Tconverter-6<<std::endl;
                //hPulse->Fit("pol1","QR","",tend-6, tend);
                //double p0 = hPulse->GetFunction("pol1")->GetParameter(0);
                //double p1 = hPulse->GetFunction("pol1")->GetParameter(1);
                //pulsetimefit = (baseline - p0)/p1;
                //frootout->cd();
                dirWaveform[k]->cd(); 
                hPulse->Write();

                //pulse time difference with repect to Ch0 in each board
                if(i==0){
                   pulsetimefitc0= pulsetimefit;//myfitfunc->GetParameter(3); //pulseST*Tconverter;//// pulsetimefitc0 needs to be 0 before entering the loop
                   pulsechargec0 = pulsecharge;
                 }
                else{
                    hTimeDiffWrsC0[k][i-1]->Fill(pulsetimefit - pulsetimefitc0);
                    hPulseTimeCorr[k][i-1]->Fill(pulsetimefitc0, pulsetimefit);
                    hTimeDiffVsChargeDiff[k][i-1]->Fill(pulsetimefit - pulsetimefitc0, pulsecharge - pulsechargec0);
                    //std::cout<<"pulse time "<<pulsetimefit<<"\t"<<pulsetimefitc0<<std::endl;
                }

                for(int j=0; j<n1; j++){
                    hPulse->SetBinContent(j+1, 0);
                }

                //std::cout<<"test"<<std::endl;
                // clear the variables
                baseline=0; baselinew=0; pulseST=0; pulsewidth=0; pulsecharge=0; pulseamp=0; pulsepeakt=0;
                //pulsetimefitc0=0;
                pulsetimefit=0;
            } // end channel loop
        } // end digitizer board loop
    } // end trigger number loop
    
    //f->Close();
    //fout.close();
    //TCanvas* canvas=new TCanvas();
    //hBaseMean[0][0]->Draw();
    
    // save file
    frootout->cd();

    for(int k=0; k<n0; k++){
        for(int i=0; i<n; i++){ // I don't need the last channel
            hBaseMean[k][i]->Write();
            hBaseWidth[k][i]->Write();

            //hPulseAmp[k][i]->Write();

            hPulsePeakTimeBin[k][i]->Write();
            hPulseStartTimeBin[k][i]->Write();
            hPulseStartTimeBinFromFit[k][i]->Write();
            hPulseWidth[k][i]->Write();
            hPulseCharge[k][i]->Write();
            hPulseRiseTime[k][i]->Write();
            hPulseFallTime[k][i]->Write();
            hPulseChargeVsRiseTime[k][i]->Write();
            hPulseChargeVsFallTime[k][i]->Write();
            if(i>0){
                hTimeDiffWrsC0[k][i-1]->Write();
                hPulseTimeCorr[k][i-1]->Write();
                hTimeDiffVsChargeDiff[k][i-1]->Write();
            }
        }
    }



    TCanvas* cAmpDis[n0];//
    TPaveStats *st[n0];// = (TPaveStats*)h->FindObject("stats")
    //int myindex=3;
    for(int k=0; k<n0; k++){
        sprintf(name,"ampDist_board%d",k);
        cAmpDis[k] = new TCanvas(name,name,1200,960);
        cAmpDis[k]->Divide(4,4);
        dirPulseAmp[k]->cd();
        for(int i=0; i<n; i++){
            hPulseAmp[k][i]->Write();
            cAmpDis[k]->cd(i+1);
            hPulseAmp[k][i]->GetXaxis()->SetLabelSize(0.06);
            hPulseAmp[k][i]->GetYaxis()->SetLabelSize(0.05);
            hPulseAmp[k][i]->Draw();
            cAmpDis[k]->Update();
            st[k] = (TPaveStats*)hPulseAmp[k][i]->FindObject("stats");
            st[k]->SetX1NDC(0.6); //new x start position
            st[k]->SetX2NDC(0.9); //new x end position
            st[k]->SetY1NDC(0.4);
            cAmpDis[k]->Update();
            fout << i <<"\t" << hPulseAmp[k][i]->GetMean() <<"\t" << hPulseAmp[k][i]->GetRMS() <<endl;
        }
        cAmpDis[k]->Write();
        //cAmpDis[k]->Close();
    }
    
    fout.close();
    //frootout->Close();

    return 0;
    //exit(0);
}
