
void process_test_step2_multiBoards(){

    TFile* f = new TFile("../data/data_dl2_run1077_23_20200205T183344_dl4.root","read");//processed_data_result.root","read");
    TTree* t = (TTree*)f->Get("caenv1730dump/events");

    std::vector<std::vector<uint16_t> > *data=0; //unsigned short
    t->SetBranchAddress("fWvfmsVec", &data);

    //std::cout<<"Tree entries (# of evts): \t"<<t->GetEntries()<<std::endl;

    t->GetEntry(10);

    int chN = 4;
    const int n = 16;// (*data).size();

    const int n0 = (*data).size()/n;
    std::cout<<"number of fragements (boards): "<<n0<<endl;
    std::cout<<"number of channels per board : "<<n<<endl;
    const int n1 = (*data)[0].size();
    std::cout<<"number of samples: "<<n1<<std::endl;

    TH1F* h[n0][n];// = new TH1F("waveform","",n1+1,0,n1+1);
    TH1F* hBaseline[n0][n];
    char name[100];
    for(int i=0;i<n0;i++){
        for(int j=0;j<n;j++){
            sprintf(name,"board%d_ch%d",i,j);
            h[i][j] = new TH1F(name, name, n1, 0, n1);
            sprintf(name,"board%d_ch%d_baseline",i,j);

            hBaseline[i][j] = new TH1F(name, name, 17000, 0, 17000);
            hBaseline[i][j]->SetXTitle("ADC counts");
            hBaseline[i][j]->SetYTitle("Counts");
        }
    }

    TH1F* h1 = new TH1F("dist","",100,15000,15100);
    //TGraph* g = new TGraph();
    for(int i=0; i<n0; i++){
        for(int j=0; j<n; j++){
           for(int k=0; k<n1; k++){
               h[i][j]->SetBinContent(k+1, 1.0*((*data)[j+n*i][k]));

               hBaseline[i][j]->Fill( 1.0*((*data)[j+n*i][k]) );
           }
        }
    }

    //f->Close();

    TCanvas * c[n0];
    int index=4;//index+1
    for(int i=index; i<index+1; i++){
        sprintf(name,"canvas_board_%d", i);
        c[i] = new TCanvas(name,name,1600,900);
        c[i]->Divide(4,n/4);
        for(int j=0; j<n; j++){
           c[i]->cd(j+1);
           gStyle->SetOptStat(0);
           h[i][j]->Draw();
           h[i][j]->GetXaxis()->SetLabelSize(0.06);
           h[i][j]->GetYaxis()->SetLabelSize(0.05);
        }
        c[i]->Modified();
        c[i]->Update();
    }

    //c->SaveAs("20191001_waveform_1.root");


/*
    fstream fout("baseline_measure20191001_1.txt",ios::out);
    TCanvas* c1 = new TCanvas("canvas1", "canvas1", 1600, 1200);
    c1->Divide(4,4);
    for(int i=0; i<16; i++){
        //if(i!=0) continue;
        c1->cd(i+1);
        gStyle->SetOptFit(1111);
        hBaseline[i]->Draw();
        hBaseline[i]->GetXaxis()->SetLabelSize(0.06);
        hBaseline[i]->GetYaxis()->SetLabelSize(0.05);
        double mean = hBaseline[i]->GetMean();
        double rms = hBaseline[i]->GetRMS();
        hBaseline[i]->Fit("gaus","","R",mean-3.0*rms, mean+3.0*rms);
        fout<<hBaseline[i]->GetFunction("gaus")->GetParameter(1)<<"\t"
            <<hBaseline[i]->GetFunction("gaus")->GetParError(1)<<"\t"
            <<hBaseline[i]->GetFunction("gaus")->GetParameter(2)<<"\t"
            <<hBaseline[i]->GetFunction("gaus")->GetParError(2)<<"\t"
            <<endl;
    }
    c1->SaveAs("20191001_waveform_baseline_1.root");

    cout<<"flag 1"<<endl;
*/
    //c->SaveAs("test.pdf");
    //c->SaveAs("test.root");
    //exit(0);
//    return 0;
}
