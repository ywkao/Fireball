void ana_MPV(){
    TFile *file = new TFile("output.root","READ");
    TH1D *hist;// = (TH1D*) file->Get("hist_1000");

    double MPV, MPVerrorLow, MPVerrorHigh;
    for(int i=0; i<11; i++){
        if(i!=0) continue;
        hist = (TH1D*) file->Get(Form("hist_%d",100*i+1000));
        hist->GetXaxis()->SetLimits(-0.1,0.3);
        hist->Rebin(1);
        FindMPV(hist, MPV, MPVerrorLow, MPVerrorHigh, 2, 1);
        PlotPosterior( 0, 0, "Cross Section [pb]", Form("yield/posteior_%d.png",100*i+1000), hist, MPV, MPVerrorLow, MPVerrorHigh);
    }
}


void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo, int nSigma){
    if(MPValgo==1){
        MPV=PosteriorDist->GetMean();
        MPVerrorLow=PosteriorDist->GetRMS();
        MPVerrorHigh=PosteriorDist->GetRMS();
    }
    
    if(MPValgo==2||MPValgo==3){
        int nBins = PosteriorDist->GetNbinsX();
        int maxbin_PosteriorDist = PosteriorDist->GetMaximumBin();
        double PosteriorDist_initial = PosteriorDist->GetBinCenter(maxbin_PosteriorDist);
        double err_PosteriorDist_initial=PosteriorDist->GetRMS();
        double PosteriorDist_par [3];
        
        TF1 *gauss;
        
        int nMaxFits=1;
        if(MPValgo==3) nMaxFits=20;
        for(int iFits=0;iFits<nMaxFits;iFits++){
            gauss = new TF1("f1", "gaus", PosteriorDist_initial-err_PosteriorDist_initial, PosteriorDist_initial+err_PosteriorDist_initial);
            gauss->SetParameters(PosteriorDist_initial,err_PosteriorDist_initial);
            PosteriorDist->Fit(gauss, "R");
            gauss->GetParameters(PosteriorDist_par);
            double ndof = 2*err_PosteriorDist_initial/PosteriorDist->GetBinWidth(1)-3;
            cout<<"chi2/ndf = "<<gauss->GetChisquare()/ndof<<endl;
            PosteriorDist_initial=PosteriorDist_par[1];
            err_PosteriorDist_initial=err_PosteriorDist_initial/2;
            if(gauss->GetChisquare()/ndof<5) break;
        }
        MPV=PosteriorDist_par[1];
    
        double OneSigmaCL;
        if(nSigma==1) OneSigmaCL=0.682689492137;
        if(nSigma==2) OneSigmaCL=0.954499736104;
        if(nSigma==3) OneSigmaCL=0.997300203937;
        double fullInt=PosteriorDist->Integral(1,nBins);
        cout<<(1-OneSigmaCL)/2.<<endl;
        
        for(int i = 1; i < nBins+1; i++){
            //cout<<i<<" "<<PosteriorDist->Integral(1,i)/fullInt<<endl;
            if(PosteriorDist->Integral(1,i)/fullInt > (1-OneSigmaCL)/2.) {MPVerrorLow=MPV-PosteriorDist->GetBinCenter(i-1); break;}
        }
        for(int i = 1; i < nBins+1; i++){
            //cout<<i<<" "<<PosteriorDist->Integral(nBins+1-i,nBins)/fullInt<<endl;
            if(PosteriorDist->Integral(nBins+1-i,nBins)/fullInt > (1-OneSigmaCL)/2.) {MPVerrorHigh=PosteriorDist->GetBinCenter(nBins-i)-MPV; break;}
        }
    }
    
    return;
}

void PlotPosterior(int ptBin, int rapBin, char xAxisTitle[200], char filename[200], TH1* histo, double MPV, double MPVerrLow, double MPVerrHigh){
    gStyle->SetPalette(1,0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.15);
    
    gStyle->SetTickLength(-0.02, "xyz");
    gStyle->SetLabelOffset(0.02, "x");
    gStyle->SetLabelOffset(0.02, "y");
    gStyle->SetTitleOffset(1.3, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleFillColor(kWhite);
    
    TLegend* plotLegend=new TLegend(0.7,0.7,0.95,0.9);
    plotLegend->SetFillColor(kWhite);
    plotLegend->SetTextFont(72);
    plotLegend->SetTextSize(0.03);
    plotLegend->SetBorderSize(1);
    
    TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",800,600);
    
    plotCanvas->SetFillColor(kWhite);
    plotCanvas->SetGrid();
    plotCanvas->GetFrame()->SetFillColor(kWhite);
    plotCanvas->GetFrame()->SetBorderSize(0);
    plotCanvas->SetRightMargin(0.05) ;
    
    histo->SetStats(kFALSE);
    histo->SetLineColor(kBlack);
    histo->SetYTitle("Posterior Probability");
    histo->SetXTitle(xAxisTitle);
    histo->GetYaxis()->SetTitleOffset(1.5);
    histo->Draw();
    
    
    maxbin_PosteriorDist = histo->GetMaximumBin();
    double maxval = histo->GetBinContent(maxbin_PosteriorDist);
    
    TLine* MeanLine = new TLine( histo->GetMean(), 0, histo->GetMean(), maxval );
    MeanLine->SetLineWidth( 1 );
    MeanLine->SetLineStyle( 1 );
    MeanLine->SetLineColor( kGreen+2 );
    MeanLine->Draw( "same" );
    TLine* MPVLine = new TLine( MPV, 0, MPV, maxval );
    MPVLine->SetLineWidth( 2.5 );
    MPVLine->SetLineStyle( 1 );
    MPVLine->SetLineColor( kRed );
    MPVLine->Draw( "same" );
    TLine* MPVerrLowLine = new TLine( MPV-MPVerrLow, 0, MPV-MPVerrLow, maxval );
    MPVerrLowLine->SetLineWidth( 2.5 );
    MPVerrLowLine->SetLineStyle( 2 );
    MPVerrLowLine->SetLineColor( kRed );
    MPVerrLowLine->Draw( "same" );
    TLine* MPVerrHighLine = new TLine( MPV+MPVerrHigh, 0, MPV+MPVerrHigh, maxval );
    MPVerrHighLine->SetLineWidth( 2.5 );
    MPVerrHighLine->SetLineStyle( 2 );
    MPVerrHighLine->SetLineColor( kRed );
    MPVerrHighLine->Draw( "same" );
    plotLegend->AddEntry(MeanLine,"Mean","l");
    plotLegend->AddEntry(MPVLine,"MPV","l");
    plotLegend->AddEntry(MPVerrLowLine,"1#sigma high/low","l");
    
    plotLegend->Draw("same");
    plotCanvas->SaveAs(filename);
    plotCanvas->Close();
    
    delete plotCanvas;
}
