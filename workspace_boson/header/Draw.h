TFile   * file;
TH1D    * hist, *hist_W, *hist_Z;
TCanvas * can;
TLatex  * latex = new TLatex(0,0,"");
int n_column=5, row=0, num_mass=15;
double offset=10, length=275, width=200;
double mean, sigma, MASS, ENERGY = 320;// MC boson multiplicity, Heavy quark mass, Goldstone boson energy in GeV 
double expected_mean, expected_sigma;// Predicted boson multiplicity
//Below two lines are not allowed
//double prediction_mean[num_mass]  = {3.75, 4.38, 5.00, 5.62, 6.25, 6.88, 7.50, 8.12, 8.75, 9.38, 10.00, 10.63, 11.25, 11.88, 12.50};
//double prediction_sigma[num_mass] = {0.97, 1.05, 1.12, 1.19, 1.25, 1.31, 1.37, 1.43, 1.48, 1.53, 1.58, 1.63, 1.68, 1.72, 1.77};

void SetLatex(){
    latex->SetNDC();
	latex->SetTextFont(43);
    latex->SetTextSize(16);
}

void PlotCombine(bool DOPLOTTING, const char* Variable){
    printf("Mass &   Prediction   &    MC Bosons   &    W Bosons    &   Z Bosons\\\\ \n");
    for(int i=6; i<21; i++){
        row=(i-6)/n_column;
        MASS = (double)i/10.;//Mass in TeV
        if(i==10 || i==20) file = new TFile(Form("output/result_run_fireball_%dTeV.root",i/10),"READ");
        else               file = new TFile(Form("output/result_run_fireball_%.1fTeV.root",MASS),"READ");
        if(DOPLOTTING)     can  = new TCanvas(Form("c%d",i),"",offset+length*((i-6)%n_column),offset+(width+25)*row,length,width);
        //===== Input Hist =====/
        hist   = (TH1D*) file->Get(Form("hist_VB_%s",Variable));
        hist_W = (TH1D*) file->Get(Form("hist_W_%s" ,Variable));
        hist_Z = (TH1D*) file->Get(Form("hist_Z_%s" ,Variable));
        //===== Hist Setting =====/
        hist   -> GetXaxis() -> SetTitle(GetXtitle(Variable));
        hist_W -> GetXaxis() -> SetTitle(GetXtitle(Variable));
        hist_Z -> GetXaxis() -> SetTitle(GetXtitle(Variable));
        hist   -> SetLineColor(kBlack);
        hist_W -> SetLineColor(kBlue);
        hist_Z -> SetLineColor(kRed);
        //===== Drawing =====/
        if(DOPLOTTING) hist   -> Draw();
        if(DOPLOTTING) hist_W -> Draw("same");
        if(DOPLOTTING) hist_Z -> Draw("same");
        //===== DrawLatex =====/
        mean = hist->GetMean();  sigma = hist->GetStdDev();
        if(DOPLOTTING) latex -> DrawLatex(0.60, 0.8, Form("%.1f TeV",MASS));
        if(DOPLOTTING) latex -> DrawLatex(0.60, 0.7, Form("Mean:%.2f",mean));
        if(DOPLOTTING) latex -> DrawLatex(0.60, 0.6, Form("RMS: %.2f",sigma));

        expected_mean  = 2*MASS*1000/ENERGY;
        expected_sigma = sqrt(expected_mean)/2.;
        printf("%4.0f & %5.2f \\pm %.2f & %5.2f \\pm %.2f & ", MASS*1000, expected_mean, expected_sigma, mean, sigma);
        mean = hist_W->GetMean();  sigma = hist_W->GetStdDev(); printf("%5.2f \\pm %.2f & ", mean, sigma);
        mean = hist_Z->GetMean();  sigma = hist_Z->GetStdDev(); printf("%5.2f \\pm %.2f ", mean, sigma);
        printf("\\\\ \n");
        //printf("%5.2f \\pm %.2f & %5.2f \\pm %.2f \\\\ \n", MASS*1000, expected_mean, expected_sigma, mean, sigma);
    }
}

void RegisterMeanAndSignal(TH1D *hist, double &mean, double &sigma){
    mean = hist->GetMean();  sigma = hist->GetStdDev();
}


void PlotIndividualy(const char* histName){
    for(int i=6; i<21; i++){
        row=(i-6)/n_column;
        MASS = (double)i/10.;
        if(i==10 || i==20) file = new TFile(Form("output/result_run_fireball_%dTeV.root",i/10),"READ");
        else               file = new TFile(Form("output/result_run_fireball_%.1fTeV.root",MASS),"READ");
        can   = new TCanvas(Form("c%d",i),"",offset+length*((i-6)%n_column),offset+(width+25)*row,length,width);
        hist  = (TH1D*) file->Get(histName);
        mean  = hist->GetMean();
        sigma = hist->GetStdDev();
        printf("mean = %.2f, sigma = %.2f\n", mean, sigma);
        hist  -> GetXaxis() -> SetTitle(GetXtitle(histName));
        hist  -> Draw();
        latex -> DrawLatex(0.60, 0.8, Form("%.1f TeV",MASS));
        latex -> DrawLatex(0.60, 0.7, Form("Mean:%.2f",mean));
        latex -> DrawLatex(0.60, 0.6, Form("RMS: %.2f",sigma));
    }
}

const char* GetXtitle(const char* histName){
    if(histName == "Multiplicity" || histName == "hist_VB_Multiplicity") return "multiplicity";
    if(histName == "Momentum"     || histName == "hist_VB_Momentum"    ) return "[GeV]";
    if(histName == "Energy"       || histName == "hist_VB_Energy"      ) return "[GeV]";
    if(histName == "Multiplicity" || histName == "hist_W_Multiplicity" ) return "multiplicity";
    if(histName == "Momentum"     || histName == "hist_W_Momentum"     ) return "[GeV]";
    if(histName == "Energy"       || histName == "hist_W_Energy"       ) return "[GeV]";
    if(histName == "Multiplicity" || histName == "hist_Z_Multiplicity" ) return "multiplicity";
    if(histName == "Momentum"     || histName == "hist_Z_Momentum"     ) return "[GeV]";
    if(histName == "Energy"       || histName == "hist_Z_Energy"       ) return "[GeV]";
}
