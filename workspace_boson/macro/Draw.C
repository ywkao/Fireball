#include "header/Draw.h"
bool VB           = false, W        = false, Z      = false;
bool Multiplicity = false, Momentum = false, Energy = false;

bool PLOT=true;

void Draw(){
    gStyle->SetOptStat(0);
    SetLatex();
    PlotCombine(PLOT,"Multiplicity");
    if( VB && Multiplicity ) PlotIndividualy ("hist_VB_Multiplicity");
    if( VB && Momentum )     PlotIndividualy ("hist_VB_Momentum"    );
    if( VB && Energy )       PlotIndividualy ("hist_VB_Energy"      );
    if( W  && Multiplicity ) PlotIndividualy ("hist_W_Multiplicity" );
    if( W  && Momentum )     PlotIndividualy ("hist_W_Momentum"     );
    if( W  && Energy )       PlotIndividualy ("hist_W_Energy"       );
    if( Z  && Multiplicity ) PlotIndividualy ("hist_Z_Multiplicity" );
    if( Z  && Momentum )     PlotIndividualy ("hist_Z_Momentum"     );
    if( Z  && Energy )       PlotIndividualy ("hist_Z_Energy"       );
}

    /*
    can = new TCanvas("c1","",10,10,300,225);
    file = new TFile("output/result_run_fireball_2TeV.root","READ");
    hist = (TH1D*) file->Get("hist_VB_Multiplicity");
    mean  = hist->GetMean();
    sigma = hist->GetStdDev();
    printf("mean = %.2f, sigma = %.2f\n", mean, sigma);
    hist->Draw();

    can = new TCanvas("c2","",310,10,300,225);
    file = new TFile("output/result_run_fireball_1TeV.root","READ");
    hist = (TH1D*) file->Get("hist_VB_Energy");
    mean  = hist->GetMean();
    sigma = hist->GetStdDev();
    printf("mean = %.2f, sigma = %.2f\n", mean, sigma);
    hist->Draw();
    */
