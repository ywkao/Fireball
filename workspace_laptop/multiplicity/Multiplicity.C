#include <TF1.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>

void Multiplicity(){
    TH1D *hist = new TH1D("hist",";# of Goldstone bosons; Probability (%)",50,0,25);
    TF1 *f1 = new TF1("f1","31.9*exp(-(x-6.25)^2/3.13)",0,25);
    TF1 *f2 = new TF1("f2","22.6*exp(-(x-12.5)^2/6.25)",0,25);
    f1->SetLineColor(kBlue);
    f2->SetLineColor(kRed);

    gStyle->SetOptStat(0);
    hist->SetMaximum(50);
    hist->Draw();
    f1->Draw("same");
    f2->Draw("same");

    TLegend *legend = new TLegend(0.7,0.6,0.85,0.85);
    TLatex *text = new TLatex(0,0,"");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(43);
    legend->SetTextSize(25);
    legend->AddEntry(f1,"1TeV","l");
    legend->AddEntry(f2,"2TeV","l");
    legend->Draw("same");
	text->SetTextFont(43);
	text->SetTextSize(18);
    text->DrawLatex(2,42,"P(n_{G})= #frac{1}{#sqrt{2#pi}#sigma}exp#(){#minus(n_{G}#minus<n_{G}>)^{2}/2#sigma^{2}}");

    const int n=18;
    Double_t Num[n] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
    Double_t prob1[n] = {0,0.1,1,6,19,31,27,12,3,0.4,0,0,0,0,0,0,0};
    Double_t prob2[n] = {0,0,0,0,0,0,0.2,0.9,3,8,16,22,22,16,8,3,0.9,0.2};
    TGraph *graph1 = new TGraph(n, Num, prob1);
    TGraph *graph2 = new TGraph(n, Num, prob2);
    graph1 -> SetMarkerStyle(20);
    graph1 -> SetMarkerColor(kBlue);
    graph1 -> SetLineColor(kBlue);
    graph2 -> SetMarkerStyle(20);
    graph2 -> SetMarkerColor(kRed);
    graph2 -> SetLineColor(kRed);
    graph1 -> Draw("p,same");
    graph2 -> Draw("p,same");

    const int NUM=20;
    double Energy = 320;//GeV for Goldstone boson
    //double MASS[6] = {1,1.2,1.4,1.6,1.8,2.0};
    double MASS, Mean, Sigma;
    for(int i=0; i<18; i++){
        MASS = 0.3+0.1*(double)i;
        Mean = 2*MASS*1000/Energy; Sigma = sqrt(Mean)/2.;
        for(int j=0; j<NUM; j++){
            if(j==0)        printf ("MASS = %2.1f TeV, prob = [%3.1f, ",MASS,Probability(Mean,Sigma,j));
            else if(j!=NUM-1) printf ("%3.1f, ",Probability(Mean,Sigma,j));
            else            printf ("%3.1f]\n",Probability(Mean,Sigma,j));
        }
    }


}

double Probability(double mean, double sigma, double n){
    return 1/(sqrt(2*TMath::Pi())*sigma)*exp( -pow( (n-mean)/(sqrt(2)*sigma) ,2) )*100;
}
