#include <stdio.h>
#include <TCanvas.h>
#include "TGaxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <TH2D.h>
#include "TLatex.h"
#include "TLegend.h"
#include <TPolyLine.h>
#include <TPad.h>
#include "signalStrength.h"
using namespace std;

const int UP=1, LOW=0;
const double UPPER_BOUND=10, LOWER_BOUND=1e-04;
//const double UPPER_BOUND=1.0e-01, LOWER_BOUND=-5e-03;
//const double UPPER_BOUND=2.5e-02, LOWER_BOUND=-1e-03;
//const double UPPER_BOUND=100, LOWER_BOUND=1e-04;

void ExclusionXsec(){
	TCanvas *can = new TCanvas("can","",800,600);
	can -> SetLogy();
	can -> SetLeftMargin(0.12);
	can -> SetBottomMargin(0.12);

    //# MASS, Xsec, Err_Xsec, signal strength which are recorded in the header file
    
    //cross section
    double Xmean[NUM], Xs1[2][NUM], Xs2[2][NUM], Xobs[NUM];
    for(int i=0; i<NUM; i++){
        Xmean[i]=Xsec[i]*mean[i]*K;
        Xobs[i] =Xsec[i]*obs[i]*K;
        for(int j=0; j<2; j++){//0:lower, 1:upper
            Xs1[j][i] = (s1[j][i]-mean[i])/mean[i]*Xmean[i]+Xmean[i];
            Xs2[j][i] = (s2[j][i]-mean[i])/mean[i]*Xmean[i]+Xmean[i];
        }
    }

    //TGraphErrors *graph_theory = new TGraphErrors(NUM, MASS, Xmeasured, Err_MASS, Err_Xmeasured);
    TH2D   *frame        = new TH2D("frame",";M_{Q} [TeV];95\%  CL  #sigma(Q#bar{Q} #rightarrow nV) [pb]",10, MASS_min, MASS_max, 50, LOWER_BOUND, UPPER_BOUND);
    TGraph *graph        = new TGraph(NUM, MASS, Xmean);
    TGraph *graph_obs    = new TGraph(NUM, MASS, Xobs);
    TGraph *graph_theory = new TGraphErrors(NUM, MASS, Xmeasured);
    TPolyLine *P1        = new TPolyLine(NUM*2);
    TPolyLine *P2        = new TPolyLine(NUM*2);

    for(int i=0; i<NUM*2; i++){
        if(i<NUM){
            P1->SetNextPoint(MASS[i], Xs1[LOW][i]);
            P2->SetNextPoint(MASS[i], Xs2[LOW][i]);
        } else{
            P1->SetNextPoint(MASS[2*NUM-1-i], Xs1[UP][2*NUM-1-i]);
            P2->SetNextPoint(MASS[2*NUM-1-i], Xs2[UP][2*NUM-1-i]);
        }
    }

	frame -> SetStats(0);
    frame -> GetYaxis() -> SetTitleOffset(1.5);
    P1->SetFillColor(kGreen);
    P2->SetFillColor(kYellow);
    //P1->SetLineColor(kGreen);
    //P2->SetLineColor(kYellow);
    P1->SetFillStyle(1001);
    P2->SetFillStyle(1001);
	graph -> SetLineStyle(2);
	graph -> SetLineWidth(2);
	graph_obs -> SetLineStyle(1);
	graph_obs -> SetLineWidth(2);
    graph_theory -> SetMarkerStyle(20);
    graph_theory -> SetMarkerColor(1);

	frame -> Draw();
    P2->Draw("f,same");
    P1->Draw("f,same");
    graph->Draw("pl,same");
    graph_obs->Draw("pl,same");
    graph_theory->Draw("pl,same");

	TLegend *legend = new TLegend(0.62,0.66,0.87,0.87);
	legend-> AddEntry(graph_theory," Theoretical value","lp");
	legend-> AddEntry(graph_obs," Observed Limit","l");
	legend-> AddEntry(graph," Expected Limit","l");
	legend-> AddEntry(P1," #pm 1 #sigma Expected","f");
	legend-> AddEntry(P2," #pm 2 #sigma Expected","f");
	legend-> SetLineColor(0);
	legend-> Draw("same");

	TLatex *text = new TLatex(0,0,"");
    text->SetTextFont(43);
    text->SetTextSize(18);
	text->DrawLatex(0.425,UPPER_BOUND*1.2,"MC Simulation from MadGraph5");
	text->DrawLatex(1.750,UPPER_BOUND*1.2,"L = 100 fb^{-1} (13TeV)");
    text->SetTextSize(28);
    //text->DrawLatex(1.200, 0.1, "#sigma = #frac{N-N_{B}}{#epsilon #times L}");
	//text->DrawLatex(1.700,0.0255,"#int L dt = 100 fb^{-1} (#sqrt{s}=13TeV)");
    
    text->DrawLatex(1,1,Form("N_{obs} = %d",Nobs));
    can->SaveAs("UpperLimitCrossSection.png");
}
