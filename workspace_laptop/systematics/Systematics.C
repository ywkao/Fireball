#include <stdio.h>
#include <TCanvas.h>
#include "TGaxis.h"
#include "TGraphErrors.h"
#include <TH2D.h>
#include "TLatex.h"
#include "TLegend.h"
#include <TPolyLine.h>
#include <TPad.h>
using namespace std;
const int NUM = 11;
const double UPPER_BOUND=1.0e-00, LOWER_BOUND=1e-04;

void Systematics(){
	TCanvas *can = new TCanvas("can","",800,600);
	can -> SetLogy();
	can -> SetLeftMargin(0.12);
	can -> SetBottomMargin(0.12);

    double MASS[NUM]          = {1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};//TeV
    double Err_MASS[NUM]      = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double Xmeasured[NUM]     = {5.48e-02, 2.82e-02, 1.50e-02, 8.24e-03, 4.62e-03, 2.63e-03, 1.53e-03, 8.97e-04, 5.32e-04, 3.18e-04, 1.91e-04};
    double Err_Xmeasured[NUM] = {9.12e-03, 5.96e-03, 4.09e-03, 3.12e-03, 2.37e-03, 2.01e-03, 1.84e-03, 1.66e-03, 1.49e-03, 1.43e-03, 1.40e-03};

    TH2D *frame = new TH2D("frame",";M_{Q} [TeV];#sigma(Q#bar{Q} #rightarrow nV) [pb]",10, 0.8, 2.2, 50, LOWER_BOUND, UPPER_BOUND);
	frame->SetStats(0);
    frame->GetYaxis()->SetTitleOffset(1.5);
	frame->Draw();

    TGraphErrors *graph_theory = new TGraphErrors(NUM, MASS, Xmeasured, Err_MASS, Err_Xmeasured);
    graph_theory->SetMarkerStyle(20);
    graph_theory->SetMarkerColor(4);
    graph_theory->SetLineColor(4);
    graph_theory->Draw("pl,same");

	TLegend *legend = new TLegend(0.62,0.66,0.87,0.87);
	legend->AddEntry(graph_theory," Theoretical value","lp");
	legend->SetLineColor(0);
	legend->Draw("same");

	TLatex *text = new TLatex(0,0,"");
    text->SetTextFont(43);
    text->SetTextSize(18);
	text->DrawLatex(0.825,UPPER_BOUND+0.1,"MC Simulation from MadGraph5");
	text->DrawLatex(1.850,UPPER_BOUND+0.1,"L = 100 fb^{-1} (13TeV)");
    text->SetTextSize(24);
    //text->DrawLatex(1.200, 0.1, "#sigma = #frac{N-N_{B}}{#epsilon #times L}");
	//text->DrawLatex(1.700,0.0255,"#int L dt = 100 fb^{-1} (#sqrt{s}=13TeV)");
    
    can->SaveAs("SystematicsCrossSection.png");
}
