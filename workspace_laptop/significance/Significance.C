#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TH2D.h>
#include <TMath.h>
#include <TPad.h>
#include <Math/DistFuncMathCore.h>
#include "Significance.h"

void Significance(){
    //=== P-Value ===//
	TCanvas *c1  = new TCanvas("c1","",10,10,800,600); c1 -> SetLogy();
    TH2D *frame1 = new TH2D("frame1",";M_{Q} [TeV]; Local p-value",10, NUM_min, NUM_max, 50, LOWER_BOUND, UPPER_BOUND);
    //TH2D *frame1 = new TH2D("frame1",";M_{Q} [TeV]; p-value for Q#bar{Q} #rightarrow nV",10, NUM_min, NUM_max, 50, LOWER_BOUND, UPPER_BOUND);
	frame1->SetStats(0);
    frame1->GetYaxis()->SetTitleOffset(1.2);
    frame1->Draw();

    TGraph *graph_pvalue = new TGraph(NUM, MASS, p_value);
    graph_pvalue->SetMarkerStyle(20);
    graph_pvalue->SetLineWidth(2);
    graph_pvalue->SetLineColor(kBlue);
    graph_pvalue->Draw("pl,same");

    TGraph *graph_pvalue_pkg = new TGraph(NUM, MASS, p_value_pkg);
    graph_pvalue_pkg->SetMarkerStyle(20);
    graph_pvalue_pkg->SetLineWidth(2);
    graph_pvalue_pkg->SetLineColor(kRed);
    graph_pvalue_pkg->Draw("pl,same");

    legend->AddEntry(graph_pvalue,"Cowan","lp");
    legend->AddEntry(graph_pvalue_pkg,"package","lp");
    legend->SetTextSize(0.04);
    legend->Draw("same");

    text->SetTextFont(43);
    text->SetTextSize(18);
    text->SetTextColor(kBlack);
	text->DrawLatex(0.425,UPPER_BOUND*1.5,"MC Simulation from MadGraph5");
	text->DrawLatex(1.750,UPPER_BOUND*1.5,"L = 100 fb^{-1} (13TeV)");
    text->SetTextSize(24);
    for(int k=0; k<NUM_Z; k++){
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(Color[k]);
        line->DrawLine(NUM_min,Z_p_value[k],NUM_max,Z_p_value[k]);
        text->SetTextColor(Color[k]);
        text->DrawLatex(2.2, Z_p_value[k]/3, Form("%.0f #sigma",Z[k]));
    }

    c1->SaveAs("pvalue.png");

    //=== Significance ===//
	TCanvas *c2  = new TCanvas("c2","",810,10,800,600);
    //TH2D *frame2 = new TH2D("frame2",";M_{Q} [TeV]; significance for Q#bar{Q} #rightarrow nV",10, NUM_min, NUM_max, 50, -1, 25);
    TH2D *frame2 = new TH2D("frame2",";M_{Q} [TeV]; Significance",10, NUM_min, NUM_max, 50, -1, 25);
	frame2->SetStats(0);
    frame2->Draw();

    TGraph *graph_significance = new TGraph(NUM, MASS, significance);
    graph_significance->SetMarkerStyle(20);
    graph_significance->SetLineWidth(2);
    graph_significance->SetLineColor(kBlue);
    graph_significance->Draw("pl,same");

    TGraph *graph_significance_pkg = new TGraph(NUM, MASS, significance_pkg);
    graph_significance_pkg->SetMarkerStyle(20);
    graph_significance_pkg->SetLineWidth(2);
    graph_significance_pkg->SetLineColor(kRed);
    graph_significance_pkg->Draw("pl,same");

    legend->Clear();
    legend->AddEntry(graph_significance,"Cowan","lp");
    legend->AddEntry(graph_significance_pkg,"package","lp");
    legend->SetTextSize(0.04);
    legend->Draw("same");

    text->SetTextFont(43);
    text->SetTextSize(18);
    text->SetTextColor(kBlack);
	text->DrawLatex(0.425,25+0.3,"MC Simulation from MadGraph5");
	text->DrawLatex(1.750,25+0.3,"L = 100 fb^{-1} (13TeV)");
    text->SetTextSize(24);
    for(int k=0; k<NUM_Z; k++){
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(Color[k]);
        line->DrawLine(NUM_min,(double)k+1.,NUM_max,(double)k+1.);
        text->SetTextColor(Color[k]);
        text->DrawLatex(2.1, (double)k+1.-.8, Form("%.0f #sigma",Z[k]));
    }

    c2->SaveAs("significance.png");

    ////scale graph_significance to the pad coordinates
    //Float_t rightmax = ROOT::Math::normal_quantile_c(UPPER_BOUND);
    //for(int i=0; i<NUM; i++) graph_significance->SetPoint(i, MASS[i], ROOT::Math::normal_cdf_c(significance,1.0,0));
    //graph_significance->SetLineWidth(2);
    //graph_significance->SetLineColor(kRed);
    //graph_significance->Draw("pl,same");

    ////draw an axis on the right side
    //TGaxis *axis = new TGaxis(frame->GetXaxis()->GetXmax(), LOWER_BOUND, frame->GetXaxis()->GetXmax(), UPPER_BOUND, 0, rightmax, 510, "+L");
    //axis->SetLineColor(kRed);
    //axis->Draw();
}
