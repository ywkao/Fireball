#ifndef _LEAFDRAW_C_
#define _LEAFDRAW_C_
#include <string>
#include <iostream>
#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "../header/LeafDraw.h"
#include "../header/ExtractParticleID.h"
using namespace std;

//=== Delare Functions and Varibable ===//
string GetXtitle(const char*);
void FillHist(TH1D* , double );
void DrawHist(TH1D* );
void DrawHistTogether(TH1D*, TH1D* , TH1D* );
void DrawHistSetting(TH1D*&, int, int, int, int, int);
void Group_SetBranchAddress(TChain *, MyGenParticle &, const char*);
void Group_FillHist(MyGenParticle &);
void Group_DrawHist(MyGenParticle &);
void Group_DrawHistTogether(MyGenParticle &, MyGenParticle &, MyGenParticle &);
const char * DirOutput = "/home/xiaokao/Desktop/MG5_aMC_v2_2_3/Delphes-3.2.0/workspace/tmp";
TCanvas    * c1        = new TCanvas("c1","", 800, 600);
TPad       * pad       = new TPad("pad","",0,0,1,1);
//--------------------------------------------------//
//=== Main Function ===//
void LeafDraw(){
    gStyle->SetOptStat(0);

    TChain *chain = new TChain("mytree");
    chain->Add("/home/xiaokao/Desktop/MG5_aMC_v2_2_3/Delphes-3.2.0/workspace/result.root");

    MyGenParticle VB("VB"), W("W"), Z("Z"), Jet("Jet"), Lep("Lep");// histogram is set well
    Group_SetBranchAddress(chain, VB , "VB");
    Group_SetBranchAddress(chain, W  , "W" );
    Group_SetBranchAddress(chain, Z  , "Z" );
    Group_SetBranchAddress(chain, Jet, "Jet");
    Group_SetBranchAddress(chain, Lep, "Lep");
    
    for (int i=0; i<chain->GetEntries(); i++) {
        chain->GetEntry(i);
        Group_FillHist(VB);
        Group_FillHist(Z);
        Group_FillHist(W);
        Group_FillHist(Lep);
        Group_FillHist(Jet);
    }

    Group_DrawHist(VB);
    Group_DrawHist(Z);
    Group_DrawHist(W);
    Group_DrawHist(Lep);
    Group_DrawHist(Jet);

    Group_DrawHistTogether(VB, Z, W);

}
//--------------------------------------------------//
//===  Functions ===//
void DrawHist(TH1D* hist){
    //=== Pad Setting ===//
    c1->cd();
    pad->Draw();
    pad->cd();
    pad->SetTopMargin(0.1);
    pad->SetLeftMargin(0.12);
    //=== Title Setting ===//
    const char *Title = hist->GetTitle();
    hist->GetXaxis()->SetTitle( GetXtitle(Title).c_str() );
    if( (string)Title == "Multiplicity" || (string)Title == "Eta" || (string)Title == "Phi" ) hist->GetYaxis()->SetTitle( "Entries" );
    else hist->GetYaxis()->SetTitle( Form("Entries / %.0f [GeV]",hist->GetXaxis()->GetXmax()/(double)hist->GetNbinsX()) );
    hist->SetTitleOffset(1.6, "Y");
    //=== Drawing ===//
    hist->Draw();
    TLatex * latex = new TLatex(0,0,"");
    latex -> SetNDC();
	latex -> SetTextFont(43);
    latex -> SetTextSize(22);
    latex -> DrawText(0.68 , 0.72 , Form("Entries: %.0f" , hist->GetEntries()));
    latex -> DrawText(0.68 , 0.66 , Form("Mean:    %.2f" , hist->GetMean()));
    latex -> DrawText(0.68 , 0.60 , Form("RMS:     %.2f" , hist->GetStdDev()));
    latex -> SetTextSize(32);
    latex -> DrawLatex(0.68 , 0.80 , Form("#font[22]{#color[4]{%s}}",ExtractParticleID(hist->GetName()).c_str()));
    latex -> SetTextSize(16);
    latex -> DrawText (0.13 , 0.91,"Delphes v3.3 GenParticleInfo");
    latex -> DrawLatex(0.72 , 0.91,"L = 100 fb^{-1} (13TeV)");

    //hist->SetTitle("");
    c1->SaveAs(Form("%s/%s.png", DirOutput, hist->GetName()));
}
void DrawHistTogether(TH1D* hist, TH1D* hist_1, TH1D* hist_2){
    //=== Pad Setting ===//
    c1->cd();
    pad->Draw();
    pad->cd();
    pad->SetTopMargin(0.1);
    pad->SetLeftMargin(0.12);
    //=== Title Setting ===//
    printf("title = %s\n", hist->GetTitle());
    const char *Title = hist->GetTitle();
    hist->GetXaxis()->SetTitle( GetXtitle(Title).c_str() );
    if( (string)Title == "Multiplicity" || (string)Title == "Eta" || (string)Title == "Phi" ) hist->GetYaxis()->SetTitle( "Entries" );
    else hist->GetYaxis()->SetTitle( Form("Entries / %.0f [GeV]",hist->GetXaxis()->GetXmax()/(double)hist->GetNbinsX()) );
    hist->SetTitleOffset(1.6, "Y");
    //=== Drawing ===//
    DrawHistSetting(hist   , kRed     , 1 , 2 , kRed     , 3004);
    DrawHistSetting(hist_1 , kGreen+2 , 1 , 2 , kGreen+2 , 3004);
    DrawHistSetting(hist_2 , kBlue+1  , 1 , 2 , kBlue+1  , 3004);
    if( (string)Title == "Multiplicity")
    hist   -> SetMaximum(4500);
    hist   -> Draw();
    hist_1 -> Draw("same");
    hist_2 -> Draw("same");
    TLegend * legend = new TLegend(0.65,0.60,0.80,0.85);
    legend->SetTextFont  ( 43 );
    legend->SetTextSize  ( 20 );
    legend->SetFillStyle ( 0  );
    legend->SetFillColor ( 0  );
    legend->SetBorderSize( 0  );
    legend->AddEntry(hist  , ExtractParticleID(hist  ->GetName()).c_str(), "l");
    legend->AddEntry(hist_1, ExtractParticleID(hist_1->GetName()).c_str(), "l");
    legend->AddEntry(hist_2, ExtractParticleID(hist_2->GetName()).c_str(), "l");
    legend->Draw("same");
    TLatex  * latex  = new TLatex(0,0,"");
    latex -> SetNDC();
	latex -> SetTextFont(43);
    latex -> SetTextSize(22);
    latex -> SetTextSize(32);
    latex -> SetTextSize(16);
    latex -> DrawText (0.13 , 0.91,"Delphes v3.3 GenParticleInfo");
    latex -> DrawLatex(0.72 , 0.91,"L = 100 fb^{-1} (13TeV)");

    hist->SetTitle("");
    c1->SaveAs(Form("%s/Combined_%s.png", DirOutput, hist->GetName()));
}
string GetXtitle(const char* Title){
    string xtitle;
    if((string)Title == "Multiplicity") xtitle = "Multiplicity";
    if((string)Title == "Energy"      ) xtitle = "Energy [GeV]";
    if((string)Title == "Momentum"    ) xtitle = "Momentum [GeV]";
    if((string)Title == "Momentum2"   ) xtitle = "Momentum [GeV]";
    if((string)Title == "PT"          ) xtitle = "PT [GeV]";
    if((string)Title == "Eta"         ) xtitle = "Eta";
    if((string)Title == "Phi"         ) xtitle = "Phi [rad]";
    return xtitle;
}
void FillHist(TH1D* hist, double value){
    if(value==-999){}
    else hist->Fill(value);
}
void Group_SetBranchAddress(TChain *chain, MyGenParticle &Particle, const char* Type){
    chain->SetBranchAddress(Form("%s_Multiplicity", Type), &Particle.Multiplicity );
    chain->SetBranchAddress(Form("%s_Energy"      , Type), &Particle.Energy       );
    chain->SetBranchAddress(Form("%s_Momentum"    , Type), &Particle.Momentum     );
    chain->SetBranchAddress(Form("%s_Momentum2"   , Type), &Particle.Momentum2    );
    chain->SetBranchAddress(Form("%s_PT"          , Type), &Particle.PT           );
    chain->SetBranchAddress(Form("%s_Eta"         , Type), &Particle.Eta          );
    chain->SetBranchAddress(Form("%s_Phi"         , Type), &Particle.Phi          );
}
void Group_FillHist(MyGenParticle &Particle){
    FillHist( Particle.hist_Multiplicity, Particle.Multiplicity );
    FillHist( Particle.hist_Energy      , Particle.Energy       );
    FillHist( Particle.hist_Momentum    , Particle.Momentum     );
    FillHist( Particle.hist_Momentum2   , Particle.Momentum2    );
    FillHist( Particle.hist_PT          , Particle.PT           );
    FillHist( Particle.hist_Eta         , Particle.Eta          );
    FillHist( Particle.hist_Phi         , Particle.Phi          );
}
void Group_DrawHist(MyGenParticle &Particle){
    DrawHist( Particle.hist_Multiplicity );
    DrawHist( Particle.hist_Energy       );
    DrawHist( Particle.hist_Momentum     );
    DrawHist( Particle.hist_Momentum2    );
    DrawHist( Particle.hist_PT           );
    DrawHist( Particle.hist_Eta          );
    DrawHist( Particle.hist_Phi          );
}
void Group_DrawHistTogether(MyGenParticle &Particle, MyGenParticle &Particle1, MyGenParticle &Particle2){
    DrawHistTogether( Particle.hist_Multiplicity, Particle1.hist_Multiplicity, Particle2.hist_Multiplicity);
    DrawHistTogether( Particle.hist_Energy      , Particle1.hist_Energy      , Particle2.hist_Energy      );
    DrawHistTogether( Particle.hist_Momentum    , Particle1.hist_Momentum    , Particle2.hist_Momentum    );
    DrawHistTogether( Particle.hist_Momentum2   , Particle1.hist_Momentum2   , Particle2.hist_Momentum2   );
    DrawHistTogether( Particle.hist_PT          , Particle1.hist_PT          , Particle2.hist_PT          );
    DrawHistTogether( Particle.hist_Eta         , Particle1.hist_Eta         , Particle2.hist_Eta         );
    DrawHistTogether( Particle.hist_Phi         , Particle1.hist_Phi         , Particle2.hist_Phi         );
}
void DrawHistSetting(TH1D *&hist, int lcolor, int lstyle, int lwidth, int fcolor, int style){
	hist->SetLineColor(lcolor);
	hist->SetLineStyle(lstyle);
	hist->SetLineWidth(lwidth);
	hist->SetFillColor(fcolor);
	hist->SetFillStyle(style);
}
#endif
