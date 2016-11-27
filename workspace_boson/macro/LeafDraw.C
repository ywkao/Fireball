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
#include "../header/ExtractMassFromFilename.h"
#include "../header/FileGenerator.h"
#include "../header/Normalization.h"
#include "../header/BosonMultiplicityTable.h"
using namespace std;

//=== Delare Functions and Global Varibable ===//
string GetXtitle(const char*);
void FillHist(TH1D* , double );
void DrawHist(TH1D* );
void DrawHistTogether(TH1D*, TH1D* , TH1D* );
void DrawHistSetting(TH1D*, int, int, int, int, int);
void Group_SetBranchAddress(TChain *, MyGenParticle &, const char*);
void Group_FillHist(MyGenParticle &);
void Group_NormalizeHist(MyGenParticle &);
void Group_DrawHist(MyGenParticle &);
void Group_DrawHistTogether(MyGenParticle &, MyGenParticle &, MyGenParticle &);
TCanvas    * c1        = new TCanvas("c1","", 800, 600);
TPad       * pad       = new TPad("pad","",0,0,1,1);
double factor, relative_error_xsec, total_generated_entries;

string DirOutput;
const char * DirInput  = "/home/xiaokao/Desktop/MG5_aMC_v2_2_3/Delphes-3.2.0/workspace/skimmed";

bool SavePlots = false;
bool TestOnly  = false;

//=================================================================//
//========================= Main Function =========================//
//=================================================================//
void LeafDraw(const char* filename){
    gStyle->SetOptStat(0);

    TChain *chain = new TChain("mytree");
    if((string)filename=="/home/xiaokao/Desktop/MG5_aMC_v2_2_3/Delphes-3.2.0/workspace/result.root"){
        DirOutput = "/home/xiaokao/Desktop/MG5_aMC_v2_2_3/Delphes-3.2.0/workspace/tmp";
        chain->Add(filename);
        TestOnly = true;
    } else{ //=== For simulated fireball process w/ vaious masses ===//
        DirOutput = Form("/home/xiaokao/Desktop/MG5_aMC_v2_2_3/Delphes-3.2.0/workspace/output/%s",filename);
        chain->Add(Form("%s/%s.root"  ,DirInput,filename));
        chain->Add(Form("%s/%s_x.root",DirInput,filename));
        if((string)filename=="result_run_fireball_1TeV") SavePlots = true;
        if((string)filename=="result_run_fireball_2TeV") SavePlots = true;
    }


    MyGenParticle VB("VB"), W("W"), Z("Z"), GenJet("GenJet"), GenLep("GenLep");// histogram is set well
    Group_SetBranchAddress(chain, VB , "VB");
    Group_SetBranchAddress(chain, W  , "W" );
    Group_SetBranchAddress(chain, Z  , "Z" );
    Group_SetBranchAddress(chain, GenJet, "GenJet");
    Group_SetBranchAddress(chain, GenLep, "GenLep");
    
    for (int i=0; i<chain->GetEntries(); i++) {
        chain->GetEntry(i);
        Group_FillHist(VB);
        Group_FillHist(Z);
        Group_FillHist(W);
        Group_FillHist(GenLep);
        Group_FillHist(GenJet);
    }

    //=== Make Boson Multiplicity Table ===//
    if(!TestOnly)
    BosonMultiplicity_SingleMass(ExtractMassFromFilename(filename), VB.hist_Multiplicity, Z.hist_Multiplicity, W.hist_Multiplicity);

    //=== Normalization ===//
    total_generated_entries = GenLep.hist_Multiplicity->GetEntries();
    factor = L*pb2fb*Xsec.GetCrossSection("result", 1)/total_generated_entries;
    relative_error_xsec = Xsec.GetCrossSectionErr( "result", 1 ) / Xsec.GetCrossSection( "result", 1 );// 1: apply K-factor
    printf("factor = %f \n", factor);

    //The function MC_Normalization would force hist increase 50 entries!!?
    Group_NormalizeHist(VB);
    Group_NormalizeHist(Z);
    Group_NormalizeHist(W);
    Group_NormalizeHist(GenLep);
    Group_NormalizeHist(GenJet);

    printf("Before Nor. = %f, After Nor. = %f\n", total_generated_entries, GenLep.hist_Multiplicity->GetEntries());

    //=== Make Plots ===//
    Group_DrawHist(VB);
    Group_DrawHist(Z);
    Group_DrawHist(W);
    Group_DrawHist(GenLep);
    Group_DrawHist(GenJet);

    Group_DrawHistTogether(VB, Z, W);

}


//=================================================================//
//========================= Aux. Function =========================//
//=================================================================//
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
    hist->Draw("hist");
    TLatex * latex = new TLatex(0,0,"");
    latex -> SetNDC();
	latex -> SetTextFont(43);
    latex -> SetTextSize(22);
    //printf("GetMeanValue = %f\n", GetMeanValue(hist));
    latex -> DrawText(0.68 , 0.72 , Form("Entries: %.0f" , hist->GetEntries()));
    latex -> DrawText(0.68 , 0.66 , Form("Mean:    %.2f" , hist->GetMean()));
    latex -> DrawText(0.68 , 0.60 , Form("RMS:     %.2f" , hist->GetStdDev()));
    latex -> SetTextSize(32);
    latex -> DrawLatex(0.68 , 0.80 , Form("#font[22]{#color[4]{%s}}",ExtractParticleID(hist->GetName()).c_str()));
    latex -> SetTextSize(16);
    latex -> DrawText (0.13 , 0.91,"Delphes v3.3 GenParticleInfo");
    latex -> DrawLatex(0.72 , 0.91,"L = 100 fb^{-1} (13TeV)");

    //hist->SetTitle("");
    if(SavePlots) c1->SaveAs(Form("%s/%s.png", DirOutput.c_str(), hist->GetName()));
}
void DrawHistTogether(TH1D* hist, TH1D* hist_1, TH1D* hist_2){
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
    DrawHistSetting(hist   , kRed     , 1 , 2 , kRed     , 3004);
    DrawHistSetting(hist_1 , kGreen+2 , 1 , 2 , kGreen+2 , 3004);
    DrawHistSetting(hist_2 , kBlue+1  , 1 , 2 , kBlue+1  , 3004);
    if( (string)Title == "Multiplicity"){
        hist_2 -> Draw("hist"); hist_2->SetTitle("");//W
        hist   -> Draw("hist,same");
        hist_1 -> Draw("hist,same");
        hist_2 -> Draw("hist,same");
    } else{
        hist   -> Draw("hist");
        hist_1 -> Draw("hist,same");
        hist_2 -> Draw("hist,same");
    }
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
    if(SavePlots) c1->SaveAs(Form("%s/Combined_%s.png", DirOutput.c_str(), hist->GetName()));
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
void Group_NormalizeHist(MyGenParticle &Particle){//Normalize each historgram
    MC_Normalization( Particle.hist_Multiplicity, factor, relative_error_xsec );
    MC_Normalization( Particle.hist_Energy      , factor, relative_error_xsec );
    MC_Normalization( Particle.hist_Momentum    , factor, relative_error_xsec );
    MC_Normalization( Particle.hist_Momentum2   , factor, relative_error_xsec );
    MC_Normalization( Particle.hist_PT          , factor, relative_error_xsec );
    MC_Normalization( Particle.hist_Eta         , factor, relative_error_xsec );
    MC_Normalization( Particle.hist_Phi         , factor, relative_error_xsec );
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
void DrawHistSetting(TH1D *hist, int lcolor, int lstyle, int lwidth, int fcolor, int style){
	hist->SetLineColor(lcolor);
	hist->SetLineStyle(lstyle);
	hist->SetLineWidth(lwidth);
	hist->SetFillColor(fcolor);
	hist->SetFillStyle(style);
}
#endif
