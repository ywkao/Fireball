#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TText.h>
#include "ana_draw.h"
using namespace std;

const int	 NBins=120;
const int	 RebinFactor=4;
const double XMax=3000;

const int 	 NUM=5;
const double MAX=1e+06, MIN=1, ScopeMax=7000;
const int 	 values[NUM]={500,1000,1500,2000,2500};
const double position[NUM]={4e+02,2e+02,1e+02,50,25};
Int_t 		 color[NUM]={kBlack,kOrange+2,kBlue,kGreen,kRed};

const double L=100., pb2fb=1000., K=1.69;
const double xsec[NUM]={9.292,0.391,0.033,0.00395,0.0005875};
const double xerr[NUM]={0.096,0.0031,0.00022,2.2e-05,2.8e-06};
double yield, err, factor;
int nentries;


//void Draw_ihtmin_GenHT(const char* TagName, const int STATUS){
void Draw_ihtmin_GenHT(const char* TagName, const int STATUS){
	gStyle->SetOptStat(0);
	TCanvas *c1 = new TCanvas("can","",800,600); c1->SetLogy();
	TLegend *legend = new TLegend(0.6,0.6,0.9,0.9);
	//TH1D 	*HIST = new TH1D("hist",Form("status=3, %s-level Particle HT distribution w/ different \"ihtmin\" values;HT [GeV];Entries",TagName),56,0,7000);
	TH1D 	*HIST = new TH1D("hist",Form(";HT [GeV];Entries / %.0f GeV",XMax/(double)(NBins/RebinFactor)),240,0,7000);
	TH1D 	*hist[NUM];
	TNtuple *ntuple;

	TFile 	*file;
	TLine 	*line = new TLine(0,0,0,0); 
	TLatex  *text = new TLatex(0,0,"");
	text->SetTextFont(43);
	text->SetTextSize(18);

	for(int i=0; i<NUM; i++ ){
		if(STATUS==-1) file = new TFile(Form("skimmed_%s_status_m1/result_pptt_%s_ihtmin_%d.root",TagName,TagName,values[i]),"READ");
		else file = new TFile(Form("skimmed_%s_status_%d/result_pptt_%s_ihtmin_%d.root",TagName,STATUS,TagName,values[i]),"READ");
		//file = new TFile(Form("skimmed_iht/result_pptt_%s_ihtmin_%d.root",TagName,values[i]),"READ");
		//file = new TFile(Form("skimmed_pythia_status_1/result_pptt_%s_ihtmin_%d.root",TagName,values[i]),"READ");
		//file = new TFile(Form("skimmed_iht/result_pptt_%s_ht_min_%d.root",TagName,values[i]),"READ");
		//if(TagName=="gen") hist[i]  = (TH1D*) file->Get(Form("hist_%s_HT",TagName));
		//else 			   hist[i]  = (TH1D*) file->Get(Form("HT",TagName));
		hist[i]  = (TH1D*) file->Get(Form("HT",TagName));
		nentries = hist[i]->GetEntries();
		factor   = L*pb2fb*xsec[i]/nentries;
		//ntuple = (TNtuple*) file->Get("ntuple_gen_HT");
		//ntuple->SetLineColor(color[i]);
		//ntuple->SetMarkerStyle(20);
		//nentries = ntuple->GetEntries();

		if(i==0){
			HIST->SetMaximum(MAX);
			HIST->SetMinimum(MIN);
			HIST->GetXaxis()->SetRangeUser(0.,ScopeMax);
			HIST->Draw();

			line->SetLineColor(kYellow-9);
			line->SetLineWidth(3);
			line->DrawLine(500,0,500,MAX);
			line->DrawLine(1500,0,1500,MAX);
			line->DrawLine(2500,0,2500,MAX);
			text->SetTextSize(16);
			text->DrawLatex(50,1.3e+06,"MC Simulation from MadGraph5");
			text->DrawLatex(5500,1.3e+06,"L = 100 fb^{-1} (13TeV)");
			text->SetTextSize(24);
			text->DrawLatex(2800,1e+05,Form("#font[22]{Status = %d}",STATUS));
		}

		ScaleAndDrawSetting(hist[i], yield, err, factor, xsec[i], xerr[i], color[i], 1, 1, 0, 0);
		printf("%s: yield = %f +/- %f, xsec = %f\n",hist[i]->GetName(), yield, err, xsec);

		legend->AddEntry(hist[i],Form("ihtmin = %d GeV",values[i]),"l");
		text->SetTextSize(16);
		if(i<3) text->DrawLatex(3500,position[i],Form("#color[%d]{Entries = %d; Yield = %d #pm %d}"	 ,color[i], nentries, (int)yield, (int)err));
		else	text->DrawLatex(3500,position[i],Form("#color[%d]{Entries = %d; Yield = %d #pm %.2f}",color[i], nentries, (int)yield, err));
		//hist[i]->SetFillColor(color[i]);
		//hist[i]->SetFillStyle(3004);
		hist[i]->Rebin(RebinFactor);
		hist[i]->Draw("he,same");
		//if(i<2) hist[i]->Draw("he,same");

		//legend->AddEntry(ntuple,Form("ihtmin = %d GeV",values[i]),"l");
		//text->DrawLatex(4000,position[i],Form("#color[%d]{Entries = %d,%03d}",color[i],(nentries/1000)%1000,nentries%1000));
		//ntuple->Draw("quantity","quantity>0","same");
	}
	legend->Draw("same");
	

	if(STATUS==-1){
		c1->SaveAs(Form("Comparison_ihtmin_%s_HT_status_m1.pdf",TagName,STATUS));
		c1->SaveAs(Form("Comparison_ihtmin_%s_HT_status_m1.png",TagName,STATUS));
		c1->SaveAs(Form("Comparison_ihtmin_%s_HT_status_m1.root",TagName,STATUS));
	} else{
		c1->SaveAs(Form("Comparison_ihtmin_%s_HT_status%d.pdf",TagName,STATUS));
		c1->SaveAs(Form("Comparison_ihtmin_%s_HT_status%d.png",TagName,STATUS));
		c1->SaveAs(Form("Comparison_ihtmin_%s_HT_status%d.root",TagName,STATUS));
	}

	//c1->Draw();
	//gROOT->ProcessLine();
}

