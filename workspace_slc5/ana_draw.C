#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include "ana_draw.h"

using namespace std;
const double CUT[NUM] = {2,6,30,40,20,1200,50,1400,0};

//void ana_draw(int SIG, int XERR){
void ana_draw(){
	SIG=2; XERR=1;

	if(SIG==1) importHist("fireball_bp_1TeV");
	if(SIG==2) importHist("fireball_bp_2TeV");
	importHist("ppvvvtt");
	importHist("ppvvvv");
	importHist("pptttt");
	importHist("ppvvtt");
	importHist("ppvvv");
	importHist("ppvtt");
	importHist("ppvv");
	importHist("pptt");

	int Nbin; double binWidth, binMin, binMax;
	for(int k=0; k<NUM; k++){
		Nbin   = vec.at(0).hist_Ori[k]->GetXaxis()->GetNbins();
		binMin = vec.at(0).hist_Ori[k]->GetXaxis()->GetXmin();
		binMax = vec.at(0).hist_Ori[k]->GetXaxis()->GetXmax();
		binWidth = (binMax-binMin)/(double)Nbin;
		if(k==2 || k==3 || k==4 || k==5 || k==6 || k==7){
		stack_ori[k] =new THStack(Form("Ori_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("Entries / %.0f %s",binWidth,units[k])));
		stack_nm1[k] =new THStack(Form("Cut_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("Entries / %.0f %s",binWidth,units[k])));
		} else{
		stack_ori[k] =new THStack(Form("Ori_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],"Entries"));
		stack_nm1[k] =new THStack(Form("Cut_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],"Entries"));
		}
	}

	for(int k=0; k<NUM; k++){
		for(std::vector<MyProcess>::iterator it=vec.begin(); it!=vec.end(); it++){
			if(it==vec.begin()) continue;
			stack_ori[k]->Add(it->hist_Ori[k]);
			stack_nm1[k]->Add(it->hist_Nm1[k]);
		}
	}

	line->SetLineColor(kBlack);
	line->SetLineWidth(3);

	
	for(std::vector<MyProcess>::iterator it=vec.begin(); it!=vec.end(); it++){
		cout<<it->Name<<"\tL*X/N = "<<it->Factor<<endl;
	}
	TCanvas *can = new TCanvas("can","",800,600);
	for(int k=0; k<NUM; k++){
		can->SetLogy();
		stack_ori[k]->Draw("hist");
		vec.at(0).hist_Ori[k]->Draw("hf,same");
		legend->Draw("same");
		can->SaveAs(Form("./output/stack_ori_%s.root",histName[k]));
		can->SaveAs(Form("./output/stack_ori_%s.png",histName[k]));

		stack_nm1[k]->Draw("hist");
		vec.at(0).hist_Nm1[k]->Draw("hf,same");
		legend->Draw("same");
		cout<<"Max="<<stack_nm1[k]->GetMaximum()<<endl;
		line->DrawLine(CUT[k],0,CUT[k],stack_nm1[k]->GetMaximum());
		can->SaveAs(Form("./output/stack_nm1_%s.root",histName[k]));
		can->SaveAs(Form("./output/stack_nm1_%s.png",histName[k]));
	}
	
}
