#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>

#include <algorithm>
#include <iterator>
#include <vector>
#include <iomanip>

#include "ana_draw.h"

using namespace std;
const double CUT[NUM] = {2,6,30,40,20,1200,50,1400,0};

//void ana_draw(int SIG, int XERR){
void ana_draw(){
SIG=2; XERR=1;

	if(SIG==1) importHist("fireball_bp_1TeV");
	if(SIG==2) importHist("fireball_bp_2TeV");
	importHist("pptt");
	importHist("ppvv");
	importHist("ppvtt");
	importHist("ppvvv");
	importHist("ppvvtt");
	importHist("pptttt");
	importHist("ppvvvv");
	importHist("ppvvvtt");

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

	//stack background!
	for(int k=0; k<NUM; k++){
		for(unsigned int i=0; i<vec.size(); i++){
			if(i+1==vec.size()) continue;//signal
			stack_ori[k]->Add(vec.at(vec.size()-i-1).hist_Ori[k]);
			stack_nm1[k]->Add(vec.at(vec.size()-i-1).hist_Nm1[k]);
		}
	}

	line->SetLineColor(kBlack);
	line->SetLineWidth(3);

	//calculate significance
	double yield_sig = 0, yield_bg=0, significance;
	for(std::vector<MyProcess>::iterator it=vec.begin(); it!=vec.end(); it++){
		cout<<setw(18)<<it->Name<<"\t"<<"L*X/N = "<<setw(8)<<it->Factor<<"\t"<<"N = "<<setw(8)<<it->Entries<<endl;
		if(it!=vec.begin()) yield_bg += it->Yield;
	}
	yield_sig = vec.at(0).Yield;
	significance = sqrt( 2*( (yield_sig + yield_bg)*TMath::Log(1+(yield_sig/yield_bg)) - yield_sig ) );
	cout<<"significance = "<<significance<<endl;

	//make plot
	TCanvas *can = new TCanvas("can","",800,600);
	for(int k=0; k<NUM; k++){
		can->SetLogy();
		stack_ori[k]->Draw("hist");
		vec.at(0).hist_Ori[k]->Draw("hist,same");//"hf,same"
		legend->Draw("same");
		can->SaveAs(Form("./output/stack_ori_%s.root",histName[k]));
		can->SaveAs(Form("./output/stack_ori_%s.png",histName[k]));

		stack_nm1[k]->Draw("hist");
		vec.at(0).hist_Nm1[k]->Draw("hist,same");
		legend->Draw("same");
		cout<<"Max="<<stack_nm1[k]->GetMaximum()<<endl;
		line->DrawLine(CUT[k],0,CUT[k],stack_nm1[k]->GetMaximum());
		can->SaveAs(Form("./output/stack_nm1_%s.root",histName[k]));
		can->SaveAs(Form("./output/stack_nm1_%s.png",histName[k]));
	}

}
