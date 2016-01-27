#include <iostream>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include "xsec.h"
using namespace std;

const int NUM = 9;
const double L = 100, pb2fb = 1000;
const char *histName[NUM]={"Num_lep","Num_jet","PT_lep","PT_jet","tot_Lep_PT","HT","MET","ST","Num_boson"};
const char *units[NUM] = {"# of lep", "# of jet", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "# of boson"};
//const char* PREFIX = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace/skimmed/result";
const char* PREFIX = "~/Desktop/workspace/skimmed/result";
const char* SUFFIX = "root";
const char* TARGET;
int SIG, XERR;

class MyProcess{
public:
	const char* Name;
	double Factor, Entries, X, Err_X, Yield;
	TH1D *hist_Ori[NUM];
	TH1D *hist_Cut[NUM];
	TH1D *hist_Nm1[NUM];
};

struct CrossSection {
	//std::pair<double,double> pptt = std::make_pair(815.96,45.51); WHY DOES IT NOT WORK??
	double GetCrossSection(const char* Name){
		if((string)Name=="fireball_bp_1TeV") 	return fireball_bp_1TeV;
		if((string)Name=="fireball_bp_2TeV") 	return fireball_bp_2TeV;
		if((string)Name=="pptt")    return pptt;
		if((string)Name=="ppvv")    return ppvv;
		if((string)Name=="ppvtt")   return ppvtt;
		if((string)Name=="ppvvv")   return ppvvv;
		if((string)Name=="ppvvtt")  return ppvvtt;
		if((string)Name=="pptttt")  return pptttt;
		if((string)Name=="ppvvvv")  return ppvvvv;
		if((string)Name=="ppvvvtt") return ppvvvtt;
	}
	double GetCrossSectionErr(const char* Name){
		if((string)Name=="fireball_bp_1TeV") 	return Err_fireball_bp_1TeV;
		if((string)Name=="fireball_bp_2TeV") 	return Err_fireball_bp_2TeV;
		if((string)Name=="pptt") 	return Err_pptt;
		if((string)Name=="ppvv")    return Err_ppvv;
		if((string)Name=="ppvtt")   return Err_ppvtt;
		if((string)Name=="ppvvv")   return Err_ppvvv;
		if((string)Name=="ppvvtt")  return Err_ppvvtt;
		if((string)Name=="pptttt")  return Err_pptttt;
		if((string)Name=="ppvvvv")  return Err_ppvvvv;
		if((string)Name=="ppvvvtt") return Err_ppvvvtt;
	}
}; CrossSection Xsec;

int GetColor(const char* Name){
	if((string)Name=="fireball_bp_1TeV") return kRed;//return kRed;
	if((string)Name=="fireball_bp_2TeV") return kRed;//return kRed;
	if((string)Name=="pptt") 			 return kPink-4;//return kOrange+3;
	if((string)Name=="ppvv")    		 return kOrange-2;//kOrange-4;//return kOrange+1;
	if((string)Name=="ppvtt")   		 return kYellow-7;//kYellow-7;//return kBlue;
	if((string)Name=="ppvvv")   		 return kGreen-3;//return kCyan;
	if((string)Name=="ppvvtt")  		 return kGreen+2;//return kGreen+2;
	if((string)Name=="pptttt")  		 return kViolet+6;//kViolet+8;//return kGreen;
	if((string)Name=="ppvvvv")  		 return kMagenta-6;//return kMagenta+2;
	if((string)Name=="ppvvvtt") 		 return kViolet-6;//return kMagenta;
}

void ScaleAndDrawSetting(TH1D *&hist, double &Yield, double factor, double X, double Err_X, int lcolor, int lstyle, int lwidth, int fcolor, int style){
	double content, error, yield=0;
	for(int i=0; i<hist->GetNbinsX(); i++){
		content = hist->GetBinContent(i+1);
		error   = hist->GetBinError(i+1);
		yield  += content*factor;

		hist->SetBinContent(i+1,content*factor);
		if(XERR==0) hist->SetBinError(i+1,error*factor);
		else		hist->SetBinError(i+1, (pow(Err_X/X,2)+ pow(error/content,2))*content*factor );
	}
	Yield = yield;
	hist->SetLineColor(lcolor);
	hist->SetLineStyle(lstyle);
	hist->SetLineWidth(lwidth);
	hist->SetFillColor(fcolor);
	hist->SetFillStyle(style);
}

std::vector<MyProcess> vec;
THStack *stack_ori[NUM], *stack_nm1[NUM];//stack background
TLegend *legend = new TLegend(0.74,0.5,0.9,0.9);
TLine 	*line	= new TLine(0,0,0,0);

void importHist(const char *NAME){
	TARGET = Form("%s_%s.%s", PREFIX, NAME, SUFFIX);
	TFile *fin = new TFile(TARGET,"READ");

	MyProcess process;
	process.Name = NAME;
	process.Entries = ( (TH1D*)fin->Get(Form("Ori_%s",histName[0])) )->GetEntries();//get the # of events of the process
	process.X = Xsec.GetCrossSection(NAME);
	process.Err_X = Xsec.GetCrossSectionErr(NAME);
	process.Factor = L*pb2fb*process.X/process.Entries;

	int lcolor, lstyle, lwidth, style;
	if ((string)NAME=="fireball_bp_1TeV" || (string)NAME=="fireball_bp_2TeV") {lcolor=2; lstyle=1; lwidth=2; style=3004;}
	else {lcolor=0; lstyle=0; lwidth=0; style=1001;}

	for(int k=0; k<NUM; k++){
		process.hist_Ori[k] = (TH1D*) fin->Get(Form("Ori_%s",histName[k]));
		process.hist_Cut[k] = (TH1D*) fin->Get(Form("Cut_%s",histName[k]));
		process.hist_Nm1[k] = (TH1D*) fin->Get(Form("N-1_%s",histName[k]));
		ScaleAndDrawSetting(process.hist_Ori[k], process.Yield, process.Factor, process.X, process.Err_X, GetColor(NAME), lstyle, lwidth, GetColor(NAME), style);
		ScaleAndDrawSetting(process.hist_Cut[k], process.Yield, process.Factor, process.X, process.Err_X, GetColor(NAME), lstyle, lwidth, GetColor(NAME), style);
		ScaleAndDrawSetting(process.hist_Nm1[k], process.Yield, process.Factor, process.X, process.Err_X, GetColor(NAME), lstyle, lwidth, GetColor(NAME), style);
	}

	legend->AddEntry(process.hist_Ori[0],NAME,"f");
	vec.push_back(process);
}

//###for finding segmentation error
//std::cout<<"[INFO] SignalHandler"<<std::endl;
//signal(SIGSEGV,sighandler);
void sighandler(int sig) {while(1);}
