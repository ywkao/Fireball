#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <signal.h>
#include "xsec.h"
#define NUM      9
#define L        100
#define pb2fb    1000
#define KFACTOR1 1.56 //ppvv-like process
#define KFACTOR2 1.69 //pptt-like process
#define PREFIX   "~/Desktop/workspace/skimmed_new/result"
#define SUFFIX   "root"

//============================================================//
//==================== Function Prototype ====================//
//============================================================//
class MyProcess{
public:
	const char* Name;
	double Factor, Entries, X, Err_X, X_ori, Err_X_ori, KFactor;
	double Yield[3][NUM], Err_Yield[3][NUM], Eff[3][NUM], Err_Eff[3][NUM];
	TH1D *hist_Ori[NUM];
	TH1D *hist_Cut[NUM];
	TH1D *hist_Nm1[NUM];
};

struct CrossSection {
	double GetCrossSection(const char*, int);
	double GetCrossSectionErr(const char*, int);
}; CrossSection Xsec;

int GetColor(const char*);
const char* GetLegendName(const char*);
void MC_Normalization(const char*, TH1D*&, double&, double&, double&, double&, double, double, double, double);
void HistDrawSetting(TH1D*&, int, int, int, int, int);
void SetLegend(TLegend *&, int, int, int, int, int, int);
void importHist(const char*);

//============================================================//
//===================== Global Variables =====================//
//============================================================//
int SIG, XERR=1;
const char *TARGET;
const char *histName[NUM] = {"Num_lep","Num_jet","PT_lep","PT_jet","LPT","HT","MET","ST","Num_boson"};
const char *units[NUM]    = {"# of lep", "# of jet", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "# of boson"};
//=== To Aggregate All Processes ===// 
std::vector<MyProcess> vec;
//=== To Stack bkg ===// 
THStack *stack_ori[NUM];
THStack *stack_nm1[NUM];
TH1D *hist_ori_ppttv[NUM];
TH1D *hist_ori_ppvvv[NUM];
TH1D *hist_nm1_ppttv[NUM];
TH1D *hist_nm1_ppvvv[NUM];
TH1D *hist_err_ori[NUM];
TH1D *hist_err_nm1[NUM];
//=== To Beautify Plots ===// 
TGaxis  *axis   = new TGaxis(0,0,0,0,0,0,510,"");//xmin,ymin,xmax,ymax,wmin,wmax
TLatex  *text   = new TLatex(0,0,"");
TLine 	*line   = new TLine(0,0,0,0);
TLegend *legend_MC_sig = new TLegend(0.62,0.72,0.87,0.88); //leg->SetBorderSize(0);
TLegend *legend_MC_bkg = new TLegend(0.62,0.54,0.9,0.70); //legend->SetBorderSize(0); legend->SetNColumns(2);
const double MAX_ori=1e+07;
const double MAX_nm1=1e+04;
const double MIN = 0.1;

//============================================================//
//=================== Import MC Histograms ===================//
//============================================================//
void importHist(const char *NAME){
	TARGET = Form("%s_%s.%s", PREFIX, NAME, SUFFIX);
	TFile *fin = new TFile(TARGET,"READ");

	MyProcess process;
	process.Name      = NAME;
	process.Entries   = ( (TH1D*)fin->Get(Form("Ori_%s",histName[0])) )->GetEntries();//get the # of events of the process
	process.X         = Xsec.GetCrossSection(NAME,1);
	process.Err_X     = Xsec.GetCrossSectionErr(NAME,1);
	process.X_ori     = Xsec.GetCrossSection(NAME,0);   // w/o applying k-factor
	process.Err_X_ori = Xsec.GetCrossSectionErr(NAME,0);// w/o applying k-factor
    process.KFactor   = process.X/process.X_ori;
	process.Factor    = L*pb2fb*process.X/process.Entries;
	for(int i=0; i<3; i++) for(int j=0; j<NUM; j++) process.Yield[i][j]     = 0;
	for(int i=0; i<3; i++) for(int j=0; j<NUM; j++) process.Err_Yield[i][j] = 0;

	int lcolor, lstyle, lwidth, style;
	if      ((string)NAME=="fireball_1TeV") {lcolor=GetColor(NAME); lstyle=1; lwidth=2; style=3004;}
    else if ((string)NAME=="fireball_2TeV") {lcolor=GetColor(NAME); lstyle=1; lwidth=2; style=3004;}
	else                                    {lcolor=GetColor(NAME); lstyle=0; lwidth=0; style=1001;}

	for(int k=0; k<NUM; k++){
		process.hist_Ori[k] = (TH1D*) fin->Get(Form("Ori_%s",histName[k]));
		process.hist_Cut[k] = (TH1D*) fin->Get(Form("Cut_%s",histName[k]));
		process.hist_Nm1[k] = (TH1D*) fin->Get(Form("N-1_%s",histName[k]));
		MC_Normalization("Ori",process.hist_Ori[k], process.Yield[0][k], process.Err_Yield[0][k], process.Eff[0][k], process.Err_Eff[0][k],
                            process.Entries, process.Factor, process.X, process.Err_X);
		MC_Normalization("Nm1",process.hist_Nm1[k], process.Yield[1][k], process.Err_Yield[1][k], process.Eff[1][k], process.Err_Eff[1][k],
                            process.Entries, process.Factor, process.X, process.Err_X);
		MC_Normalization("Cut",process.hist_Cut[k], process.Yield[2][k], process.Err_Yield[2][k], process.Eff[2][k], process.Err_Eff[2][k],
                            process.Entries, process.Factor, process.X, process.Err_X);

        HistDrawSetting(process.hist_Ori[k], lcolor, lstyle, lwidth, GetColor(NAME), style);
        HistDrawSetting(process.hist_Nm1[k], lcolor, lstyle, lwidth, GetColor(NAME), style);
        HistDrawSetting(process.hist_Cut[k], lcolor, lstyle, lwidth, GetColor(NAME), style);
	}

	if(NAME=="fireball_1TeV" || NAME=="fireball_2TeV")                 legend_MC_sig->AddEntry(process.hist_Ori[0],GetLegendName(NAME),"f");
	if(NAME=="pptt" || NAME=="ppvv" || NAME=="ppvtt" || NAME=="ppvvv") legend_MC_bkg->AddEntry(process.hist_Ori[0],GetLegendName(NAME),"f");
	vec.push_back(process);
}

//============================================================//
//=================== Normalize MC by Lumi ===================//
//============================================================//
void MC_Normalization(const char* TAG, TH1D *&hist, double &Yield, double &Err_Yield, double &Eff, double &Err_Eff, 
                         double entries, double factor, double X, double Err_X){
	double content, error, yield=0, err_yield=0;
	for(int i=0; i<hist->GetNbinsX(); i++){
		content    = hist->GetBinContent(i+1);
		error      = hist->GetBinError(i+1);
		yield     += content*factor;
        err_yield += pow(error*factor,2);

		hist->SetBinContent(i+1,content*factor);
		if(XERR==0) hist->SetBinError(i+1,error*factor);
		else		hist->SetBinError(i+1, sqrt(pow(Err_X/X,2)+ pow(error/content,2))*content*factor );
	}
	Yield     = yield;
    err_yield = sqrt(err_yield);
    Err_Yield = sqrt(pow(Err_X/X,2)+ pow(err_yield/yield,2))*yield;
    Eff       = yield/(factor*entries);
    Err_Eff   = sqrt(Eff*(1-Eff)/entries);
    //Err_Eff   = err_yield/(factor*entries);
}

//============================================================//
//================= Miscellaneous Functions ==================//
//============================================================//
//void SetLegend(TLegend *&legend, int NColumns=1, int TextFont=43, int TextSize=16, int FillColor=0, int FillStyle=0, int BorderSize=0){ //INVALID
void SetLegend(TLegend *&legend, int NColumns, int TextFont, int TextSize, int FillColor, int FillStyle, int BorderSize){
    legend->SetNColumns  ( NColumns   );
    legend->SetTextFont  ( TextFont   );
    legend->SetTextSize  ( TextSize   );
    legend->SetFillStyle ( FillStyle  );
    legend->SetFillColor ( FillColor  );
    legend->SetBorderSize( BorderSize );
}
void HistDrawSetting(TH1D *&hist, int lcolor, int lstyle, int lwidth, int fcolor, int style){
	hist->SetLineColor(lcolor);
	hist->SetLineStyle(lstyle);
	hist->SetLineWidth(lwidth);
	hist->SetFillColor(fcolor);
	hist->SetFillStyle(style);
}
double CrossSection::GetCrossSection(const char* Name, int K){
    if((string)Name=="fireball_bp_1TeV") return fireball_bp_1TeV;
    if((string)Name=="fireball_bp_2TeV") return fireball_bp_2TeV;
    if((string)Name=="fireball_1TeV"  )  return (K==0 ? fireball_1TeV   : KFACTOR2*fireball_1TeV  );
    if((string)Name=="fireball_2TeV"  )  return (K==0 ? fireball_2TeV   : KFACTOR2*fireball_2TeV  );
    if((string)Name=="fireball_1.1TeV")  return (K==0 ? fireball_1p1TeV : KFACTOR2*fireball_1p1TeV);
    if((string)Name=="fireball_1.2TeV")  return (K==0 ? fireball_1p2TeV : KFACTOR2*fireball_1p2TeV);
    if((string)Name=="fireball_1.3TeV")  return (K==0 ? fireball_1p3TeV : KFACTOR2*fireball_1p3TeV);
    if((string)Name=="fireball_1.4TeV")  return (K==0 ? fireball_1p4TeV : KFACTOR2*fireball_1p4TeV);
    if((string)Name=="fireball_1.5TeV")  return (K==0 ? fireball_1p5TeV : KFACTOR2*fireball_1p5TeV);
    if((string)Name=="fireball_1.6TeV")  return (K==0 ? fireball_1p6TeV : KFACTOR2*fireball_1p6TeV);
    if((string)Name=="fireball_1.7TeV")  return (K==0 ? fireball_1p7TeV : KFACTOR2*fireball_1p7TeV);
    if((string)Name=="fireball_1.8TeV")  return (K==0 ? fireball_1p8TeV : KFACTOR2*fireball_1p8TeV);
    if((string)Name=="fireball_1.9TeV")  return (K==0 ? fireball_1p9TeV : KFACTOR2*fireball_1p9TeV);
	if((string)Name=="pptt") 	         return (K==0 ? pptt: KFACTOR2*pptt);
	if((string)Name=="ppvv")             return (K==0 ? ppvv: KFACTOR1*ppvv);
	if((string)Name=="ppvtt")            return (K==0 ? ppvtt: KFACTOR2*ppvtt);
	if((string)Name=="ppvvv")            return (K==0 ? ppvvv: KFACTOR1*ppvvv);
	if((string)Name=="ppvvtt")           return (K==0 ? ppvvtt: KFACTOR2*ppvvtt);
	if((string)Name=="pptttt")           return (K==0 ? pptttt: KFACTOR2*pptttt);
	if((string)Name=="ppvvvv")           return (K==0 ? ppvvvv: KFACTOR1*ppvvvv);
	if((string)Name=="ppvvvtt")          return (K==0 ? ppvvvtt: KFACTOR2*ppvvvtt);
}
double CrossSection::GetCrossSectionErr(const char* Name, int K){
    if((string)Name=="fireball_bp_1TeV") return Err_fireball_bp_1TeV;
    if((string)Name=="fireball_bp_2TeV") return Err_fireball_bp_2TeV;
    if((string)Name=="fireball_1TeV"  )  return (K==0 ? Err_fireball_1TeV   : KFACTOR2*Err_fireball_1TeV  );
    if((string)Name=="fireball_2TeV"  )  return (K==0 ? Err_fireball_2TeV   : KFACTOR2*Err_fireball_2TeV  );
    if((string)Name=="fireball_1.1TeV")  return (K==0 ? Err_fireball_1p1TeV : KFACTOR2*Err_fireball_1p1TeV);
    if((string)Name=="fireball_1.2TeV")  return (K==0 ? Err_fireball_1p2TeV : KFACTOR2*Err_fireball_1p2TeV);
    if((string)Name=="fireball_1.3TeV")  return (K==0 ? Err_fireball_1p3TeV : KFACTOR2*Err_fireball_1p3TeV);
    if((string)Name=="fireball_1.4TeV")  return (K==0 ? Err_fireball_1p4TeV : KFACTOR2*Err_fireball_1p4TeV);
    if((string)Name=="fireball_1.5TeV")  return (K==0 ? Err_fireball_1p5TeV : KFACTOR2*Err_fireball_1p5TeV);
    if((string)Name=="fireball_1.6TeV")  return (K==0 ? Err_fireball_1p6TeV : KFACTOR2*Err_fireball_1p6TeV);
    if((string)Name=="fireball_1.7TeV")  return (K==0 ? Err_fireball_1p7TeV : KFACTOR2*Err_fireball_1p7TeV);
    if((string)Name=="fireball_1.8TeV")  return (K==0 ? Err_fireball_1p8TeV : KFACTOR2*Err_fireball_1p8TeV);
    if((string)Name=="fireball_1.9TeV")  return (K==0 ? Err_fireball_1p9TeV : KFACTOR2*Err_fireball_1p9TeV);
	if((string)Name=="pptt") 	         return (K==0 ? Err_pptt: KFACTOR2*Err_pptt);
	if((string)Name=="ppvv")             return (K==0 ? Err_ppvv: KFACTOR1*Err_ppvv);
	if((string)Name=="ppvtt")            return (K==0 ? Err_ppvtt: KFACTOR2*Err_ppvtt);
	if((string)Name=="ppvvv")            return (K==0 ? Err_ppvvv: KFACTOR1*Err_ppvvv);
	if((string)Name=="ppvvtt")           return (K==1 ? Err_ppvvtt: KFACTOR2*Err_ppvvtt);
	if((string)Name=="pptttt")           return (K==0 ? Err_pptttt: KFACTOR2*Err_pptttt);
	if((string)Name=="ppvvvv")           return (K==0 ? Err_ppvvvv: KFACTOR1*Err_ppvvvv);
	if((string)Name=="ppvvvtt")          return (K==0 ? Err_ppvvvtt: KFACTOR2*Err_ppvvvtt);
}
int GetColor(const char* Name){
    if((string)Name=="fireball_bp_1TeV") return kBlue+1;
    if((string)Name=="fireball_bp_2TeV") return kRed;
	if((string)Name=="fireball_1TeV")    return kBlue+1;//kRed;//return kRed;
	if((string)Name=="fireball_2TeV")    return kRed;//kBlue+1;//kRed+2;//kOrange+8;//kPink-2;//kViolet+8;//return kRed;
    if((string)Name=="fireball_1.1TeV")  return kBlack;
    if((string)Name=="fireball_1.2TeV")  return kBlack;
    if((string)Name=="fireball_1.3TeV")  return kBlack;
    if((string)Name=="fireball_1.4TeV")  return kBlack;
    if((string)Name=="fireball_1.5TeV")  return kBlack;
    if((string)Name=="fireball_1.6TeV")  return kBlack;
    if((string)Name=="fireball_1.7TeV")  return kBlack;
    if((string)Name=="fireball_1.8TeV")  return kBlack;
    if((string)Name=="fireball_1.9TeV")  return kBlack;
	if((string)Name=="pptt") 			 return kPink-4;//return kOrange+3;
	if((string)Name=="ppvv")    		 return kViolet+6;//kOrange-2;//kOrange-4;//return kOrange+1;
	if((string)Name=="ppvtt")   		 return kYellow-7;//kYellow-7;//return kBlue;
	if((string)Name=="ppvvv")   		 return kGreen-3;//return kCyan;
	if((string)Name=="ppvvtt")  		 return kGreen+2;//return kGreen+2;
	if((string)Name=="pptttt")  		 return kViolet+6;//kViolet+8;//return kGreen;
	if((string)Name=="ppvvvv")  		 return kMagenta-6;//return kMagenta+2;
	if((string)Name=="ppvvvtt") 		 return kViolet-6;//return kMagenta;
}
const char* GetLegendName(const char* Name){
	if((string)Name=="fireball_bp_1TeV") return "fireball_bp_1TeV";
	if((string)Name=="fireball_bp_2TeV") return "fireball_bp_2TeV";
	if((string)Name=="fireball_1TeV") 	 return "fireball_bp_1TeV";
	if((string)Name=="fireball_2TeV") 	 return "fireball_bp_2TeV";
	if((string)Name=="pptt")             return "pp #rightarrow t#bar{t}";
	if((string)Name=="ppvtt")            return "pp #rightarrow vt#bar{t}";
	if((string)Name=="ppvvtt")           return "pp #rightarrow vvt#bar{t}";
	if((string)Name=="pptttt")           return "pp #rightarrow tt#bar{t}#bar{t}";
	if((string)Name=="ppvvvtt")          return "pp #rightarrow vvvt#bar{t}";
	if((string)Name=="ppvv")             return "pp #rightarrow vv";
	if((string)Name=="ppvvv")            return "pp #rightarrow vvv";
	if((string)Name=="ppvvvv")           return "pp #rightarrow vvvv";
}

//======================//
//=== SIGNAL HANDLER ===//
//======================//
void sighandler(int sig) {while(1);}
//### This is meant for finding segmentation error
//### by putting the below two lines to proper place
//std::cout<<"[INFO] SignalHandler"<<std::endl;
//signal(SIGSEGV,sighandler);
