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
#include "xsec.h"
using namespace std;

const int NUM = 9;
const double L = 100, pb2fb = 1000;
const double Err_L = 0.05;// deltaL/L
const double KFACTOR1=1.56;//ppvv-like process
const double KFACTOR2=1.69;//pptt-like process
const double Err_KFACTOR1=0.03;//ppvv-like process
const double Err_KFACTOR2=0.10;//pptt-like process
const char *histName[NUM] = {"Num_lep","Num_jet","PT_lep","PT_jet","LPT","HT","MET","ST","Num_boson"};
const char *units[NUM]    = {"# of lep", "# of jet", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "# of boson"};
const char* PREFIX = "~/Desktop/workspace/skimmed_new/result";
const char* SUFFIX = "root";
const char* TARGET;
int SIG, XERR=1;

double TotalYield, Err_ExpectedSigYield, MeasuredCrossSection, Err_MeasuredCrossSection, significance, p_value;
double MeasuredCrossSection_up, MeasuredCrossSection_down;
double First, Err_First, Second, Err_Second;
double yield_sig=0, yield_bg=0, Err_yield_bg=0, N=0, NB=0;
double ybg1, ybg2, ybg3, ybg4;

class MyProcess{
public:
	const char* Name;
	double Factor, Entries, X, Err_X, X_ori, Err_X_ori, KFactor, Err_KFactor;
	double Yield[3][NUM], Err_Yield[3][NUM], Eff[3][NUM], Err_Eff[3][NUM];
	TH1D *hist_Ori[NUM];
	TH1D *hist_Cut[NUM];
	TH1D *hist_Nm1[NUM];
    void Init(MyProcess process){
	    process.Factor    = 0;
        process.Entries   = 0;
        process.X         = 0;
        process.Err_X     = 0;
        process.X_ori     = 0;
        process.Err_X_ori = 0;
        process.KFactor   = 0;
        for(int i=0; i<3; i++){ for(int j=0; j<NUM; j++){
	        process.Yield[i][j]     = 0;
            process.Err_Yield[i][j] = 0;
            process.Eff[i][j]       = 0;
            process.Err_Eff[i][j]   = 0;
        }}
    }
};

struct CrossSection {
	double GetCrossSection(const char* Name, int K){
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
		if((string)Name=="pptt") 	return (K==0 ? pptt: KFACTOR2*pptt);
		if((string)Name=="ppvv")    return (K==0 ? ppvv: KFACTOR1*ppvv);
		if((string)Name=="ppvtt")   return (K==0 ? ppvtt: KFACTOR2*ppvtt);
		if((string)Name=="ppvvv")   return (K==0 ? ppvvv: KFACTOR1*ppvvv);
		if((string)Name=="ppvvtt")  return (K==0 ? ppvvtt: KFACTOR2*ppvvtt);
		if((string)Name=="pptttt")  return (K==0 ? pptttt: KFACTOR2*pptttt);
		if((string)Name=="ppvvvv")  return (K==0 ? ppvvvv: KFACTOR1*ppvvvv);
		if((string)Name=="ppvvvtt") return (K==0 ? ppvvvtt: KFACTOR2*ppvvvtt);
	}
	double GetCrossSectionErr(const char* Name, int K){
        if((string)Name=="fireball_bp_1TeV") return Err_fireball_bp_1TeV;
        if((string)Name=="fireball_bp_2TeV") return Err_fireball_bp_2TeV;
        if((string)Name=="fireball_1TeV"  ) return (K==0 ? Err_fireball_1TeV   : sqrt(pow(Err_KFACTOR2*fireball_1TeV  ,2)+pow(KFACTOR2*Err_fireball_1TeV  ,2)));
        if((string)Name=="fireball_2TeV"  ) return (K==0 ? Err_fireball_2TeV   : sqrt(pow(Err_KFACTOR2*fireball_2TeV  ,2)+pow(KFACTOR2*Err_fireball_2TeV  ,2)));
        if((string)Name=="fireball_1.1TeV") return (K==0 ? Err_fireball_1p1TeV : sqrt(pow(Err_KFACTOR2*fireball_1p1TeV,2)+pow(KFACTOR2*Err_fireball_1p1TeV,2)));
        if((string)Name=="fireball_1.2TeV") return (K==0 ? Err_fireball_1p2TeV : sqrt(pow(Err_KFACTOR2*fireball_1p2TeV,2)+pow(KFACTOR2*Err_fireball_1p2TeV,2)));
        if((string)Name=="fireball_1.3TeV") return (K==0 ? Err_fireball_1p3TeV : sqrt(pow(Err_KFACTOR2*fireball_1p3TeV,2)+pow(KFACTOR2*Err_fireball_1p3TeV,2)));
        if((string)Name=="fireball_1.4TeV") return (K==0 ? Err_fireball_1p4TeV : sqrt(pow(Err_KFACTOR2*fireball_1p4TeV,2)+pow(KFACTOR2*Err_fireball_1p4TeV,2)));
        if((string)Name=="fireball_1.5TeV") return (K==0 ? Err_fireball_1p5TeV : sqrt(pow(Err_KFACTOR2*fireball_1p5TeV,2)+pow(KFACTOR2*Err_fireball_1p5TeV,2)));
        if((string)Name=="fireball_1.6TeV") return (K==0 ? Err_fireball_1p6TeV : sqrt(pow(Err_KFACTOR2*fireball_1p6TeV,2)+pow(KFACTOR2*Err_fireball_1p6TeV,2)));
        if((string)Name=="fireball_1.7TeV") return (K==0 ? Err_fireball_1p7TeV : sqrt(pow(Err_KFACTOR2*fireball_1p7TeV,2)+pow(KFACTOR2*Err_fireball_1p7TeV,2)));
        if((string)Name=="fireball_1.8TeV") return (K==0 ? Err_fireball_1p8TeV : sqrt(pow(Err_KFACTOR2*fireball_1p8TeV,2)+pow(KFACTOR2*Err_fireball_1p8TeV,2)));
        if((string)Name=="fireball_1.9TeV") return (K==0 ? Err_fireball_1p9TeV : sqrt(pow(Err_KFACTOR2*fireball_1p9TeV,2)+pow(KFACTOR2*Err_fireball_1p9TeV,2)));
        if((string)Name=="pptt")            return (K==0 ? Err_pptt            : sqrt(pow(Err_KFACTOR2*pptt           ,2)+pow(KFACTOR2*Err_pptt           ,2)));
        if((string)Name=="ppvv")            return (K==0 ? Err_ppvv            : sqrt(pow(Err_KFACTOR1*ppvv           ,2)+pow(KFACTOR1*Err_ppvv           ,2)));
        if((string)Name=="ppvtt")           return (K==0 ? Err_ppvtt           : sqrt(pow(Err_KFACTOR2*ppvtt          ,2)+pow(KFACTOR2*Err_ppvtt          ,2)));
        if((string)Name=="ppvvv")           return (K==0 ? Err_ppvvv           : sqrt(pow(Err_KFACTOR1*ppvvv          ,2)+pow(KFACTOR1*Err_ppvvv          ,2)));
        if((string)Name=="ppvvtt")          return (K==0 ? Err_ppvvtt          : sqrt(pow(Err_KFACTOR2*ppvvtt         ,2)+pow(KFACTOR2*Err_ppvvtt         ,2)));
        if((string)Name=="pptttt")          return (K==0 ? Err_pptttt          : sqrt(pow(Err_KFACTOR2*pptttt         ,2)+pow(KFACTOR2*Err_pptttt         ,2)));
        if((string)Name=="ppvvvv")          return (K==0 ? Err_ppvvvv          : sqrt(pow(Err_KFACTOR1*ppvvvv         ,2)+pow(KFACTOR1*Err_ppvvvv         ,2)));
        if((string)Name=="ppvvvtt")         return (K==0 ? Err_ppvvvtt         : sqrt(pow(Err_KFACTOR2*ppvvvtt        ,2)+pow(KFACTOR2*Err_ppvvvtt        ,2)));
	}
}; CrossSection Xsec;
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
	if((string)Name=="fireball_bp_1TeV") 	return "fireball_bp_1TeV";
	if((string)Name=="fireball_bp_2TeV") 	return "fireball_bp_2TeV";
	if((string)Name=="fireball_1TeV") 	return "fireball_bp_1TeV";
	if((string)Name=="fireball_2TeV") 	return "fireball_bp_2TeV";
	if((string)Name=="pptt")    return "pp #rightarrow t#bar{t}";
	if((string)Name=="ppvtt")   return "pp #rightarrow vt#bar{t}";
	if((string)Name=="ppvvtt")  return "pp #rightarrow vvt#bar{t}";
	if((string)Name=="pptttt")  return "pp #rightarrow tt#bar{t}#bar{t}";
	if((string)Name=="ppvvvtt") return "pp #rightarrow vvvt#bar{t}";
	if((string)Name=="ppvv")    return "pp #rightarrow vv";
	if((string)Name=="ppvvv")   return "pp #rightarrow vvv";
	if((string)Name=="ppvvvv")  return "pp #rightarrow vvvv";
}

void ScaleAndDrawSetting(const char* TAG, TH1D *&hist, double &Yield, double &Err_Yield, double &Eff, double &Err_Eff, 
                         double entries, double factor, double X, double Err_X, int lcolor, int lstyle, int lwidth, int fcolor, int style){
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
    Err_Yield = sqrt(pow(Err_L,2) + pow(Err_X/X,2) + pow(err_yield/yield,2))*yield;
    Eff       = yield/(factor*entries);
    Err_Eff   = sqrt(Eff*(1-Eff)/entries);
	hist->SetLineColor(lcolor);
	hist->SetLineStyle(lstyle);
	hist->SetLineWidth(lwidth);
	hist->SetFillColor(fcolor);
	hist->SetFillStyle(style);
}

std::vector<MyProcess> vec, vec_bkg, vec_combine;
MyProcess Combined_ppvvv, Combined_ppttv;

void importHist(std::vector<MyProcess> &vec, const char *NAME){
	TARGET = Form("%s_%s.%s", PREFIX, NAME, SUFFIX);
	TFile *fin = new TFile(TARGET,"READ");
	MyProcess process;
	process.Name        =  NAME;
	process.Entries     =  ( (TH1D*)fin->Get(Form("Ori_%s",histName[0])) )->GetEntries();//get the # of events of the process
	process.X           =  Xsec.GetCrossSection(NAME,1);
	process.Err_X       =  Xsec.GetCrossSectionErr(NAME,1);
	process.X_ori       =  Xsec.GetCrossSection(NAME,0);   // w/o applying k-factor
	process.Err_X_ori   =  Xsec.GetCrossSectionErr(NAME,0);// w/o applying k-factor
    process.KFactor     =  process.X/process.X_ori;
    process.Err_KFactor =  (process.KFactor==KFACTOR1) ? Err_KFACTOR1 : Err_KFACTOR2;
	process.Factor      =  L*pb2fb*process.X/process.Entries;
	for(int i=0; i<3; i++) for(int j=0; j<NUM; j++) process.Yield[i][j]     = 0;
	for(int i=0; i<3; i++) for(int j=0; j<NUM; j++) process.Err_Yield[i][j] = 0;

	int lcolor, lstyle, lwidth, style;
	if      ((string)NAME=="fireball_1TeV") {lcolor=GetColor(NAME); lstyle=1; lwidth=2; style=3004;}
    else if ((string)NAME=="fireball_2TeV") {lcolor=GetColor(NAME); lstyle=1; lwidth=2; style=3004;}
	else {lcolor=GetColor(NAME); lstyle=0; lwidth=0; style=1001;}

	for(int k=0; k<NUM; k++){
		process.hist_Ori[k] = (TH1D*) fin->Get(Form("Ori_%s",histName[k]));
		process.hist_Cut[k] = (TH1D*) fin->Get(Form("Cut_%s",histName[k]));
		process.hist_Nm1[k] = (TH1D*) fin->Get(Form("N-1_%s",histName[k]));
		ScaleAndDrawSetting("Ori",process.hist_Ori[k], process.Yield[0][k], process.Err_Yield[0][k], process.Eff[0][k], process.Err_Eff[0][k],
                            process.Entries, process.Factor, process.X, process.Err_X, lcolor, lstyle, lwidth, GetColor(NAME), style);
		ScaleAndDrawSetting("Nm1",process.hist_Nm1[k], process.Yield[1][k], process.Err_Yield[1][k], process.Eff[1][k], process.Err_Eff[1][k],
                            process.Entries, process.Factor, process.X, process.Err_X, lcolor, lstyle, lwidth, GetColor(NAME), style);
		ScaleAndDrawSetting("Cut",process.hist_Cut[k], process.Yield[2][k], process.Err_Yield[2][k], process.Eff[2][k], process.Err_Eff[2][k],
                            process.Entries, process.Factor, process.X, process.Err_X, lcolor, lstyle, lwidth, GetColor(NAME), style);
	}

	vec.push_back(process);
}

//##################################//
//### calculate the combined bkg ###//
//##################################//
void CombineBackgroundProcess(std::vector<MyProcess> vec_bkg){
    Combined_ppvvv.Name="Combined_ppvvv"; Combined_ppttv.Name="Combined_ppttv";
    Combined_ppvvv.Init(Combined_ppvvv); Combined_ppttv.Init(Combined_ppttv);
	for(unsigned int i=0; i<vec_bkg.size(); i++){
        if(i==3||i==6)            {Combined_ppvvv.Entries         += vec_bkg.at(i).Entries;
                                   Combined_ppvvv.X               += vec_bkg.at(i).X;
                                   Combined_ppvvv.Err_X           += pow(vec_bkg.at(i).Err_X,2);
                                   Combined_ppvvv.Yield[2][0]     += vec_bkg.at(i).Yield[2][0];
                                   Combined_ppvvv.Err_Yield[2][0] += pow(vec_bkg.at(i).Err_Yield[2][0],2);}
        if(i==2||i==4||i==5||i==7){Combined_ppttv.Entries         += vec_bkg.at(i).Entries;
                                   Combined_ppttv.X               += vec_bkg.at(i).X;
                                   Combined_ppttv.Err_X           += pow(vec_bkg.at(i).Err_X,2);
                                   Combined_ppttv.Yield[2][0]     += vec_bkg.at(i).Yield[2][0];
                                   Combined_ppttv.Err_Yield[2][0] += pow(vec_bkg.at(i).Err_Yield[2][0],2);}
        if(i==vec_bkg.size()-1){Combined_ppvvv.Err_X = sqrt(Combined_ppvvv.Err_X); Combined_ppvvv.Err_Yield[2][0] = sqrt(Combined_ppvvv.Err_Yield[2][0]);
                                Combined_ppttv.Err_X = sqrt(Combined_ppttv.Err_X); Combined_ppttv.Err_Yield[2][0] = sqrt(Combined_ppttv.Err_Yield[2][0]);
                                vec_combine.push_back(Combined_ppttv); vec_combine.push_back(Combined_ppvvv);}
	}
    for(int i=0; i<2; i++){
        vec_combine.at(i).Eff[2][0]    = vec_combine.at(i).Yield[2][0]/(L*pb2fb*vec_combine.at(i).X);
        vec_combine.at(i).Err_Eff[2][0]= sqrt( vec_combine.at(i).Eff[2][0] * (1-vec_combine.at(i).Eff[2][0]) / vec_combine.at(i).Entries );
        //printf("%16s: X-sec = %7.4f \u00B1 %8.6f (%6.4f), yield = %5.2f \u00B1 %5.2f (%6.4f), Eff = %8.7f \u00B1 %8.7f (%8.7f)\n",
        //        vec_combine.at(i).Name, vec_combine.at(i).X, vec_combine.at(i).Err_X, vec_combine.at(i).Err_X/vec_combine.at(i).X, 
        //        vec_combine.at(i).Yield[2][0], vec_combine.at(i).Err_Yield[2][0], vec_combine.at(i).Err_Yield[2][0]/vec_combine.at(i).Yield[2][0],
        //        vec_combine.at(i).Eff[2][0], vec_combine.at(i).Err_Eff[2][0], vec_combine.at(i).Err_Eff[2][0]/vec_combine.at(i).Eff[2][0]);
    }
    ybg1 = vec_bkg.at(0).Yield[2][0];
    ybg2 = vec_bkg.at(1).Yield[2][0];
    ybg3 = vec_combine.at(0).Yield[2][0];
    ybg4 = vec_combine.at(1).Yield[2][0];
    yield_bg = ybg1+ybg2+ybg3+ybg4;
    for(int i=0; i<2; i++){
        Err_yield_bg += pow(vec_bkg.at(i).Err_Yield[2][0],2); 
        Err_yield_bg += pow(vec_combine.at(i).Err_Yield[2][0],2); 
        if(i==1) Err_yield_bg = sqrt(Err_yield_bg);
    }
    //printf("%5.2f \u00B1 %5.2f = %5.2f + %5.2f + %5.2f + %5.2f\n\n", yield_bg, Err_yield_bg, ybg1, ybg2, ybg3, ybg4);
}
