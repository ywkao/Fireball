#include <stdio.h> //FILE
#include <iostream> //cout
#include <TCollection.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include "path.h"

const int Nbin = 20;
const double L = 100;
const double pb2fb	 = 1000;
const double pptt    = 815.96   ,Err_pptt    = 45.51;
const double ppvv    = 192.4    ,Err_ppvv    = 4.1;
const double ppttv   = 1.77     ,Err_ppttv   = 0.10;
const double ppvvv   = 0.621    ,Err_ppvvv   = 0.073;
const double ppttvv  = 0.0219   ,Err_ppttvv  = 0.0013;
const double pptttt  = 0.0164   ,Err_pptttt  = 9e-04;
const double ppvvvv  = 0.00405  ,Err_ppvvvv  = 9e-05;
const double ppttvvv = 0.000249 ,Err_ppttvvv = 3.4e-05;
//const char  *PREFIX	 = "/home/xiaokao/Desktop/work_ABCD/skimmed/simulation_delphes";
//const char  *SUFFIX	 = "skimmed.root";
const char  *oriTAG	 = "original";
const char  *TAG	 = "PTcut";
FILE *LOGFile;

int color[8] = {kOrange+3, kOrange+1, kBlue, kCyan, kGreen+2, kGreen, kMagenta+2, kMagenta};
using namespace std;

//###for finding segmentation error
//std::cout<<"[INFO] SignalHandler"<<std::endl;
//signal(SIGSEGV,sighandler);
void sighandler(int sig) {while(1);}

//#############################################
//### Structure & Classes for Storing Info
//#############################################
struct CrossSection {
	//std::pair<double,double> pptt = std::make_pair(815.96,45.51); WHY DOES IT NOT WORK??
	double GetCrossSection(const char* Name){
		if(Name=="pptt") 	return pptt;
		if(Name=="ppvv")    return ppvv;
		if(Name=="ppttv")   return ppttv;
		if(Name=="ppvvv")   return ppvvv;
		if(Name=="ppttvv")  return ppttvv;
		if(Name=="pptttt")  return pptttt;
		if(Name=="ppvvvv")  return ppvvvv;
		if(Name=="ppttvvv") return ppttvvv;
	}
	double GetCrossSectionErr(const char* Name){
		if(Name=="pptt") 	return Err_pptt;
		if(Name=="ppvv")    return Err_ppvv;
		if(Name=="ppttv")   return Err_ppttv;
		if(Name=="ppvvv")   return Err_ppvvv;
		if(Name=="ppttvv")  return Err_ppttvv;
		if(Name=="pptttt")  return Err_pptttt;
		if(Name=="ppvvvv")  return Err_ppvvvv;
		if(Name=="ppttvvv") return Err_ppttvvv;
	}
}; CrossSection X;
//-----
class MyYield {
public:
	double entryA, entryB, entryC, entryD, tot;//entries
	double A, B, C, D, Ex;//yield = entries * (L*X/N)
	double A2, B2, C2, D2;//yield = eff * (L*X)
	double errhA, errhB, errhC, errhD;
	double errlA, errlB, errlC, errlD;
	double biErrA, biErrB, biErrC, biErrD, biErrEx;//BINOMIAL ERROR: sqrt(N*eff*(1-eff))=sqrt(A*(1-eff))
	double pErrA, pErrB, pErrC, pErrD, pErrEx;//Poisson ERROR: sqrt(Yield)
};
class MySpectrum {
public:
	MySpectrum(){}
	const char *Name;
	TH1D *yA, *yB, *yC, *yD;//yield for each region
	TGraphAsymmErrors *gA, *gB, *gC, *gD;//efficiency for each region A,B,C,D
	MyYield yield;//calculate tot yield for each region
	~MySpectrum(){}
private:
	TH1D *A, *B, *C, *D, *ABCD;
	TH1D *oriA, *oriB, *oriC, *oriD, *oriABCD;
};

class MyProcess {
public:
	const char *Name;
	double factor;
	MySpectrum ST, HT, LPT, MET, NumJet;
	TH2D *NjvsNl, *HTvsMET;// HAVE TO BE SCALED
};
//-----
class MySpectrumCombined : public MySpectrum {
public:
	THStack *stackA, *stackB, *stackC, *stackD;
	TH1D *yExpected;
	TGraphAsymmErrors *gExpected;
};
class MyProcessCombined {
public:
	MySpectrumCombined ST, HT, LPT, MET, NumJet;
	TH2D *NjvsNl, *HTvsMET;// HAVE TO BE SCALED
};

//#################################################
//### Functions for individual & combined processes
//#################################################
void INFOLOG(const char *message){
	LOGFile = fopen("INFOLOG","a");
	fputs(Form("%s\n",message),LOGFile);
	fclose(LOGFile);
}
double GetMax(const char* Name){
	if(Name=="ST"	 ) return 5000;
	if(Name=="HT"	 ) return 5000;
	if(Name=="LPT"	 ) return 1200;
	if(Name=="MET"	 ) return 1200;
	if(Name=="NumJet") return 25;
}

const char* GetUnit(const char* Name){
	if(Name=="ST") return "ST [GeV]";
	if(Name=="HT") return "HT [GeV]";
	if(Name=="LPT") return "LPT [GeV]";
	if(Name=="MET") return "MET [GeV]";
	if(Name=="NumJet") return "# of Jets";
}

void RecordYield(const char* SpectrumName, MyYield yield){
	INFOLOG(Form("%s Entries = %.0f",SpectrumName,yield.tot));
	char *message; double errh, errl, biErr, pErr, ratio;

	message=Form("Entries A = %10.3f, yield A = %10.3f \u00B1 %10.5f \u00B1 %10.5f", yield.entryA, yield.A, yield.pErrA, yield.biErrA); INFOLOG(message);
	message=Form("Entries B = %10.3f, yield B = %10.3f \u00B1 %10.5f \u00B1 %10.5f", yield.entryB, yield.B, yield.pErrB, yield.biErrB); INFOLOG(message);
	message=Form("Entries C = %10.3f, yield C = %10.3f \u00B1 %10.5f \u00B1 %10.5f", yield.entryC, yield.C, yield.pErrC, yield.biErrC); INFOLOG(message);
	message=Form("Entries D = %10.3f, yield D = %10.3f \u00B1 %10.5f \u00B1 %10.5f", yield.entryD, yield.D, yield.pErrD, yield.biErrD); INFOLOG(message);
	message=Form("Expected yield                  = %10.3f \u00B1 %10.5f \u00B1 %10.5f", yield.Ex, yield.pErrEx, yield.biErrEx); INFOLOG(message);

	ratio = yield.A/yield.B;
	errh = sqrt( pow(yield.errhA/yield.B,2) + pow(yield.errhB*yield.A/pow(yield.B,2),2) );
	errl = sqrt( pow(yield.errlA/yield.B,2) + pow(yield.errlB*yield.A/pow(yield.B,2),2) );
	biErr= sqrt( pow(yield.biErrA/yield.B,2) + pow(yield.biErrB*yield.A/pow(yield.B,2),2) + pow(ratio,2)/yield.tot );
	pErr = sqrt( pow(yield.pErrA/yield.B,2) + pow(yield.pErrB*yield.A/pow(yield.B,2),2) );
	if(yield.B==0) {ratio=0; pErr=0; biErr=0;}
	message=Form("A/B = %7.3f \u00B1 %10.5f \u00B1 %10.5f", ratio, pErr, biErr); INFOLOG(message);

	ratio = yield.D/yield.C;
	errh = sqrt( pow(yield.errhD/yield.C,2) + pow(yield.errhC*yield.D/pow(yield.C,2),2) );
	errl = sqrt( pow(yield.errlD/yield.C,2) + pow(yield.errlC*yield.D/pow(yield.C,2),2) );
	biErr= sqrt( pow(yield.biErrD/yield.C,2) + pow(yield.biErrC*yield.D/pow(yield.C,2),2) + pow(ratio,2)/yield.tot );
	pErr = sqrt( pow(yield.pErrD/yield.C,2) + pow(yield.pErrC*yield.D/pow(yield.C,2),2) );
	if(yield.C==0) {ratio=0; pErr=0; biErr=0;}
	message=Form("D/C = %7.3f \u00B1 %10.5f \u00B1 %10.5f\n", ratio, pErr, biErr); INFOLOG(message);
}

//#############################################
//### Functions for importing bkg processes
//#############################################
void Initialization(const char* Name, MySpectrum &spectrum){
	spectrum.Name = Name;
	spectrum.yA = new TH1D(Form("hist%s_A",Name),"",Nbin,0,GetMax(Name)); spectrum.yA->Sumw2(); 
	spectrum.yB = new TH1D(Form("hist%s_B",Name),"",Nbin,0,GetMax(Name)); spectrum.yB->Sumw2();
	spectrum.yC = new TH1D(Form("hist%s_C",Name),"",Nbin,0,GetMax(Name)); spectrum.yC->Sumw2();
	spectrum.yD = new TH1D(Form("hist%s_D",Name),"",Nbin,0,GetMax(Name)); spectrum.yD->Sumw2();
	spectrum.ABCD = new TH1D(Form("hist%s_ABCD",Name),"",Nbin,0,GetMax(Name)); spectrum.ABCD->Sumw2();
	spectrum.oriABCD = new TH1D(Form("hist%s_oriABCD",Name),"",Nbin,0,GetMax(Name)); spectrum.oriABCD->Sumw2();
	spectrum.gA = new TGraphAsymmErrors(Nbin);
	spectrum.gB = new TGraphAsymmErrors(Nbin);
	spectrum.gC = new TGraphAsymmErrors(Nbin);
	spectrum.gD = new TGraphAsymmErrors(Nbin);
	spectrum.yield.tot=0.;
	spectrum.yield.entryA=0.; spectrum.yield.errhA=0.; spectrum.yield.errlA=0.; spectrum.yield.biErrA=0.; spectrum.yield.pErrA=0.;
	spectrum.yield.entryB=0.; spectrum.yield.errhB=0.; spectrum.yield.errlB=0.; spectrum.yield.biErrB=0.; spectrum.yield.pErrB=0.;
	spectrum.yield.entryC=0.; spectrum.yield.errhC=0.; spectrum.yield.errlC=0.; spectrum.yield.biErrC=0.; spectrum.yield.pErrC=0.;
	spectrum.yield.entryD=0.; spectrum.yield.errhD=0.; spectrum.yield.errlD=0.; spectrum.yield.biErrD=0.; spectrum.yield.pErrD=0.;
	spectrum.yield.A=0.; spectrum.yield.A2=0.;
	spectrum.yield.B=0.; spectrum.yield.B2=0.;
	spectrum.yield.C=0.; spectrum.yield.C2=0.;
	spectrum.yield.D=0.; spectrum.yield.D2=0.;
}

//void CalculateYieldAndBayesError(TGraphAsymmErrors *&graph, TH1D *&hist, double &yield, double &errl, double &errh, double factor){
void CalculateEfficiency(TGraphAsymmErrors *&graph, double &yield, double &errl, double &errh, double factor){
	for(int i=0; i < graph->GetN();i++){
		double x,y;
		graph->GetPoint(i,x,y);
		//graph->SetPoint(i,x,y*factor);
		//graph->SetPointError(i,graph->GetErrorXlow(i),graph->GetErrorXhigh(i),graph->GetErrorYlow(i)*factor,graph->GetErrorYhigh(i)*factor);
		//hist->SetBinContent(i+1, y*factor);
		yield += y*factor;
		errl  += pow(graph->GetErrorYlow(i)*factor,2);
		errh  += pow(graph->GetErrorYhigh(i)*factor,2);
	}
		errl  = sqrt(errl);
		errh  = sqrt(errh);
}
//void CalculateErrors(TH1D *hist, double &entries, double &yield, double &biErr, double &pErr, double N, double factor){
void CalculateErrors(TH1D *hist, TH1D *&histYield, double &entries, double &yield, double &biErr, double &pErr, double N, double factor){
	for(int bin=0; bin<hist->GetXaxis()->GetNbins(); bin++){
		entries += hist->GetBinContent(bin);
		histYield->SetBinContent(bin,hist->GetBinContent(bin)*factor/N);
		histYield->SetBinError(bin, sqrt(hist->GetBinContent(bin))*factor/N);
	}
	double eff = entries/N;
	biErr = sqrt(N*eff*(1-eff))*factor/N;
	pErr  = sqrt(entries)*factor/N;
	yield = entries*factor/N;
}

void CalculateEachRegion(const char* SpectrumName, MySpectrum &spectrum, double factor){
	Initialization(SpectrumName,spectrum);
	spectrum.oriABCD->Add(spectrum.oriA);
	spectrum.oriABCD->Add(spectrum.oriB);
	spectrum.oriABCD->Add(spectrum.oriC);
	spectrum.oriABCD->Add(spectrum.oriD);

	//###Calculate Poisson error & Binomial error
	//double tot = 0.; for(int bin=0; bin<spectrum.ABCD->GetXaxis()->GetNbins(); bin++) tot += spectrum.ABCD->GetBinContent(bin);
	for(int bin=0; bin<spectrum.oriABCD->GetXaxis()->GetNbins(); bin++) spectrum.yield.tot += spectrum.oriABCD->GetBinContent(bin);
	CalculateErrors(spectrum.A, spectrum.yA, spectrum.yield.entryA, spectrum.yield.A, spectrum.yield.biErrA, spectrum.yield.pErrA, spectrum.yield.tot, factor);
	CalculateErrors(spectrum.B, spectrum.yB, spectrum.yield.entryB, spectrum.yield.B, spectrum.yield.biErrB, spectrum.yield.pErrB, spectrum.yield.tot, factor);
	CalculateErrors(spectrum.C, spectrum.yC, spectrum.yield.entryC, spectrum.yield.C, spectrum.yield.biErrC, spectrum.yield.pErrC, spectrum.yield.tot, factor);
	CalculateErrors(spectrum.D, spectrum.yD, spectrum.yield.entryD, spectrum.yield.D, spectrum.yield.biErrD, spectrum.yield.pErrD, spectrum.yield.tot, factor);
	spectrum.yield.Ex = (spectrum.yield.C/spectrum.yield.B)*spectrum.yield.A;
	spectrum.yield.biErrEx = spectrum.yield.Ex
		*sqrt(pow(spectrum.yield.biErrA/spectrum.yield.A,2) + pow(spectrum.yield.biErrB/spectrum.yield.B,2) + pow(spectrum.yield.biErrC/spectrum.yield.C,2));
	spectrum.yield.pErrEx = spectrum.yield.Ex
		*sqrt(pow(spectrum.yield.pErrA/spectrum.yield.A,2) + pow(spectrum.yield.pErrB/spectrum.yield.B,2) + pow(spectrum.yield.pErrC/spectrum.yield.C,2));

	//###Calculate error with biult-in function
	spectrum.gA->BayesDivide(spectrum.A,spectrum.oriABCD);
	spectrum.gB->BayesDivide(spectrum.B,spectrum.oriABCD);
	spectrum.gC->BayesDivide(spectrum.C,spectrum.oriABCD);
	spectrum.gD->BayesDivide(spectrum.D,spectrum.oriABCD);
	CalculateEfficiency(spectrum.gA, spectrum.yield.A2, spectrum.yield.errlA, spectrum.yield.errhA, factor);
	CalculateEfficiency(spectrum.gB, spectrum.yield.B2, spectrum.yield.errlB, spectrum.yield.errhB, factor);
	CalculateEfficiency(spectrum.gC, spectrum.yield.C2, spectrum.yield.errlC, spectrum.yield.errhC, factor);
	CalculateEfficiency(spectrum.gD, spectrum.yield.D2, spectrum.yield.errlD, spectrum.yield.errhD, factor);
	//CalculateYieldAndBayesError(spectrum.gA, spectrum.yA, spectrum.yield.A2, spectrum.yield.errlA, spectrum.yield.errhA, factor);
	//CalculateYieldAndBayesError(spectrum.gB, spectrum.yB, spectrum.yield.B2, spectrum.yield.errlB, spectrum.yield.errhB, factor);
	//CalculateYieldAndBayesError(spectrum.gC, spectrum.yC, spectrum.yield.C2, spectrum.yield.errlC, spectrum.yield.errhC, factor);
	//CalculateYieldAndBayesError(spectrum.gD, spectrum.yD, spectrum.yield.D2, spectrum.yield.errlD, spectrum.yield.errhD, factor);
	//spectrum.gA->Divide(spectrum.A,spectrum.ABCD,"pois");
	//spectrum.gB->Divide(spectrum.B,spectrum.ABCD,"pois");
	//spectrum.gC->Divide(spectrum.C,spectrum.ABCD,"pois");
	//spectrum.gD->Divide(spectrum.D,spectrum.ABCD,"pois");

	RecordYield(SpectrumName, spectrum.yield);
}

void importEvents(std::vector<MyProcess> &vec, const char* NAME){
	const char *inputFile = Form("%s_%s_%s",PREFIX,NAME,SUFFIX);
	TFile *fin = new TFile(inputFile);

	MyProcess process;
	process.Name = NAME; INFOLOG(Form("### %s ###",NAME));
	process.factor = L*pb2fb*X.GetCrossSection(NAME); INFOLOG(Form("factor = %.0f",process.factor));
	//NEED TO CHECK IF SUMW2 WORKS
	process.ST.A = (TH1D*) fin->Get(Form("histST_%s01",TAG)); process.LPT.A = (TH1D*) fin->Get(Form("histLPT_%s01",TAG));
	process.ST.B = (TH1D*) fin->Get(Form("histST_%s00",TAG)); process.LPT.B = (TH1D*) fin->Get(Form("histLPT_%s00",TAG));
	process.ST.C = (TH1D*) fin->Get(Form("histST_%s10",TAG)); process.LPT.C = (TH1D*) fin->Get(Form("histLPT_%s10",TAG));
	process.ST.D = (TH1D*) fin->Get(Form("histST_%s11",TAG)); process.LPT.D = (TH1D*) fin->Get(Form("histLPT_%s11",TAG));
	process.HT.A = (TH1D*) fin->Get(Form("histHT_%s01",TAG)); process.MET.A = (TH1D*) fin->Get(Form("histMET_%s01",TAG));
	process.HT.B = (TH1D*) fin->Get(Form("histHT_%s00",TAG)); process.MET.B = (TH1D*) fin->Get(Form("histMET_%s00",TAG));
	process.HT.C = (TH1D*) fin->Get(Form("histHT_%s10",TAG)); process.MET.C = (TH1D*) fin->Get(Form("histMET_%s10",TAG));
	process.HT.D = (TH1D*) fin->Get(Form("histHT_%s11",TAG)); process.MET.D = (TH1D*) fin->Get(Form("histMET_%s11",TAG));
	process.NumJet.A = (TH1D*) fin->Get(Form("histNumJet_%s01",TAG)); process.NumJet.C = (TH1D*) fin->Get(Form("histNumJet_%s10",TAG));
	process.NumJet.B = (TH1D*) fin->Get(Form("histNumJet_%s00",TAG)); process.NumJet.D = (TH1D*) fin->Get(Form("histNumJet_%s11",TAG));
	process.NjvsNl  = (TH2D*) fin->Get(Form("hist_NjvsNl_%s",TAG));
	process.HTvsMET = (TH2D*) fin->Get(Form("hist_HTvsMET_%s",TAG));

	process.ST.oriA = (TH1D*) fin->Get(Form("histST_%s01",oriTAG)); process.LPT.oriA = (TH1D*) fin->Get(Form("histLPT_%s01",oriTAG));
	process.ST.oriB = (TH1D*) fin->Get(Form("histST_%s00",oriTAG)); process.LPT.oriB = (TH1D*) fin->Get(Form("histLPT_%s00",oriTAG));
	process.ST.oriC = (TH1D*) fin->Get(Form("histST_%s10",oriTAG)); process.LPT.oriC = (TH1D*) fin->Get(Form("histLPT_%s10",oriTAG));
	process.ST.oriD = (TH1D*) fin->Get(Form("histST_%s11",oriTAG)); process.LPT.oriD = (TH1D*) fin->Get(Form("histLPT_%s11",oriTAG));
	process.HT.oriA = (TH1D*) fin->Get(Form("histHT_%s01",oriTAG)); process.MET.oriA = (TH1D*) fin->Get(Form("histMET_%s01",oriTAG));
	process.HT.oriB = (TH1D*) fin->Get(Form("histHT_%s00",oriTAG)); process.MET.oriB = (TH1D*) fin->Get(Form("histMET_%s00",oriTAG));
	process.HT.oriC = (TH1D*) fin->Get(Form("histHT_%s10",oriTAG)); process.MET.oriC = (TH1D*) fin->Get(Form("histMET_%s10",oriTAG));
	process.HT.oriD = (TH1D*) fin->Get(Form("histHT_%s11",oriTAG)); process.MET.oriD = (TH1D*) fin->Get(Form("histMET_%s11",oriTAG));
	process.NumJet.oriA = (TH1D*) fin->Get(Form("histNumJet_%s01",oriTAG)); process.NumJet.oriC = (TH1D*) fin->Get(Form("histNumJet_%s10",oriTAG));
	process.NumJet.oriB = (TH1D*) fin->Get(Form("histNumJet_%s00",oriTAG)); process.NumJet.oriD = (TH1D*) fin->Get(Form("histNumJet_%s11",oriTAG));

	CalculateEachRegion("ST",process.ST,process.factor);
	CalculateEachRegion("HT",process.HT,process.factor);
	CalculateEachRegion("LPT",process.LPT,process.factor);
	CalculateEachRegion("MET",process.MET,process.factor);
	CalculateEachRegion("NumJet",process.NumJet,process.factor);

	vec.push_back(process);
}
//#############################################
//### Functions for Combining Processes
//#############################################
void InitializationCombined(const char* Name, MySpectrumCombined &spectrum){
	spectrum.Name = Name;
	spectrum.stackA = new THStack(Form("%s_A",Name),Form(";%s;Entries",GetUnit(Name)));
	spectrum.stackB = new THStack(Form("%s_B",Name),Form(";%s;Entries",GetUnit(Name)));
	spectrum.stackC = new THStack(Form("%s_C",Name),Form(";%s;Entries",GetUnit(Name)));
	spectrum.stackD = new THStack(Form("%s_D",Name),Form(";%s;Entries",GetUnit(Name)));

	spectrum.yA = new TH1D(Form("hist%s_Combined_A",Name),"",Nbin,0,GetMax(Name)); spectrum.yA->Sumw2();
	spectrum.yB = new TH1D(Form("hist%s_Combined_B",Name),"",Nbin,0,GetMax(Name)); spectrum.yB->Sumw2();
	spectrum.yC = new TH1D(Form("hist%s_Combined_C",Name),"",Nbin,0,GetMax(Name)); spectrum.yC->Sumw2();
	spectrum.yD = new TH1D(Form("hist%s_Combined_D",Name),"",Nbin,0,GetMax(Name)); spectrum.yD->Sumw2();
	spectrum.yExpected = new TH1D(Form("hist%s_Combined_Expected",Name),"",Nbin,0,GetMax(Name)); spectrum.yExpected->Sumw2();

	spectrum.gA = new TGraphAsymmErrors(Nbin);
	spectrum.gB = new TGraphAsymmErrors(Nbin);
	spectrum.gC = new TGraphAsymmErrors(Nbin);
	spectrum.gD = new TGraphAsymmErrors(Nbin);
	spectrum.gExpected = new TGraphAsymmErrors(Nbin);

	spectrum.yield.tot=0.;
	spectrum.yield.entryA=0.; spectrum.yield.errhA=0.; spectrum.yield.errlA=0.; spectrum.yield.biErrA=0.; spectrum.yield.pErrA=0.;
	spectrum.yield.entryB=0.; spectrum.yield.errhB=0.; spectrum.yield.errlB=0.; spectrum.yield.biErrB=0.; spectrum.yield.pErrB=0.;
	spectrum.yield.entryC=0.; spectrum.yield.errhC=0.; spectrum.yield.errlC=0.; spectrum.yield.biErrC=0.; spectrum.yield.pErrC=0.;
	spectrum.yield.entryD=0.; spectrum.yield.errhD=0.; spectrum.yield.errlD=0.; spectrum.yield.biErrD=0.; spectrum.yield.pErrD=0.;

	spectrum.yield.Ex=0.; spectrum.yield.biErrEx=0.; spectrum.yield.pErrEx=0.;

	spectrum.yield.A=0.; spectrum.yield.A2=0.;
	spectrum.yield.B=0.; spectrum.yield.B2=0.;
	spectrum.yield.C=0.; spectrum.yield.C2=0.;
	spectrum.yield.D=0.; spectrum.yield.D2=0.;
}
void HistCombineIndividualSpectrum(MySpectrum &itSpectrum, MySpectrumCombined &sumSpectrum){
	sumSpectrum.stackA->Add(itSpectrum.yA); sumSpectrum.yA->Add(itSpectrum.yA); sumSpectrum.gA->DoMerge(itSpectrum.gA);
	sumSpectrum.stackB->Add(itSpectrum.yB); sumSpectrum.yB->Add(itSpectrum.yB); sumSpectrum.gB->DoMerge(itSpectrum.gB);
	sumSpectrum.stackC->Add(itSpectrum.yC); sumSpectrum.yC->Add(itSpectrum.yC); sumSpectrum.gC->DoMerge(itSpectrum.gC);
	sumSpectrum.stackD->Add(itSpectrum.yD); sumSpectrum.yD->Add(itSpectrum.yD); sumSpectrum.gD->DoMerge(itSpectrum.gD);

	sumSpectrum.yield.tot+=itSpectrum.yield.tot;
	sumSpectrum.yield.A+=itSpectrum.yield.A; sumSpectrum.yield.errhA+=pow(itSpectrum.yield.errhA,2); sumSpectrum.yield.errlA+=pow(itSpectrum.yield.errlA,2);
	sumSpectrum.yield.B+=itSpectrum.yield.B; sumSpectrum.yield.errhB+=pow(itSpectrum.yield.errhB,2); sumSpectrum.yield.errlB+=pow(itSpectrum.yield.errlB,2);
	sumSpectrum.yield.C+=itSpectrum.yield.C; sumSpectrum.yield.errhC+=pow(itSpectrum.yield.errhC,2); sumSpectrum.yield.errlC+=pow(itSpectrum.yield.errlC,2);
	sumSpectrum.yield.D+=itSpectrum.yield.D; sumSpectrum.yield.errhD+=pow(itSpectrum.yield.errhD,2); sumSpectrum.yield.errlD+=pow(itSpectrum.yield.errlD,2);
	sumSpectrum.yield.Ex+=itSpectrum.yield.D; 

	sumSpectrum.yield.biErrA+=pow(itSpectrum.yield.biErrA,2); sumSpectrum.yield.pErrA+=pow(itSpectrum.yield.pErrA,2);
	sumSpectrum.yield.biErrB+=pow(itSpectrum.yield.biErrB,2); sumSpectrum.yield.pErrB+=pow(itSpectrum.yield.pErrB,2);
	sumSpectrum.yield.biErrC+=pow(itSpectrum.yield.biErrC,2); sumSpectrum.yield.pErrC+=pow(itSpectrum.yield.pErrC,2);
	sumSpectrum.yield.biErrD+=pow(itSpectrum.yield.biErrD,2); sumSpectrum.yield.pErrD+=pow(itSpectrum.yield.pErrD,2);
	sumSpectrum.yield.biErrEx+=pow(itSpectrum.yield.biErrEx,2); sumSpectrum.yield.pErrEx+=pow(itSpectrum.yield.pErrEx,2);

	sumSpectrum.yield.entryA+=itSpectrum.yield.entryA; sumSpectrum.yield.A2+=itSpectrum.yield.A2;
	sumSpectrum.yield.entryB+=itSpectrum.yield.entryB; sumSpectrum.yield.B2+=itSpectrum.yield.B2;
	sumSpectrum.yield.entryC+=itSpectrum.yield.entryC; sumSpectrum.yield.C2+=itSpectrum.yield.C2;
	sumSpectrum.yield.entryD+=itSpectrum.yield.entryD; sumSpectrum.yield.D2+=itSpectrum.yield.D2;
}
void HistCombinedErrorCalculation(MyYield &yield){
	yield.errhA = sqrt(yield.errhA); yield.errlA = sqrt(yield.errlA); yield.biErrA = sqrt(yield.biErrA); yield.pErrA = sqrt(yield.pErrA);
	yield.errhB = sqrt(yield.errhB); yield.errlB = sqrt(yield.errlB); yield.biErrB = sqrt(yield.biErrB); yield.pErrB = sqrt(yield.pErrB);
	yield.errhC = sqrt(yield.errhC); yield.errlC = sqrt(yield.errlC); yield.biErrC = sqrt(yield.biErrC); yield.pErrC = sqrt(yield.pErrC);
	yield.errhD = sqrt(yield.errhD); yield.errlD = sqrt(yield.errlD); yield.biErrD = sqrt(yield.biErrD); yield.pErrD = sqrt(yield.pErrD);
	yield.biErrEx = sqrt(yield.biErrEx); yield.pErrEx = sqrt(yield.pErrEx);
}
void SetExpectedHist(double bin, TH1D *&hist, TH1D *yA, double factor){
	hist->SetBinContent(bin, yA->GetBinContent(bin)*factor);
	hist->SetBinError(bin, yA->GetBinError(bin)*factor);
}
void HistCombined(std::vector<MyProcess> vec, MyProcessCombined &sum){
	InitializationCombined("ST", sum.ST);
	InitializationCombined("HT", sum.HT);
	InitializationCombined("LPT", sum.LPT);
	InitializationCombined("MET", sum.MET);
	InitializationCombined("NumJet", sum.NumJet);

	sum.NjvsNl  = new TH2D(Form("hist_NjvsNl_combined_%s",TAG) ,";# of lep;# of jet",7,0,7,25,0,25);                sum.NjvsNl ->Sumw2();
	sum.HTvsMET = new TH2D(Form("hist_HTvsMET_combined_%s",TAG),";MET [60 GeV]; HT [250 GeV]",20,0,1200,20,0,5000); sum.HTvsMET->Sumw2();
	
	for(std::vector<MyProcess>::iterator it=vec.begin(); it!=vec.end(); it++){
		HistCombineIndividualSpectrum(it->ST, sum.ST);
		HistCombineIndividualSpectrum(it->HT, sum.HT);
		HistCombineIndividualSpectrum(it->LPT, sum.LPT);
		HistCombineIndividualSpectrum(it->MET, sum.MET);
		HistCombineIndividualSpectrum(it->NumJet, sum.NumJet);
		sum.NjvsNl  -> Add(it->NjvsNl); 
        sum.HTvsMET -> Add(it->HTvsMET); 
	}
	HistCombinedErrorCalculation(sum.ST.yield);
	HistCombinedErrorCalculation(sum.HT.yield);
	HistCombinedErrorCalculation(sum.LPT.yield);
	HistCombinedErrorCalculation(sum.MET.yield);
	HistCombinedErrorCalculation(sum.NumJet.yield);

	INFOLOG("#############################");
	INFOLOG("### Combine all processes ###");
	INFOLOG("#############################");
	RecordYield("ST", sum.ST.yield);
	RecordYield("HT", sum.HT.yield);
	RecordYield("LPT", sum.LPT.yield);
	RecordYield("MET", sum.MET.yield);
	RecordYield("NumJet", sum.NumJet.yield);

	//###Expected data by simply mapping with const. factor: D = (C/B)*A
	//sum.ST.yExpected  = (TH1D*)sum.ST.yA ->Clone(); 
	//sum.HT.yExpected  = (TH1D*)sum.HT.yA ->Clone(); 
	//sum.LPT.yExpected = (TH1D*)sum.LPT.yA->Clone(); 
	//sum.MET.yExpected = (TH1D*)sum.MET.yA->Clone(); 
	//sum.NumJet.yExpected = (TH1D*)sum.NumJet.yA->Clone();

	for(int bin=0; bin<Nbin; bin++) SetExpectedHist(bin, sum.ST.yExpected , sum.ST.yA, sum.ST.yield.C/sum.ST.yield.B);
	for(int bin=0; bin<Nbin; bin++) SetExpectedHist(bin, sum.HT.yExpected , sum.HT.yA, sum.HT.yield.C/sum.HT.yield.B);
	for(int bin=0; bin<Nbin; bin++) SetExpectedHist(bin, sum.LPT.yExpected, sum.LPT.yA, sum.LPT.yield.C/sum.LPT.yield.B);
	for(int bin=0; bin<Nbin; bin++) SetExpectedHist(bin, sum.MET.yExpected, sum.MET.yA, sum.MET.yield.C/sum.MET.yield.B);
	for(int bin=0; bin<Nbin; bin++) SetExpectedHist(bin, sum.NumJet.yExpected, sum.NumJet.yA, sum.NumJet.yield.C/sum.NumJet.yield.B);

	//###Expected data with template fit
	//sum.ST.yExpected  = (TH1D*)sum.ST.yA ->Clone(); sum.ST.yExpected  -> Divide(sum.ST.yB);  sum.ST.yExpected  -> Multiply(sum.ST.yC);
	//sum.HT.yExpected  = (TH1D*)sum.HT.yA ->Clone(); sum.HT.yExpected  -> Divide(sum.HT.yB);  sum.HT.yExpected  -> Multiply(sum.HT.yC);
	//sum.LPT.yExpected = (TH1D*)sum.LPT.yA->Clone(); sum.LPT.yExpected -> Divide(sum.LPT.yB); sum.LPT.yExpected -> Multiply(sum.LPT.yC);
	//sum.MET.yExpected = (TH1D*)sum.MET.yA->Clone(); sum.MET.yExpected -> Divide(sum.MET.yB); sum.MET.yExpected -> Multiply(sum.MET.yC);
	//sum.NumJet.yExpected = (TH1D*)sum.NumJet.yA->Clone(); sum.NumJet.yExpected -> Divide(sum.NumJet.yB); sum.NumJet.yExpected -> Multiply(sum.NumJet.yC);

	//NEED TO COPE WITH GRAPH
}
//#############################################
//### Functions for Plotting
//#############################################
void DrawSettingTH1D(const char* Name, TH1D *hist, int LineColor, int FillColor, int FillStyle){
	hist->SetLineColor(LineColor);
	hist->SetFillColor(FillColor);
	hist->SetFillStyle(FillStyle);
	hist->SetOption("hist");
}
void DrawIndividualSpectrumSetting(const char* Name, MySpectrum spectrum, int color, TLegend* &legend){
	DrawSettingTH1D(Name,spectrum.yA,color,color,1001);
	DrawSettingTH1D(Name,spectrum.yB,color,color,1001);
	DrawSettingTH1D(Name,spectrum.yC,color,color,1001);
	DrawSettingTH1D(Name,spectrum.yD,color,color,1001);
	legend->AddEntry(spectrum.yD,Name,"f");
}
void DrawStackExecution(TCanvas* can, MySpectrumCombined spectrum){
	can->Divide(2,2); can->SetLogy();
	can->cd(1); spectrum.stackA->Draw();
	can->cd(3); spectrum.stackB->Draw();
	can->cd(4); spectrum.stackC->Draw();
	can->cd(2); spectrum.stackD->Draw();

	can->SaveAs(Form("result_%s_%s.pdf",spectrum.Name,TAG));
	//can->Clear();
}
void DrawHistExecution(TCanvas* can, MySpectrumCombined spectrum){
	can->Divide(2,2); can->SetLogy();
	can->cd(1); spectrum.yA->Draw();
	can->cd(3); spectrum.yB->Draw();
	can->cd(4); spectrum.yC->Draw();
	can->cd(2); spectrum.yD->Draw();

	int max_bin = spectrum.yExpected->GetMaximumBin();
	double max  = spectrum.yExpected->GetBinContent(max_bin)+spectrum.yExpected->GetBinError(max_bin)+100;
	spectrum.yD->SetMaximum(max);

	spectrum.yExpected->SetLineColor(kRed);
	spectrum.yExpected->Draw("hist,same");

	can->SaveAs(Form("result_%s_%s.pdf",spectrum.Name,TAG));
	//can->Clear();
}
//###
void DrawExpectedSpectrum(TCanvas *can, MySpectrumCombined spectrum, TLegend* &legend, double Max){
	gStyle->SetOptStat(0); //can->SetLogy();
	TPad *pad1 = new TPad("pad1","",0,0.2,1,1);
	TPad *pad2 = new TPad("pad2","",0,0.05,1,0.25);
	pad1->Draw(); 
	pad2->Draw(); 
	pad1->cd();
	spectrum.stackD->Draw("hist");
	spectrum.yD->SetFillStyle(3004);
	spectrum.yD->SetFillColor(kBlue-4);
	spectrum.yD->Draw("E2,same");

	int max_bin = spectrum.yExpected->GetMaximumBin();
	double max  = spectrum.yExpected->GetBinContent(max_bin)+spectrum.yExpected->GetBinError(max_bin)+100;
	spectrum.stackD->SetMaximum(max);
	spectrum.yExpected->SetLineColor(kRed);
	spectrum.yExpected->Draw("same");
	legend->Draw("same");

	pad2->cd();
	double ratio, err;
	TLine *line = new TLine(0,1,Max,1);
	TH1D *hist = (TH1D*) spectrum.yD->Clone();
	for(int bin=0; bin<Nbin; bin++){
		if(spectrum.yExpected->GetBinContent(bin)==0){ratio=0; err=0;}
		else{
		ratio=spectrum.yD->GetBinContent(bin)/spectrum.yExpected->GetBinContent(bin);
		err = ratio*sqrt(
		pow(spectrum.yD->GetBinError(bin)/spectrum.yD->GetBinContent(bin),2) + pow(spectrum.yExpected->GetBinError(bin)/spectrum.yExpected->GetBinContent(bin),2));
		}
		hist->SetBinContent(bin,ratio);
		hist->SetBinError(bin,err);
	}
	hist->SetFillStyle(0);
	hist->SetFillColor(kBlue);
	hist->SetMaximum(1.5);
	hist->SetMinimum(-0.2);
	hist->GetYaxis()->SetTitle("#frac{Expected}{Original}");
	hist->GetYaxis()->SetTitleSize(0.12);
	hist->GetYaxis()->SetTitleOffset(0.25);
	hist->GetYaxis()->CenterTitle();
	hist->Draw();
	line->SetLineColor(2);
	line->Draw("same");
	
	can->SaveAs(Form("result_%s_%s.pdf",spectrum.Name,TAG));
}
//###
const char* Process[8]={"pptt", "ppvv", "ppvtt", "ppvvv", "ppttvv", "pptttt", "ppvvvv", "ppttvvv"};
void DrawHistSpectrum(TCanvas *can, const char* Option, const char* Name, std::vector<MyProcess> &vec, MyProcessCombined &sum, TLegend* &legend){
	if(Option=="ABCD"){
		if(Name=="ST")  DrawHistExecution(can, sum.ST);
		if(Name=="HT")  DrawHistExecution(can, sum.HT);
		if(Name=="LPT") DrawHistExecution(can, sum.LPT);
		if(Name=="MET") DrawHistExecution(can, sum.MET);
		if(Name=="NumJet") DrawHistExecution(can, sum.NumJet);
	}
	//TLegend *legend = new TLegend(0.74,0.5,0.9,0.9);
	
	if(Option=="Stack" || Option=="Expected"){
		for(int i=0; i<8; i++){
			if(Name=="ST")  DrawIndividualSpectrumSetting(Process[i],vec.at(7-i).ST,color[i],legend);  
			if(Name=="HT")  DrawIndividualSpectrumSetting(Process[i],vec.at(7-i).HT,color[i],legend);  
			if(Name=="LPT") DrawIndividualSpectrumSetting(Process[i],vec.at(7-i).LPT,color[i],legend); 
			if(Name=="MET") DrawIndividualSpectrumSetting(Process[i],vec.at(7-i).MET,color[i],legend); 
			if(Name=="NumJet") DrawIndividualSpectrumSetting(Process[i],vec.at(7-i).NumJet,color[i],legend); 
		}
	}
	if(Option=="Stack"){
		if(Name=="ST")  DrawStackExecution(can, sum.ST);
		if(Name=="HT")  DrawStackExecution(can, sum.HT);
		if(Name=="LPT") DrawStackExecution(can, sum.LPT);
		if(Name=="MET") DrawStackExecution(can, sum.MET);
		if(Name=="NumJet") DrawStackExecution(can, sum.NumJet);
	}
	if(Option=="Expected"){
		if(Name=="ST")  DrawExpectedSpectrum(can, sum.ST,  legend, GetMax(Name));
		if(Name=="HT")  DrawExpectedSpectrum(can, sum.HT,  legend, GetMax(Name));
		if(Name=="LPT") DrawExpectedSpectrum(can, sum.LPT, legend, GetMax(Name));
		if(Name=="MET") DrawExpectedSpectrum(can, sum.MET, legend, GetMax(Name));
		if(Name=="NumJet") DrawExpectedSpectrum(can, sum.NumJet, legend, GetMax(Name));
	}
}
//###
void DrawMap(TCanvas *can, MyProcessCombined sum){
	can->SetLogy(0);

	TLine *lineH1 = new TLine(0,6,7,6);
	TLine *lineV1 = new TLine(2,0,2,25);
	sum.NjvsNl ->Draw("colz"); 
	lineH1 -> Draw("same"); 
	lineV1 -> Draw("same");
	can->SaveAs(Form("result_NjvsNl_%s.pdf",TAG));
	//can->Clear();

	//TLine *lineH2 = new TLine(0,1400,1200,1400);
	//TLine *lineV2 = new TLine(50,0,50,5000);
	//sum.HTvsMET->Draw("colz");
	//lineH2 -> Draw("same"); 
	//lineV2 -> Draw("same");
	//can->SaveAs(Form("result_HTvsMET_%s.pdf",TAG));
	//can->Clear();
}
//#############################################
//### Functions for Quick Check
//#############################################
void quickCheck(MySpectrum spectrum){
	TCanvas *c5 = new TCanvas("c5","Original",0,0,650,350);
	c5->Divide(2,2);
	c5->cd(1); spectrum.A->Draw();
	c5->cd(2); spectrum.D->Draw();
	c5->cd(3); spectrum.B->Draw();
	c5->cd(4); spectrum.C->Draw();

	TCanvas *c6 = new TCanvas("c6","Original(A+B+C+D)",650,0,650,350);
	spectrum.ABCD->Draw();

	TCanvas *c7 = new TCanvas("c7","Efficiency",0,350,650,350);
	c7->Divide(2,2);
	c7->cd(1); spectrum.gA->Draw();
	c7->cd(2); spectrum.gD->Draw();
	c7->cd(3); spectrum.gB->Draw();
	c7->cd(4); spectrum.gC->Draw();

	TCanvas *c8 = new TCanvas("c8","Yield",650,350,650,350);
	c8->Divide(2,2);
	c8->cd(1); spectrum.yA->Draw();
	c8->cd(2); spectrum.yD->Draw();
	c8->cd(3); spectrum.yB->Draw();
	c8->cd(4); spectrum.yC->Draw();
}
//void DrawStackSpectrum(TCanvas *can, const char* Name, std::vector<MyProcess> &vec, MyProcessCombined &sum){
//	for(int i=0; i<8; i++){
//		if(Name=="ST")  DrawIndividualSpectrumSetting(vec.at(7-i).ST,color[i]);  
//		if(Name=="HT")  DrawIndividualSpectrumSetting(vec.at(7-i).HT,color[i]);  
//		if(Name=="LPT") DrawIndividualSpectrumSetting(vec.at(7-i).LPT,color[i]); 
//		if(Name=="MET") DrawIndividualSpectrumSetting(vec.at(7-i).MET,color[i]); 
//		if(Name=="NumJet") DrawIndividualSpectrumSetting(vec.at(7-i).NumJet,color[i]); 
//	}
//	if(Name=="ST")  DrawStackExecution(can, sum.ST);
//	if(Name=="HT")  DrawStackExecution(can, sum.HT);
//	if(Name=="LPT") DrawStackExecution(can, sum.LPT);
//	if(Name=="MET") DrawStackExecution(can, sum.MET);
//	if(Name=="NumJet") DrawStackExecution(can, sum.NumJet);
//}
