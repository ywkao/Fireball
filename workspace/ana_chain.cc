#include <stdio.h>
#include <iostream>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

using namespace std;

//const char *PREFIX = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/ppvv/Events";
//const char *SUFFIX = "tag_1_delphes_events.root";
const char *PREFIX = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/skimmed";
const char *oriTAG = "original";
const char *TAG	   = "PTcut";

//-----
class MySpectrum {
public:
	MySpectrum(){}
	TH1D *A, *B, *C, *D, *ABCD;
	TH1D *oriA, *oriB, *oriC, *oriD, *oriABCD;
	~MySpectrum(){}
};
class MyProcess {
public:
	MySpectrum ST, HT, LPT, MET, NumJet;
	TH2D *NjvsNl, *HTvsMET;// HAVE TO BE SCALED
};
//-----

void importEvents(MyProcess &sum, const char* NAME){
	const char *inputFile = Form("%s/%s",PREFIX,NAME);
	TFile *fin = new TFile(inputFile);

	//NEED TO CHECK IF SUMW2 WORKS
	sum.ST.A = (TH1D*) fin->Get(Form("histST_%s01",TAG)); sum.LPT.A = (TH1D*) fin->Get(Form("histLPT_%s01",TAG));
	sum.ST.B = (TH1D*) fin->Get(Form("histST_%s00",TAG)); sum.LPT.B = (TH1D*) fin->Get(Form("histLPT_%s00",TAG));
	sum.ST.C = (TH1D*) fin->Get(Form("histST_%s10",TAG)); sum.LPT.C = (TH1D*) fin->Get(Form("histLPT_%s10",TAG));
	sum.ST.D = (TH1D*) fin->Get(Form("histST_%s11",TAG)); sum.LPT.D = (TH1D*) fin->Get(Form("histLPT_%s11",TAG));
	sum.HT.A = (TH1D*) fin->Get(Form("histHT_%s01",TAG)); sum.MET.A = (TH1D*) fin->Get(Form("histMET_%s01",TAG));
	sum.HT.B = (TH1D*) fin->Get(Form("histHT_%s00",TAG)); sum.MET.B = (TH1D*) fin->Get(Form("histMET_%s00",TAG));
	sum.HT.C = (TH1D*) fin->Get(Form("histHT_%s10",TAG)); sum.MET.C = (TH1D*) fin->Get(Form("histMET_%s10",TAG));
	sum.HT.D = (TH1D*) fin->Get(Form("histHT_%s11",TAG)); sum.MET.D = (TH1D*) fin->Get(Form("histMET_%s11",TAG));
	sum.NumJet.A = (TH1D*) fin->Get(Form("histNumJet_%s01",TAG)); sum.NumJet.C = (TH1D*) fin->Get(Form("histNumJet_%s10",TAG));
	sum.NumJet.B = (TH1D*) fin->Get(Form("histNumJet_%s00",TAG)); sum.NumJet.D = (TH1D*) fin->Get(Form("histNumJet_%s11",TAG));
	sum.NjvsNl  = (TH2D*) fin->Get(Form("hist_NjvsNl_%s",TAG));
	sum.HTvsMET = (TH2D*) fin->Get(Form("hist_HTvsMET_%s",TAG));

	sum.ST.oriA = (TH1D*) fin->Get(Form("histST_%s01",oriTAG)); sum.LPT.oriA = (TH1D*) fin->Get(Form("histLPT_%s01",oriTAG));
	sum.ST.oriB = (TH1D*) fin->Get(Form("histST_%s00",oriTAG)); sum.LPT.oriB = (TH1D*) fin->Get(Form("histLPT_%s00",oriTAG));
	sum.ST.oriC = (TH1D*) fin->Get(Form("histST_%s10",oriTAG)); sum.LPT.oriC = (TH1D*) fin->Get(Form("histLPT_%s10",oriTAG));
	sum.ST.oriD = (TH1D*) fin->Get(Form("histST_%s11",oriTAG)); sum.LPT.oriD = (TH1D*) fin->Get(Form("histLPT_%s11",oriTAG));
	sum.HT.oriA = (TH1D*) fin->Get(Form("histHT_%s01",oriTAG)); sum.MET.oriA = (TH1D*) fin->Get(Form("histMET_%s01",oriTAG));
	sum.HT.oriB = (TH1D*) fin->Get(Form("histHT_%s00",oriTAG)); sum.MET.oriB = (TH1D*) fin->Get(Form("histMET_%s00",oriTAG));
	sum.HT.oriC = (TH1D*) fin->Get(Form("histHT_%s10",oriTAG)); sum.MET.oriC = (TH1D*) fin->Get(Form("histMET_%s10",oriTAG));
	sum.HT.oriD = (TH1D*) fin->Get(Form("histHT_%s11",oriTAG)); sum.MET.oriD = (TH1D*) fin->Get(Form("histMET_%s11",oriTAG));
	sum.NumJet.oriA = (TH1D*) fin->Get(Form("histNumJet_%s01",oriTAG)); sum.NumJet.oriC = (TH1D*) fin->Get(Form("histNumJet_%s10",oriTAG));
	sum.NumJet.oriB = (TH1D*) fin->Get(Form("histNumJet_%s00",oriTAG)); sum.NumJet.oriD = (TH1D*) fin->Get(Form("histNumJet_%s11",oriTAG));
}

void combineEvents(MyProcess &sum, const char* NAME){
	const char *inputFile = Form("%s/%s",PREFIX,NAME);
	TFile *fin = new TFile(inputFile);

	sum.ST.A -> Add( (TH1D*)fin->Get(Form("histST_%s01",TAG)) ); sum.LPT.A -> Add( (TH1D*)fin->Get(Form("histLPT_%s01",TAG)) );
	sum.ST.B -> Add( (TH1D*)fin->Get(Form("histST_%s00",TAG)) ); sum.LPT.B -> Add( (TH1D*)fin->Get(Form("histLPT_%s00",TAG)) );
	sum.ST.C -> Add( (TH1D*)fin->Get(Form("histST_%s10",TAG)) ); sum.LPT.C -> Add( (TH1D*)fin->Get(Form("histLPT_%s10",TAG)) );
	sum.ST.D -> Add( (TH1D*)fin->Get(Form("histST_%s11",TAG)) ); sum.LPT.D -> Add( (TH1D*)fin->Get(Form("histLPT_%s11",TAG)) );
	sum.HT.A -> Add( (TH1D*)fin->Get(Form("histHT_%s01",TAG)) ); sum.MET.A -> Add( (TH1D*)fin->Get(Form("histMET_%s01",TAG)) );
	sum.HT.B -> Add( (TH1D*)fin->Get(Form("histHT_%s00",TAG)) ); sum.MET.B -> Add( (TH1D*)fin->Get(Form("histMET_%s00",TAG)) );
	sum.HT.C -> Add( (TH1D*)fin->Get(Form("histHT_%s10",TAG)) ); sum.MET.C -> Add( (TH1D*)fin->Get(Form("histMET_%s10",TAG)) );
	sum.HT.D -> Add( (TH1D*)fin->Get(Form("histHT_%s11",TAG)) ); sum.MET.D -> Add( (TH1D*)fin->Get(Form("histMET_%s11",TAG)) );
	sum.NumJet.A-> Add( (TH1D*)fin->Get(Form("histNumJet_%s01",TAG)) ); sum.NumJet.C -> Add( (TH1D*)fin->Get(Form("histNumJet_%s10",TAG)) );
	sum.NumJet.B-> Add( (TH1D*)fin->Get(Form("histNumJet_%s00",TAG)) ); sum.NumJet.D -> Add( (TH1D*)fin->Get(Form("histNumJet_%s11",TAG)) );
	sum.NjvsNl  -> Add( (TH2D*)fin->Get(Form("hist_NjvsNl_%s" ,TAG)) );
	sum.HTvsMET -> Add( (TH2D*)fin->Get(Form("hist_HTvsMET_%s",TAG)) );

	sum.ST.oriA -> Add( (TH1D*)fin->Get(Form("histST_%s01",oriTAG)) ); sum.LPT.oriA -> Add( (TH1D*)fin->Get(Form("histLPT_%s01",oriTAG)) );
	sum.ST.oriB -> Add( (TH1D*)fin->Get(Form("histST_%s00",oriTAG)) ); sum.LPT.oriB -> Add( (TH1D*)fin->Get(Form("histLPT_%s00",oriTAG)) );
	sum.ST.oriC -> Add( (TH1D*)fin->Get(Form("histST_%s10",oriTAG)) ); sum.LPT.oriC -> Add( (TH1D*)fin->Get(Form("histLPT_%s10",oriTAG)) );
	sum.ST.oriD -> Add( (TH1D*)fin->Get(Form("histST_%s11",oriTAG)) ); sum.LPT.oriD -> Add( (TH1D*)fin->Get(Form("histLPT_%s11",oriTAG)) );
	sum.HT.oriA -> Add( (TH1D*)fin->Get(Form("histHT_%s01",oriTAG)) ); sum.MET.oriA -> Add( (TH1D*)fin->Get(Form("histMET_%s01",oriTAG)) );
	sum.HT.oriB -> Add( (TH1D*)fin->Get(Form("histHT_%s00",oriTAG)) ); sum.MET.oriB -> Add( (TH1D*)fin->Get(Form("histMET_%s00",oriTAG)) );
	sum.HT.oriC -> Add( (TH1D*)fin->Get(Form("histHT_%s10",oriTAG)) ); sum.MET.oriC -> Add( (TH1D*)fin->Get(Form("histMET_%s10",oriTAG)) );
	sum.HT.oriD -> Add( (TH1D*)fin->Get(Form("histHT_%s11",oriTAG)) ); sum.MET.oriD -> Add( (TH1D*)fin->Get(Form("histMET_%s11",oriTAG)) );
	sum.NumJet.oriA -> Add( (TH1D*)fin->Get(Form("histNumJet_%s01",oriTAG)) ); sum.NumJet.oriC -> Add( (TH1D*)fin->Get(Form("histNumJet_%s10",oriTAG)) );
	sum.NumJet.oriB -> Add( (TH1D*)fin->Get(Form("histNumJet_%s00",oriTAG)) ); sum.NumJet.oriD -> Add( (TH1D*)fin->Get(Form("histNumJet_%s11",oriTAG)) );
}

void writeEvents(MyProcess sum){
	sum.ST.oriA -> Write(); sum.LPT.oriA -> Write(); sum.HT.oriA -> Write(); sum.MET.oriA -> Write();
	sum.ST.oriB -> Write(); sum.LPT.oriB -> Write(); sum.HT.oriB -> Write(); sum.MET.oriB -> Write();
	sum.ST.oriC -> Write(); sum.LPT.oriC -> Write(); sum.HT.oriC -> Write(); sum.MET.oriC -> Write();
	sum.ST.oriD -> Write(); sum.LPT.oriD -> Write(); sum.HT.oriD -> Write(); sum.MET.oriD -> Write();
	sum.NumJet.oriA -> Write(); sum.NumJet.oriC -> Write();
	sum.NumJet.oriB -> Write(); sum.NumJet.oriD -> Write();

	sum.ST.A -> Write(); sum.LPT.A -> Write(); sum.HT.A -> Write(); sum.MET.A -> Write();
	sum.ST.B -> Write(); sum.LPT.B -> Write(); sum.HT.B -> Write(); sum.MET.B -> Write();
	sum.ST.C -> Write(); sum.LPT.C -> Write(); sum.HT.C -> Write(); sum.MET.C -> Write();
	sum.ST.D -> Write(); sum.LPT.D -> Write(); sum.HT.D -> Write(); sum.MET.D -> Write();
	sum.NumJet.A-> Write(); sum.NumJet.C -> Write(); sum.NjvsNl  -> Write();
	sum.NumJet.B-> Write(); sum.NumJet.D -> Write(); sum.HTvsMET -> Write();
}

int main(int argc, char* argv[]){
	const char *savingPath = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/skimmed"; 
	const char *savingName = "simulation_delphes_ppvv_skimmed.root";

	MyProcess sum;
	importEvents(sum,"ppvv_run_1M_0.root");
	for(int i=1; i<10; i++)
		combineEvents(sum,Form("ppvv_run_1M_%d.root",i));

	TFile *fout = new TFile(Form("%s/%s",savingPath,savingName),"recreate");
	writeEvents(sum);

	std::cout<<"The info of 10 separated ppvv rootfile are Chained!"<<std::endl;

	fout->Close();
	return 1;
}
