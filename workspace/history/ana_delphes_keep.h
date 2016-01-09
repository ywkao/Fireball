#include <stdio.h>//FILE
#include <iostream>
#include <vector>
#include <math.h>
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include "xsec.h"
using namespace std;

const int num_bg=8;
const int NUM=9;
const int NUM_hist=13;
const int NUM_CUT=9;
const double L=100.;
const double pb2fb=1000.;
const double L_MASS=60.;//mass lower bound for bonson candidates
const double U_MASS=110.;//mass upper bound for bonson candidates
//double X_sig[2] = {0.051, 0.000177}; //1TeV, 2TeV
//double Err_X_sig[2] = {0.001, 4e-06}; //1TeV, 2TeV
//double X_bg[8] = {815.96, 192.4, 1.77, 0.621, 0.0219, 0.0164, 0.00405, 0.000249}; //pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
//double Err_X_bg[8] = {45.51, 4.1, 0.10, 0.073, 0.0013, 9e-04, 9e-05, 3.4e-05}; //pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt

char *histName[13]={"Num_lep","Num_jet","PT_lep","PT_jet","tot_Lep_PT","HT","MET","ST","Num_boson","M_boson","PT_chosenJet","Eta_chosenJet","Phi_chosenJet"};
char *processes[NUM] = {"fireball", "pptt", "ppvv", "ppvtt", "ppvvv", "ppvvtt", "pptttt", "ppvvvv", "ppvvvtt"};
char *units[NUM_hist] = {"# of lep", "# of jet", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "# of boson", "GeV", "GeV", "eta", "rad."};

int SIG;//record the input value
char *savingPath;
//char *suffix;
double X[NUM];
double Err_X[NUM];

//#############################################
//### Set The Format to Store Particle Info
//#############################################
class Particle {
public:
	double pt, eta, phi;
};

class Particle_Jet : public Particle {
public:
	int BTag;
};

class Particle_Boson : public Particle {
public:
	std::pair<int,int> jetIds;
	double mass;
};

class MyEvent {
public:
	std::vector<Particle> Electrons;
	std::vector<Particle> Muons;
	std::vector<Particle> METs;
	std::vector<Particle_Jet> Jets;
	std::vector<Particle_Boson> H_Bosons;
};

//#############################################
//### Deal With Hadronic Bonsons Counting
//#############################################
//std::vector<std::pair<std::pair<int,int>,double>> bosonCandidates;
void Sorting(std::vector< std::pair<std::pair<int,int>,double> > &bosonCandidates){
	double MASS = 90.; //take W-boson mass as the sorting norm
	for(int i=0; i<bosonCandidates.size(); i++){
		for(int j=0; j<i; j++){
			if(fabs(bosonCandidates.at(j).second-MASS) > fabs(bosonCandidates.at(i).second-MASS)){
				std::pair<int,int> jetIds = bosonCandidates.at(i).first;
				double mass = bosonCandidates.at(i).second;
				bosonCandidates.at(i).first  = bosonCandidates.at(j).first;
				bosonCandidates.at(i).second = bosonCandidates.at(j).second;
				bosonCandidates.at(j).first  = jetIds;
				bosonCandidates.at(j).second = mass;
			}
		}
	}
}

void Removing(std::vector< std::pair<std::pair<int,int>,double> > &bosonCandidates){
	std::pair<int,int> Ids = bosonCandidates.at(0).first;
	for(int i=0; i<bosonCandidates.size(); i++){
		if(bosonCandidates.at(i).first.first==Ids.first) bosonCandidates.at(i).second = 0.;
		if(bosonCandidates.at(i).first.first==Ids.second) bosonCandidates.at(i).second = 0.;
		if(bosonCandidates.at(i).first.second==Ids.first) bosonCandidates.at(i).second = 0.;
		if(bosonCandidates.at(i).first.second==Ids.second) bosonCandidates.at(i).second = 0.;
	}
}

//#############################################
//### The Function to Load Data
//#############################################
std::vector<MyEvent> vec_sig, vec_bg[num_bg];
void importEvents(std::vector<MyEvent> &vec, TString filename){
	std::cout<<"Loading: "<<filename<<std::endl;

	TChain *chain = new TChain("Delphes");
	chain->Add(filename);

	ExRootTreeReader *TreeReader = new ExRootTreeReader(chain);
	TClonesArray *branchElectron = TreeReader->UseBranch("Electron");
	TClonesArray *branchMuon = TreeReader->UseBranch("Muon");
	TClonesArray *branchJet = TreeReader->UseBranch("Jet");
	TClonesArray *branchMissingET = TreeReader->UseBranch("MissingET");
	//TClonesArray *branchParticle = TreeReader->UseBranch("Particle");
	//TClonesArray *branchGenJet = TreeReader->UseBranch("GenJet");
	//TClonesArray *branchPhoton = TreeReader->UseBranch("Photon");

	for(int entry=0; entry<TreeReader->GetEntries(); entry++){
		TreeReader->ReadEntry(entry);
		MyEvent event;

		for(int i=0; i<branchElectron->GetEntries(); i++){
			Particle particle;
			Electron *electron = (Electron*)branchElectron->At(i);
			particle.pt = electron->PT;
			particle.eta = electron->Eta;
			particle.phi = electron->Phi;
			event.Electrons.push_back(particle);
		}

		for(int i=0; i<branchMuon->GetEntries(); i++){
			Particle particle;
			Muon *muon = (Muon*)branchMuon->At(i);
			particle.pt = muon->PT;
			particle.eta = muon->Eta;
			particle.phi = muon->Phi;
			event.Muons.push_back(particle);
		}

		for(int i=0; i<branchJet->GetEntries(); i++){
			Particle_Jet particle;
			Jet *jet = (Jet*)branchJet->At(i);
			particle.pt = jet->PT;
			particle.eta = jet->Eta;
			particle.phi = jet->Phi;
			particle.BTag = jet->BTag;
			event.Jets.push_back(particle);
		}

		for(int i=0; i<branchMissingET->GetEntries(); i++){
			Particle particle;
			MissingET *met = (MissingET*)branchMissingET->At(i);
			particle.pt = met->MET;
			particle.eta = met->Eta;
			particle.phi = met->Phi;
			event.METs.push_back(particle);
		}

		std::vector< std::pair<std::pair<int,int>,double> > bosonCandidates;
		for(int i=0; i<branchJet->GetEntries(); i++){
			for(int j=i+1; j<branchJet->GetEntries(); j++){
				Jet *jet1 = (Jet*) branchJet->At(i);
				Jet *jet2 = (Jet*) branchJet->At(i);
				double mass = (jet1->P4()+jet2->P4()).Mag();
				if(mass>L_MASS && mass<U_MASS)//m_w=80, m_z=91.2
					bosonCandidates.push_back(std::make_pair(std::make_pair(i,j),mass));
			}
		}

		Sorting(bosonCandidates);
		for(int i=0; i<bosonCandidates.size(); i++){
			if(bosonCandidates.at(0).second == 0.)	break;// the candidate jet has been used
			Particle_Boson particle;
			particle.jetIds = bosonCandidates.at(0).first;
			particle.mass = bosonCandidates.at(0).second;
			event.H_Bosons.push_back(particle);
			Removing(bosonCandidates);
			Sorting(bosonCandidates);
		}

		vec.push_back(event);

	}//end of entry

	std::cout << vec.size() << " events stored" << std::endl;
}

//#############################################
//### Initialization & Selection Cuts
//#############################################
TH1D *hist_Ori[NUM_hist][NUM];//before cuts 
TH1D *hist_Cut[NUM_hist][NUM];//after cuts 
TH1D *hist_Nm1[NUM_hist][NUM];//N-1 plots 
//THStack *hist_Stack_Ori[NUM_CUT];//before cuts 
//THStack *hist_Stack_Cut[NUM_CUT];//after cuts 
//THStack *hist_Stack_Copy[NUM_CUT];//before cuts 
TLine *CutLine[NUM_CUT];//the cut lines
int CutLineColor = kBlack;

int    n_bin[2][NUM_hist]={{10,20,20,20,20,20,20,20,12,U_MASS-L_MASS,60,16,16},
						   {10,20,20,20,20,20,20,20,12,U_MASS-L_MASS,100,16,16}};
double l_bin[2][NUM_hist]={{0,0,0,0,0,0,0,0,0,L_MASS,0,-8,-8},
						   {0,0,0,0,0,0,0,0,0,L_MASS,0,-8,-8}};
double h_bin[2][NUM_hist]={{10,40, 800,1400, 600,5000,1200,7500,12,U_MASS,1200,8,8},
						   {10,40,1000,2000,1000,7000,1500,8000,12,U_MASS,1600,8,8}};
//int    n_bin_ori[2][NUM_hist]={{10,40,50,50,50,50,50,50,12,U_MASS-L_MASS,60,16,16},
//							   {10,40,50,50,50,35,50,40,12,U_MASS-L_MASS,100,16,16}};
//int    n_bin_cut[2][NUM_hist]={{ 8,40,50,50,50,50,50,50,12,U_MASS-L_MASS,60,16,16},
//							   {12,40,50,50,50,50,50,50,12,U_MASS-L_MASS,100,16,16}};
//int    n_bin_cut[2][NUM_hist]={{ 8,20,20,20,20,20,20,20,12,U_MASS-L_MASS,60,16,16},
//							   {12,20,20,20,20,20,20,20,12,U_MASS-L_MASS,100,16,16}};
//double l_bin_cut[2][NUM_hist]={{0,0,0,0,0,0,0,0,0,L_MASS,0,-8,-8},
//							   {0,0,0,0,0,0,0,0,0,L_MASS,0,-8,-8}};
//double h_bin_cut[2][NUM_hist]={{ 8,40, 800,1400, 600,5000,1200,7500,12,U_MASS,1200,8,8},
//							   {12,40,1000,2000,1000,5500,1500,8000,12,U_MASS,1600,8,8}};
//double scaleBosonMass, scaleChosenJetEta, scaleChosenJetPhi, scaleLepSPT;//hist->SetMaximum(scale)
//if(SIG==1)	{scaleBosonMass = 20000, scaleChosenJetEta = 20000, scaleChosenJetPhi = 200000, scaleLepSPT = 1000;}
//if(SIG==2)  {scaleBosonMass = 500, scaleChosenJetEta = 10000, scaleChosenJetPhi = 5000, scaleLepSPT = 1000;}

//#############################################
//### Miscellaneous Functions
//#############################################
FILE *LOGFile;
const char *message;
void INFOLOG(const char *message){
	LOGFile = fopen("INFOLOG","a");
	fputs(Form("%s\n",message),LOGFile);
	fclose(LOGFile);
}
void ListSelectionCuts(std::vector<double> factors){
	int n = 0;
    for( std::vector<Double_t>::iterator it = factors.begin(); it<factors.end(); it++ ){
		message = Form("%s: %5.0f", histName[n], (*it)); INFOLOG(message);
        std::cout << setw(11) <<histName[n] << ": " << (*it) << std::endl;
        n++;
    }
}
struct CrossSection {
	//std::pair<double,double> pptt = std::make_pair(815.96,45.51); WHY DOES IT NOT WORK??
	double GetCrossSection(const char* Name){
		if(Name=="fireball_bp_1TeV") 	return fireball_bp_1TeV;
		if(Name=="fireball_bp_2TeV") 	return fireball_bp_2TeV;
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
		if(Name=="fireball_bp_1TeV") 	return Err_fireball_bp_1TeV;
		if(Name=="fireball_bp_2TeV") 	return Err_fireball_bp_2TeV;
		if(Name=="pptt") 	return Err_pptt;
		if(Name=="ppvv")    return Err_ppvv;
		if(Name=="ppttv")   return Err_ppttv;
		if(Name=="ppvvv")   return Err_ppvvv;
		if(Name=="ppttvv")  return Err_ppttvv;
		if(Name=="pptttt")  return Err_pptttt;
		if(Name=="ppvvvv")  return Err_ppvvvv;
		if(Name=="ppttvvv") return Err_ppttvvv;
	}
}; CrossSection Xsec;


/*
void SetStyle(){
	gStyle->SetOptStat(0);
}
void DrawSetting(TH1D *&hist, TLegend *&leg, char *text, char* title, char *xtitle, char *ytitle, int style, int color){
	hist->SetTitle(title);
	hist->SetXTitle(xtitle);
	hist->SetYTitle(ytitle);
	hist->GetYaxis()->SetTitleOffset(1.2);
	hist->SetFillStyle(style);
	hist->SetFillColor(color);
	hist->SetLineColor(1);
	hist->SetLineWidth(0);
	leg ->AddEntry(hist,text,"f");
	leg ->SetTextSize(0.042);
}
void DrawSettingGen(TNtuple *&ntuple, TLegend *&leg, char *text, Double_t markerSize, int markerStyle, int lineStyle, int lineWidth, int color){
	ntuple->SetMarkerSize(markerSize);
	ntuple->SetMarkerColor(color);
	ntuple->SetMarkerStyle(markerStyle);
	ntuple->SetLineColor(color);
	ntuple->SetLineStyle(lineStyle);
	ntuple->SetLineWidth(lineWidth);
	leg->AddEntry(ntuple,text,"lp");//lepf
}
*/

//void PrintSelectionCuts(std::vector<double> factors){
//	std::cout<<"factors are:";
//	for( std::vector<Double_t>::iterator it = factors.begin(); it!=factors.end(); it++)
//		std::cout<<' '<<*it; 
//	std::cout<<'\n';
//}

