#include <stdio.h>//FILE
#include <iostream>
#include <iomanip> //for using setw(n)
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

char *histName[13]={"Num_lep","Num_jet","PT_lep","PT_jet","tot_Lep_PT","HT","MET","ST","Num_boson","M_boson","PT_chosenJet","Eta_chosenJet","Phi_chosenJet"};
char *processes[NUM] = {"fireball", "pptt", "ppvv", "ppvtt", "ppvvv", "ppvvtt", "pptttt", "ppvvvv", "ppvvvtt"};
char *units[NUM_hist] = {"# of lep", "# of jet", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "# of boson", "GeV", "GeV", "eta", "rad."};

int TAIL;//record the input value
//char *savingPath;
//char *suffix;
//double X[NUM];
//double Err_X[NUM];
double X;
double Err_X;

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
//std::vector<MyEvent> vec_sig, vec_bg[num_bg];
std::vector<MyEvent> vec;
//void importEvents(std::vector<MyEvent> &vec, TString filename){
void importEvents(std::vector<MyEvent> &vec, ExRootTreeReader *TreeReader){
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
TH1D *hist_Ori[NUM_hist];//before cuts 
TH1D *hist_Cut[NUM_hist];//after cuts 
TH1D *hist_Nm1[NUM_hist];//N-1 plots 
//TLine *CutLine[NUM_CUT];//the cut lines
//int CutLineColor = kBlack;

int    n_bin[NUM_hist]={10,30,20,20,20,20,20,20,12,U_MASS-L_MASS,100,16,16};
double l_bin[NUM_hist]={0,0,0,0,0,0,0,0,0,L_MASS,0,-8,-8};
double h_bin[NUM_hist]={10,30, 1000,2000, 1000,7000,1500,8000,12,U_MASS,1600,8,8};

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
	double GetCrossSection(const char* Name){
		if((string)Name=="pptt") 	return pptt;
		if((string)Name=="ppvv")    return ppvv;
		if((string)Name=="ppvtt")   return ppvtt;
		if((string)Name=="ppvvv")   return ppvvv;
		if((string)Name=="ppvvtt")  return ppvvtt;
		if((string)Name=="pptttt")  return pptttt;
		if((string)Name=="ppvvvv")  return ppvvvv;
		if((string)Name=="ppvvvtt") return ppvvvtt;
	}
	double GetCrossSectionErr(const char* Name){
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
