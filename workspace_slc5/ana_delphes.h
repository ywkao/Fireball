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
#include <TNtuple.h>
#include <TMath.h>

#include "xsec.h"
using namespace std;

const int num_bg=8;
const int NUM=9;
const int NUM_hist=14;
const int NUM_CUT=9;
const double L=100.;
const double pb2fb=1000.;
const double L_MASS=60.;//mass lower bound for bonson candidates
const double U_MASS=110.;//mass upper bound for bonson candidates
char *histName[NUM_hist]={"Num_lep","Num_jet","PT_lep","PT_jet","LPT","HT","MET","ST","Num_boson","M_boson","PT_chJet","Eta_chJet","Phi_chJet","Gen_HT"};

double X;
double Err_X;
int PRESELECTION;
int TAIL;//record the input value

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
class Particle_Gen : public Particle {
public:
	int status, PID;
};

class MyEvent {
public:
	std::vector<Particle> Electrons;
	std::vector<Particle> Muons;
	std::vector<Particle> METs;
	std::vector<Particle_Jet> Jets;
	std::vector<Particle_Boson> H_Bosons;
	std::vector<Particle_Gen> GenParticles;
	double GenHT;
	bool PreSelectionPass;
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
//### Parton-level PreSelection
//#############################################
const double ihtmin=500;//Preselection ihtmin=500GeV
bool PreSelection(int GenNumLep, double GenHT){
	if( GenNumLep==0 || GenHT<ihtmin ) return false;
	else return true;
}


//#############################################
//### The Functions to Load Data
//#############################################
void ChainingEvents(string PROCESS, TChain *&chain){
	if(PROCESS=="pptt"){
	//chain->Add("/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/pptt_Go02/Events/run_ht_min_2500/tag_1_delphes_events.root");
		if(PRESELECTION==0 || PRESELECTION==1){//apply on the events w/o preselection
			//chain->Add("/raid2/w/ykao/simulation/pptt_31/Events/run_100k_31_0/tag_1_delphes_events.root");
			for(int i=1; i<30; i++){
				if( i==1 || i==10 || i==20) continue;
				else if(i<10)
					for(int j=0; j<10; j++)	chain->Add(Form("/raid1/w/ykao/simulation/pptt_0%d/Events/run_100k_0%d_%d/tag_1_delphes_events.root",i,i,j));
				else
					for(int j=0; j<10; j++)	chain->Add(Form("/raid1/w/ykao/simulation/pptt_%d/Events/run_100k_%d_%d/tag_1_delphes_events.root",i,i,j));
			}
		} else{
			//chain->Add("/raid2/w/ykao/simulation/pptt_51/Events/run_20k_51_0/tag_1_delphes_events.root");
			chain->Add("/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/pptt_52/Events/run_100k_52_0/tag_1_delphes_events.root");
			chain->Add("/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/pptt_52/Events/run_100k_52_1/tag_1_delphes_events.root");
			for(int i=41; i<50; i++){
				for(int j=0; j<10; j++)	chain->Add(Form("/raid2/w/ykao/simulation/pptt_%d/Events/run_10k_%d_%d/tag_1_delphes_events.root",i,i,j));
			}
			for(int i=51; i<60; i++){
				if( i==52 ) continue;
				for(int j=0; j<10; j++){
					if( i==55 && j==9 ) continue;
					chain->Add(Form("/raid2/w/ykao/simulation/pptt_%d/Events/run_20k_%d_%d/tag_1_delphes_events.root",i,i,j));
				}
			}
		}
	}
	if(PROCESS=="ppvv"){
		if(PRESELECTION==0 || PRESELECTION==1){//apply on the events w/o preselection
			for(int j=0; j<10; j++)	chain->Add(Form("/raid1/w/ykao/simulation/ppvv_jet_matching/Events/run_1M_%d/tag_1_delphes_events.root",j));
			for(int i=1; i<3; i++){
				for(int j=0; j<10; j++)	chain->Add(Form("/raid1/w/ykao/simulation/ppvv_0%d/Events/run_100k_0%d_%d/tag_1_delphes_events.root",i,i,j));
			}
		} else{
			for(int i=11; i<13; i++){
				for(int j=0; j<10; j++)	chain->Add(Form("/raid2/w/ykao/simulation/ppvv_%d/Events/run_10k_%d_%d/tag_1_delphes_events.root",i,i,j));
			}
		}
	}
	if(PROCESS=="ppvtt"){
		chain->Add("/raid1/w/ykao/01source/simulation_delphes_ppvtt_jet_matching.root");
		for(int j=0; j<4; j++)	chain->Add(Form("/raid1/w/ykao/simulation/ppvtt_01/Events/run_100k_01_%d/tag_1_delphes_events.root",j));
	}
	if(PROCESS=="ppvvv"){
		chain->Add("/raid1/w/ykao/01source/simulation_delphes_ppvvv_jet_matching.root");
		for(int j=0; j<4; j++)	chain->Add(Form("/raid1/w/ykao/simulation/ppvvv_01/Events/run_100k_01_%d/tag_1_delphes_events.root",j));
	}
}
std::vector<MyEvent> vec;
//void importEvents(std::vector<MyEvent> &vec, TString filename){
void importEvents(std::vector<MyEvent> &vec, ExRootTreeReader *TreeReader){
	TClonesArray *branchElectron = TreeReader->UseBranch("Electron");
	TClonesArray *branchMuon = TreeReader->UseBranch("Muon");
	TClonesArray *branchJet = TreeReader->UseBranch("Jet");
	TClonesArray *branchMissingET = TreeReader->UseBranch("MissingET");
	TClonesArray *branchGenParticle = TreeReader->UseBranch("Particle");
	//TClonesArray *branchGenJet = TreeReader->UseBranch("GenJet");
	//TClonesArray *branchPhoton = TreeReader->UseBranch("Photon");

	for(int entry=0; entry<TreeReader->GetEntries(); entry++){
		if((entry+1)%10000==0 || entry+1==TreeReader->GetEntries()) printf("Processing: %d/%d\n",entry+1,TreeReader->GetEntries());
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

		for(int i=0; i<branchGenParticle->GetEntries(); i++){
			if(i==0) event.GenHT = 0;
			GenParticle *genparticle = (GenParticle*)branchGenParticle->At(i);
			if(genparticle->Status!=3) continue;
			Particle_Gen particle;
			particle.status = genparticle->Status;
			particle.PID    = genparticle->PID;
			particle.pt     = genparticle->PT;
			
			int status = genparticle->Status;
			int PID = genparticle->PID;
			//if(i==0) cout<<"##### "<<"Entries = "<<entry+1<<" #####"<<endl;
			//printf("Status = %d, PID = %d, PT=%6.2f\n", status, PID, particle.pt);

			if((abs(PID)==1 || abs(PID)==2 || abs(PID)==3 || abs(PID)==4 || abs(PID)==5 || abs(PID)==21)) event.GenHT += genparticle->PT;

			event.GenParticles.push_back(particle);
		}
		
		event.PreSelectionPass = PreSelection(event.Electrons.size()+event.Muons.size(), event.GenHT);

		vec.push_back(event);

	}//end of entry

	std::cout << vec.size() << " events stored" << std::endl;
}


//#############################################
//### Declare Hist/Ntuple to Store Wanted INFO
//#############################################
TH1D *hist_Ori[NUM_hist];//before cuts 
TH1D *hist_Cut[NUM_hist];//after cuts 
TH1D *hist_Nm1[NUM_hist];//N-1 plots 
TNtuple *ntuple_Ori[NUM_hist];
TNtuple *ntuple_Nm1[NUM_hist];
TNtuple *ntuple_Cut[NUM_hist];
TNtuple *ntuple_Ori_HL;
TNtuple *ntuple_Cut_HL;

int    n_bin[NUM_hist]={10,30,120,120,120,280,120,320,12,U_MASS-L_MASS,80,16,16,280};
double l_bin[NUM_hist]={ 0, 0, 0, 0, 0, 0, 0, 0, 0,		  L_MASS, 0,-8,-8, 0};
double h_bin[NUM_hist]={10,30,3000,3000,3000,7000,3000,8000,12,U_MASS,2000,8,8,7000};

//#############################################
//### Miscellaneous Functions
//#############################################
FILE *LOGFile;
const char *message;
void INFOLOG(const char *message){
	LOGFile = fopen("INFOLOG","a");
	fputs(Form("%s\n",message),LOGFile);
	std::cout<<message<<std::endl;
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
