#include <signal.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>

#include <iomanip> //for using setw(n)
#include <fstream>
#include <stdlib.h>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TPad.h>

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"

using namespace std;
using namespace TMVA;

const int MaxBin = 20;
const int num_bg = 8;
const int NUM = num_bg+1;
const int NUM_CUT = 9;
const double L=100.;
const double pb2fb=1000.;
const double L_MASS=60.;//mass lower bound for bonson candidates
const double U_MASS=110.;//mass upper bound for bonson candidates

double X[NUM], Err_X[NUM];//to be determined by SIG (the Xsections of all considered processes)
double X_sig[2] = {0.051, 0.000177}; //1TeV, 2TeV
double Err_X_sig[2] = {0.001, 4e-06}; //1TeV, 2TeV
double X_bg[8] = {815.96, 192.4, 1.77, 0.621, 0.0219, 0.0164, 0.00405, 0.000249}; //pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
double Err_X_bg[8] = {45.51, 4.1, 0.10, 0.073, 0.0013, 9e-04, 9e-05, 3.4e-05}; //pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
char *processes[NUM] = {"fireball", "pptt", "ppvv", "ppvtt", "ppvvv", "ppvvtt", "pptttt", "ppvvvv", "ppvvvtt"};
char *cutName[NUM_CUT]={"Num_lep","Num_jet","PT_lep","PT_jet","tot_Lep_PT","HT","MET","ST","Num_boson"};

int SIG;//record the input value
char *savingPath;
char *suffix;

//#############################################
//### Set The Format to Store Particle Info
//#############################################
class Particle {
public:
	//static int status;
	double pt, eta, phi;
};
//int Particle::status = 0;

class Particle_Jet : public Particle {
public:
	int BTag;
};

class Particle_Boson : public Particle {
public:
	std::pair<int,int> jetIds;
	double mass;
};

class status{
public:
	static int flag_NUM_lep;
	static int flag_NUM_jet;
	static int flag_SPT_lep;
	static int flag_SPT_jet;
	static int flag_MET;
};

int flag_NUM_lep = 1;
int flag_NUM_jet = 1;
int flag_SPT_lep = 1;
int flag_SPT_jet = 1;
int flag_MET = 1;

class MyEvent {
public:
	status FLAG;
	std::vector<Particle> Electrons;
	std::vector<Particle> Muons;
	std::vector<Particle> METs;
	std::vector<Particle_Jet> Jets;
	std::vector<Particle_Boson> H_Bosons;
};

//#############################################
//### Deal With Hadronic Bonsons Counting
//#############################################
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
//### Miscellaneous Functions
//#############################################
void ListSelectionCuts(std::vector<double> factors){
	int n = 0;
    for( std::vector<Double_t>::iterator it = factors.begin(); it<factors.end(); it++ ){
        std::cout << setw(11) <<cutName[n] << ": " << (*it) << std::endl;
        n++;
    }
}

double GetST(std::vector<MyEvent>::iterator it){
	double SPT_lep=0, SPT_jet=0, TOT_met=0;
	for(int i=0; i < it->Electrons.size(); i++)	SPT_lep += it->Electrons[i].pt;
	for(int i=0; i < it->Muons.size(); i++)		SPT_lep += it->Muons[i].pt;
	for(int i=0; i < it->Jets.size(); i++)		SPT_jet += it->Jets[i].pt;
	for(int i=0; i < it->METs.size(); i++) 		TOT_met += it->METs[i].pt;
	return SPT_lep + SPT_jet + TOT_met;
}
double GetHT(std::vector<MyEvent>::iterator it){
	double SPT_jet=0;
	for(int i=0; i < it->Jets.size(); i++)		SPT_jet += it->Jets[i].pt;
	return SPT_jet;
}
double GetMET(std::vector<MyEvent>::iterator it){
	double TOT_met=0;
	for(int i=0; i < it->METs.size(); i++) 		TOT_met += it->METs[i].pt;
	return TOT_met;
}
double GetSPTLep(std::vector<MyEvent>::iterator it){
	double SPT_lep=0;
	for(int i=0; i < it->Electrons.size(); i++)	SPT_lep += it->Electrons[i].pt;
	for(int i=0; i < it->Muons.size(); i++)		SPT_lep += it->Muons[i].pt;
	return SPT_lep;
}

void RecordIndividualSpectrums(std::vector<MyEvent>::iterator it, double weighting, TH1D *&hist_ST, TH1D *&hist_HT, TH1D *&hist_MET, TH1D *&hist_LPT){
	hist_ST->Fill(GetST(it),weighting);
	hist_HT->Fill(GetHT(it),weighting);
	hist_MET->Fill(GetMET(it),weighting);
	hist_LPT->Fill(GetSPTLep(it),weighting);
}
int DecideLocation(int k){//k=2i+j
	if(k==0) return 3; //B(0,0)
	if(k==1) return 1; //A(0,1)
	if(k==2) return 4; //C(1,0)
	if(k==3) return 2; //D(1,1)
}
double GetTotContent(TH1D *hist){
	double sum = 0.;
	for(int bin=1; bin<MaxBin+1; bin++) sum += hist->GetBinContent(bin);
	return sum;
}

//#############################################
//### THE MAIN FUNCTION
//#############################################
int main(int argc, char* argv[]){
	std::cout << "The chosen bp mass is " << argv[1] << " TeV" << std::endl;
	SIG = atoi(argv[1]);
	suffix = argv[2];
	savingPath = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/ToBeProcessed"; 

	//###Load processes: pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
	if(SIG==1) importEvents(vec_sig,"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_fireball_bp_1TeV.root");
	if(SIG==2) importEvents(vec_sig,"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_fireball_bp_2TeV.root");
	importEvents(vec_bg[0],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_pptt_jet_matching.root.backup");
	importEvents(vec_bg[1],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvv_jet_matching.root.backup");
	importEvents(vec_bg[2],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvtt_jet_matching.root");
	importEvents(vec_bg[3],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvv_jet_matching.root");
	importEvents(vec_bg[4],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvtt_jet_matching.root");
	importEvents(vec_bg[5],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_pptttt.root");
	importEvents(vec_bg[6],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvvv_jet_matching.root");
	importEvents(vec_bg[7],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvvtt.root");

	//###set cross section list
	for(int i=0; i<NUM; i++){
		if(i==0) {X[i] = X_sig[SIG-1]; Err_X[i] = Err_X_sig[SIG-1];}
		else {X[i] = X_bg[i-1]; Err_X[i] = Err_X_bg[i-1];}
	}

	//###set selection cuts
	std::vector<double> factors;//to store the optimized selection cuts
	double list01[9]={2,6,30,40,20,1200,50,1400,0};
	double list02[9]={2,6,30,40,20,1200,50,1400,0};
	if(SIG==1) factors.insert(factors.begin(),list01,list01+9);
	if(SIG==2) factors.insert(factors.begin(),list02,list02+9);
    std::cout << "\nChosen factors:" << std::endl;
	ListSelectionCuts(factors);//print out the best seclectionCuts!
	double CUT_Num_lep   = factors.at(0);
	double CUT_Num_jet   = factors.at(1);
	double CUT_PT_lep    = factors.at(2);
	double CUT_PT_jet    = factors.at(3);
	double CUT_SPT_lep   = factors.at(4);
	double CUT_SPT_jet   = factors.at(5);
	double CUT_MET     	 = factors.at(6);
	double CUT_ST 		 = factors.at(7);
	double CUT_Num_boson = factors.at(8);

	//TFile *fout = new TFile(Form("%s/result_%dTeV.root",savingPath,SIG),"recreate");
	TFile *fout = new TFile(Form("%s/result.root",savingPath),"recreate");
	TCanvas *can = new TCanvas("can","",800,600);
	TH1D *hist_ST[2][2], *hist_HT[2][2], *hist_MET[2][2], *hist_LPT[2][2];
	TH2D *hist_NlvsNj  = new TH2D("hist_NlvsNj ",";# of lep;# of jet",7,0,7,25,0,25);
	TH2D *hist_HTvsMET = new TH2D("hist_HTvsMET",";MET (60 GeV); HT (250 GeV)",20,0,1200,20,0,5000);
	gStyle -> SetStatX(0.9);
	gStyle -> SetStatY(0.9);

	char *Label[4] = {"B","A","C","D"};//A(0,1) B(0,0) C(1,0) D(1,1)
	for(int i=0; i<2; i++){
		for(int j=0; j<2; j++){
			hist_ST[i][j]  = new TH1D(Form("hist_ST_%d%d",i,j) ,Form("%s;GeV;Entries / 50 GeV",Label[2*i+j]),20,0,2000); hist_ST[i][j] ->Sumw2(); 
			hist_HT[i][j]  = new TH1D(Form("hist_HT_%d%d",i,j) ,Form("%s;GeV;Entries / 50 GeV",Label[2*i+j]),20,0,2000); hist_HT[i][j] ->Sumw2(); 
			hist_MET[i][j] = new TH1D(Form("hist_MET_%d%d",i,j),Form("%s;GeV;Entries / 50 GeV",Label[2*i+j]),20,0, 800); hist_MET[i][j]->Sumw2(); 
			hist_LPT[i][j] = new TH1D(Form("hist_LPT_%d%d",i,j),Form("%s;GeV;Entries / 50 GeV",Label[2*i+j]),20,0, 800); hist_LPT[i][j]->Sumw2(); 
		}
	}
	
	//###ABCD Nethod
	//int Num_lep = 0, Num_jet = 0, A[num_bg][2][2]={0}, totA[2][2]={0};//A: event entries based on Num lep/jet
	int Num_lep = 0, Num_jet = 0;//A: event entries based on Num lep/jet
	double HT = 0., MET = 0.;
	double totA[2][2]={0.} 		   ,Err_totA[2][2]={0.};
	double totYield[2][2]={0.} 	   ,Err_totYield[2][2]={0.};
	double A[num_bg][2][2]={0.}	   ,Err_A[num_bg][2][2]={0.};
	double Eff[num_bg][2][2]={0.}  ,Err_Eff[num_bg][2][2]={0.};
	double Yield[num_bg][2][2]={0.},Err_Yield[num_bg][2][2]={0.};

	for(int N=0; N<num_bg; N++){
		//std::vector<MyEvent> vec; if(N==0) vec = vec_sig; else vec = vec_bg[N-1];
		std::vector<MyEvent> vec; vec = vec_bg[N];
		double weighting = L*pb2fb*X[N+1]/(double)vec.size();
		for( std::vector<MyEvent>::iterator it = vec.begin(); it!=vec.end(); it++){
			/*##### Num_jet vs. Num_lep #####*/
			Num_lep = it->Electrons.size()+it->Muons.size();
			Num_jet = it->Jets.size();
			hist_NlvsNj -> Fill(Num_lep,Num_jet,weighting);
			if(Num_lep<CUT_Num_lep && Num_jet<CUT_Num_jet)
				{A[N][0][0]+=1;RecordIndividualSpectrums(it,weighting,hist_ST[0][0],hist_HT[0][0],hist_MET[0][0],hist_LPT[0][0]);}
			if(!(Num_lep<CUT_Num_lep) && Num_jet<CUT_Num_jet)
				{A[N][1][0]+=1;RecordIndividualSpectrums(it,weighting,hist_ST[1][0],hist_HT[1][0],hist_MET[1][0],hist_LPT[1][0]);}
			if(Num_lep<CUT_Num_lep && !(Num_jet<CUT_Num_jet))
				{A[N][0][1]+=1;RecordIndividualSpectrums(it,weighting,hist_ST[0][1],hist_HT[0][1],hist_MET[0][1],hist_LPT[0][1]);}
			if(!(Num_lep<CUT_Num_lep) && !(Num_jet<CUT_Num_jet))
				{A[N][1][1]+=1;RecordIndividualSpectrums(it,weighting,hist_ST[1][1],hist_HT[1][1],hist_MET[1][1],hist_LPT[1][1]);}

			/*##### HT vs. MET #####
			HT  = GetHT(it);
			MET = GetMET(it);
			hist_HTvsMET -> Fill(MET,HT);
			if(MET<CUT_MET && HT<CUT_SPT_jet && !(HT<1100))
				{A[N][0][0]+=1;RecordIndividualSpectrums(it,weighting,hist_ST[0][0],hist_HT[0][0],hist_MET[0][0],hist_LPT[0][0]);}
			if(!(MET<CUT_MET) && HT<CUT_SPT_jet && !(HT<1100))
				{A[N][1][0]+=1;RecordIndividualSpectrums(it,weighting,hist_ST[1][0],hist_HT[1][0],hist_MET[1][0],hist_LPT[1][0]);}
			if(MET<CUT_MET && !(HT<CUT_SPT_jet))
				{A[N][0][1]+=1;RecordIndividualSpectrums(it,weighting,hist_ST[0][1],hist_HT[0][1],hist_MET[0][1],hist_LPT[0][1]);}
			if(!(MET<CUT_MET) && !(HT<CUT_SPT_jet))
				{A[N][1][1]+=1;RecordIndividualSpectrums(it,weighting,hist_ST[1][1],hist_HT[1][1],hist_MET[1][1],hist_LPT[1][1]);}
			*/

		}//end of iterator

		std::cout << "\n" << processes[N+1] << ":" << std::endl;
		for(int i=0; i<2; i++){//for individual A, B, C, D
			for(int j=0; j<2; j++){
				int k = 0; if(i==0) k=1-j; else k=j;// A(0,1) B(0,0) C(1,0) D(1,1)
				//Eff[k] = count[k] / (double)vec.size(); 
				//Err_Eff[k] = sqrt( Eff[k]*(1-Eff[k])/ (double)vec.size() );
				//Yield[k] = L*pb2fb*X[k]*Eff[k];
				//Err_Yield[k] = sqrt( pow(Err_X[k]/X[k],2) + pow(Err_Eff[k]/Eff[k],2) )*Yield[k];//relative_err times Yield

				Eff[N][i][k]   = A[N][i][k] / (double)vec.size();
				Yield[N][i][k] = L*pb2fb*X[N+1]*Eff[N][i][k]; 

				Err_A[N][i][k] 	   = sqrt( (double)vec.size()*Eff[N][i][k]*(1-Eff[N][i][k]) );
				Err_Eff[N][i][k]   = sqrt( Eff[N][i][k]*(1-Eff[N][i][k])/(double)vec.size() );
				Err_Yield[N][i][k] = sqrt( pow(Err_X[N+1]/X[N+1],2) + pow(Err_Eff[N][i][k]/Eff[N][i][k],2) )*Yield[N][i][k];//relative_err times Yield

				totA[i][k] += A[N][i][k];//this would be useless
				totYield[i][k] += Yield[N][i][k];
				
				std::cout << Label[2*i+k] << "   = " << Yield[N][i][k] << " +/- " << Err_Yield[N][i][k] << std::endl;
			}
		}
		for(int i=0; i<2; i++){//for comparison: whether Yield/B == D/C
			double err = sqrt( pow(Err_Yield[N][i][1]/Yield[N][i][0],2) + pow(Err_Yield[N][i][0]*Yield[N][i][1]/pow(Yield[N][i][0],2),2) );//error transfer
			std::cout << Label[2*i+1] << "/" << Label[2*i+0] << " = " << Yield[N][i][1]/ Yield[N][i][0] << " +/- " << err << std::endl;
		}

	}//end of N, which loops for all processes

	std::cout << "\n" << "tot. Entries = " << totYield[0][0]+totYield[0][1]+totYield[1][0]+totYield[1][1] << "\n";
	for(int i=0; i<2; i++){
		for(int j=0; j<2; j++){
			int k = 0; if(i==0) k=1-j; else k=j; // Yield(0,1) B(0,0) C(1,0) D(1,1)
			for(int N=0; N<NUM; N++) Err_totYield[i][k] += pow(Err_Yield[N][i][k],2);
			Err_totYield[i][k] = sqrt(Err_totYield[i][k]);
			std::cout << Label[2*i+k] << "   = " << totYield[i][k] << " +/- " << Err_totYield[i][k] << std::endl;
		}
	}
	for(int i=0; i<2; i++){
		double err = sqrt( pow(Err_totYield[i][1]/totYield[i][0],2) + pow(Err_totYield[i][0]*totYield[i][1]/pow(totYield[i][0],2),2) );//error transfer
		std::cout << Label[2*i+1] << "/" << Label[2*i+0] << " = " << totYield[i][1]/ totYield[i][0] << " +/- " << err << std::endl;
	}

	/*
	double contribution[4][2][2]={0.};
	for(int i=0; i<2; i++){
		for(int j=0; j<2; j++){
			contribution[0][i][j] = GetTotContent(hist_ST[i][j]);
			contribution[1][i][j] = GetTotContent(hist_HT[i][j]);
			contribution[2][i][j] = GetTotContent(hist_MET[i][j]);
			contribution[3][i][j] = GetTotContent(hist_LPT[i][j]);
		}
	}
	double TOT = contribution[0][0][0]+contribution[0][0][1]+contribution[0][1][0]+contribution[0][1][1];

	std::cout << "\n" << "tot. Entries = " << TOT << "\n";
	for(int i=0; i<2; i++){
		for(int j=0; j<2; j++){
			int k = 0; if(i==0) k=1-j; else k=j; // A(0,1) B(0,0) C(1,0) D(1,1)
			for(int N=0; N<NUM; N++) Err_totA[i][k] += pow(Err_A[N][i][k],2);
			Err_totA[i][k] = sqrt(Err_totA[i][k]);
			std::cout << Label[2*i+k] << "   = " << totA[i][k] << " +/- " << Err_totA[i][k] << std::endl;
		}
	}
	for(int i=0; i<2; i++){
		double err = sqrt( pow(Err_totA[i][1]/totA[i][0],2) + pow(Err_totA[i][0]*totA[i][1]/pow(totA[i][0],2),2) );//error transfer
		std::cout << Label[2*i+1] << "/" << Label[2*i+0] << " = " << totA[i][1]/ totA[i][0] << " +/- " << err << std::endl;
	}
	*/

	fout->Write();
	fout->Close();
	return 1;
}
