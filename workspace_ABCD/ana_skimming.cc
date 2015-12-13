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
#include <TCut.h>
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

const int Nbins = 20;
const int NUM_CUT = 9;
const double L=100.;
const double pb2fb=1000.;
const double L_MASS=60.;//mass lower bound for bonson candidates
const double U_MASS=110.;//mass upper bound for bonson candidates
char *cutName[NUM_CUT]={"Num_lep","Num_jet","PT_lep","PT_jet","tot_Lep_PT","HT","MET","ST","Num_boson"};

char *inputFile;
char *savingPath;
char *savingName;

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
//std::vector<MyEvent> vec;
void importEvents(std::vector<MyEvent> &vec, TString filename){
	std::cout<<"Loading: "<<filename<<"\t\t";

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

void GetParticleInfo(std::vector<MyEvent>::iterator it, char *name, double CUT_PT, double &Num, double &TotPT, double &TotPT_selected){
	if(name=="Lep"){
		for(int i=0; i < it->Electrons.size(); i++){
			if(it->Electrons[i].pt > CUT_PT) {Num += 1; TotPT_selected += it->Electrons[i].pt;}
			TotPT += it->Electrons[i].pt;
		}
		for(int i=0; i < it->Muons.size(); i++){
			if(it->Muons[i].pt > CUT_PT) {Num += 1; TotPT_selected += it->Muons[i].pt;}
			TotPT += it->Muons[i].pt;
		}
	} else if(name=="Jet"){
		for(int i=0; i < it->Jets.size(); i++){
			if(it->Jets[i].pt > CUT_PT) {Num += 1; TotPT_selected += it->Jets[i].pt;}
			TotPT += it->Jets[i].pt;
		}
	} else std::cout<<"Please check the code carefully."<<std::endl;
}

TCut DecideCut(std::pair<int,int> pair,TCut cut1, TCut cut2){
	if(pair == std::make_pair(0,0)) return cut1 && cut2;
	if(pair == std::make_pair(0,1)) return cut1 && !cut2;
	if(pair == std::make_pair(1,0)) return !cut1 && cut2;
	if(pair == std::make_pair(1,1)) return !cut1 && !cut2;
}

//##########################################################################################
//### THE MAIN FUNCTION
//##########################################################################################
int main(int argc, char* argv[]){
	inputFile  = argv[1];//eg. "/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_fireball_bp_1TeV.root"
	savingName = argv[2];//eg. "simulation_delphes_ppvv_skimmed.root"
	savingPath = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/output/skimmed"; 
	std::cout << "Processing File: " << inputFile << std::endl;

	std::vector<MyEvent> vec;
	importEvents(vec,inputFile);
	
	//###set selection cuts
	std::vector<double> factors;//to store the optimized selection cuts
	double list[9]={2,6,30,40,20,1200,50,1400,0};
	factors.insert(factors.begin(),list,list+9);
    //std::cout << "\nChosen factors:" << std::endl;
	//ListSelectionCuts(factors);//print out the best seclectionCuts!
	double CUT_Num_lep   = factors.at(0);
	double CUT_Num_jet   = factors.at(1);
	double CUT_PT_lep    = factors.at(2);
	double CUT_PT_jet    = factors.at(3);
	double CUT_TotLepPT  = factors.at(4);
	double CUT_TotJetPT  = factors.at(5);
	double CUT_MET     	 = factors.at(6);
	double CUT_ST 		 = factors.at(7);
	double CUT_Num_boson = factors.at(8);

	//###ABCD method ( w/o PTcut vs. w/ PTcut; A,B,C,D regions )
	TFile *fout = new TFile(Form("%s/%s",savingPath,savingName),"recreate");
	TH2D *hist_NjvsNl[2], *hist_HTvsMET[2];
	TH1D *histST[2][2][2], *histHT[2][2][2], *histMET[2][2][2], *histLPT[2][2][2];//Spectrum for Num_jet vs. Num_lep; [ori/cut][01][01]
	TH1D *histNumJet[2][2][2];//Spectrum for HT vs. MET
	char *CutTag[2] = {"woPTcut","PTcut"}, *Label[4] = {"B","A","C","D"}, *title;//A(0,1) B(0,0) C(1,0) D(1,1)
	for(int k=0; k<2; k++){
		for(int i=0; i<2; i++){
			for(int j=0; j<2; j++){
				title = Form("%s;GeV;Entries / 125 GeV",Label[2*i+j]);
				histST[k][i][j]	    = new TH1D(Form("histST_%s%d%d"	   ,CutTag[k],i,j),title,Nbins,0,5000); histST[k][i][j]	->Sumw2();
				histHT[k][i][j]	    = new TH1D(Form("histHT_%s%d%d"	   ,CutTag[k],i,j),title,Nbins,0,5000); histHT[k][i][j]	->Sumw2();
				histLPT[k][i][j]    = new TH1D(Form("histLPT_%s%d%d"   ,CutTag[k],i,j),title,Nbins,0,1200); histLPT[k][i][j]->Sumw2();
				histMET[k][i][j]    = new TH1D(Form("histMET_%s%d%d"   ,CutTag[k],i,j),title,Nbins,0,1200); histMET[k][i][j]->Sumw2();
				title = Form("%s;# of jets;Entries",Label[2*i+j]);
				histNumJet[k][i][j] = new TH1D(Form("histNumJet_%s%d%d",CutTag[k],i,j),title,Nbins,0,25)  ; histNumJet[k][i][j]->Sumw2();
			}
		}
		hist_NjvsNl[k]  = new TH2D(Form("hist_NjvsNl_%s",CutTag[k]) ,";# of lep;# of jet",7,0,7,25,0,25);				 hist_NjvsNl[k] ->Sumw2();
		hist_HTvsMET[k] = new TH2D(Form("hist_HTvsMET_%s",CutTag[k]),";MET (60 GeV); HT (250 GeV)",20,0,1200,20,0,5000); hist_HTvsMET[k]->Sumw2();
	}//end of k

	for( std::vector<MyEvent>::iterator it=vec.begin(); it!=vec.end(); it++){
		double Num_lep = 0, Num_jet = 0, MET = 0;
		double TotLepPT = 0, TotLepPT_selected = 0;
		double TotJetPT = 0, TotJetPT_selected = 0;
		double ST = 0, ST_selected = 0;
		
		//lep, jet
		GetParticleInfo(it, "Lep", CUT_Num_lep, Num_lep, TotLepPT, TotLepPT_selected);
		GetParticleInfo(it, "Jet", CUT_Num_jet, Num_jet, TotJetPT, TotJetPT_selected);
		//met
		MET = it->METs[0].pt;
		//ST
		ST = TotLepPT + TotJetPT + MET;
		ST_selected = TotLepPT_selected + TotJetPT_selected + MET;


		hist_NjvsNl[0]->Fill(it->Electrons.size()+it->Muons.size(),it->Jets.size());
		hist_HTvsMET[0]->Fill(MET,TotJetPT);
		hist_NjvsNl[1]->Fill(Num_lep,Num_jet);
		hist_HTvsMET[1]->Fill(MET,TotJetPT_selected);


		if(it->Electrons.size()+it->Muons.size() < CUT_Num_lep && it->Jets.size() < CUT_Num_jet)
			{histST[0][0][0]->Fill(ST); histHT[0][0][0]->Fill(TotJetPT); histLPT[0][0][0]->Fill(TotLepPT); histMET[0][0][0]->Fill(MET);}
		if(!(it->Electrons.size()+it->Muons.size() < CUT_Num_lep) && it->Jets.size() < CUT_Num_jet)
			{histST[0][1][0]->Fill(ST); histHT[0][1][0]->Fill(TotJetPT); histLPT[0][1][0]->Fill(TotLepPT); histMET[0][1][0]->Fill(MET);}
		if(it->Electrons.size()+it->Muons.size() < CUT_Num_lep && !(it->Jets.size() < CUT_Num_jet))
			{histST[0][0][1]->Fill(ST); histHT[0][0][1]->Fill(TotJetPT); histLPT[0][0][1]->Fill(TotLepPT); histMET[0][0][1]->Fill(MET);}
		if(!(it->Electrons.size()+it->Muons.size() < CUT_Num_lep) && !(it->Jets.size() < CUT_Num_jet))
			{histST[0][1][1]->Fill(ST); histHT[0][1][1]->Fill(TotJetPT); histLPT[0][1][1]->Fill(TotLepPT); histMET[0][1][1]->Fill(MET);}
			
		if(Num_lep < CUT_Num_lep && Num_jet < CUT_Num_jet)
			{histST[1][0][0]->Fill(ST_selected); histHT[1][0][0]->Fill(TotJetPT_selected); histLPT[1][0][0]->Fill(TotLepPT_selected); histMET[1][0][0]->Fill(MET);}
		if(!(Num_lep < CUT_Num_lep) && Num_jet < CUT_Num_jet)
			{histST[1][1][0]->Fill(ST_selected); histHT[1][1][0]->Fill(TotJetPT_selected); histLPT[1][1][0]->Fill(TotLepPT_selected); histMET[1][1][0]->Fill(MET);}
		if(Num_lep < CUT_Num_lep && !(Num_jet < CUT_Num_jet))
			{histST[1][0][1]->Fill(ST_selected); histHT[1][0][1]->Fill(TotJetPT_selected); histLPT[1][0][1]->Fill(TotLepPT_selected); histMET[1][0][1]->Fill(MET);}
		if(!(Num_lep < CUT_Num_lep) && !(Num_jet < CUT_Num_jet))
			{histST[1][1][1]->Fill(ST_selected); histHT[1][1][1]->Fill(TotJetPT_selected); histLPT[1][1][1]->Fill(TotLepPT_selected); histMET[1][1][1]->Fill(MET);}
			
			
		if(MET < CUT_MET && ST < CUT_ST)		histNumJet[0][0][0]->Fill(it->Jets.size());
		if(!(MET < CUT_MET) && ST < CUT_ST)		histNumJet[0][1][0]->Fill(it->Jets.size());
		if(MET < CUT_MET && !(ST < CUT_ST))		histNumJet[0][0][1]->Fill(it->Jets.size());
		if(!(MET < CUT_MET) && !(ST < CUT_ST))	histNumJet[0][1][1]->Fill(it->Jets.size());
			
		if(MET < CUT_MET && ST_selected < CUT_ST)		histNumJet[1][0][0]->Fill(Num_jet);
		if(!(MET < CUT_MET) && ST_selected < CUT_ST)	histNumJet[1][1][0]->Fill(Num_jet);
		if(MET < CUT_MET && !(ST_selected < CUT_ST))	histNumJet[1][0][1]->Fill(Num_jet);
		if(!(MET < CUT_MET) && !(ST_selected < CUT_ST))	histNumJet[1][1][1]->Fill(Num_jet);
			
		//TCut cut, cut1[2], cut2[2], cut3[2], cut4[2];
		//cut1[0] = "it->Electrons.size()+it->Muons.size() < CUT_Num_lep";
		//cut2[0] = "it->Jets.size() < CUT_Num_jet";
		//cut1[1] = "Num_lep < CUT_Num_lep";
		//cut2[1] = "Num_jet < CUT_Num_jet";
		//cut3[0] = "MET < CUT_MET";
		//cut4[0] = "ST < CUT_ST";
		//cut3[1] = "MET < CUT_MET";
		//cut4[1] = "ST_selected < CUT_ST";
		
		/*
		for(int k=0; k<2; k++){
			if(k==0) {hist_NjvsNl[k]->Fill(it->Electrons.size()+it->Muons.size(),it->Jets.size()); hist_HTvsMET[k]->Fill(MET,TotJetPT);}
			if(k==1) {hist_NjvsNl[k]->Fill(Num_lep,Num_jet); hist_HTvsMET[k]->Fill(MET,TotJetPT_selected);}
			for(int i=0; i<2; i++){
				for(int j=0; j<2; j++){
					cut = DecideCut(std::make_pair(i,j),cut1[k],cut2[k]);
					if(cut){
						//Spectrum for Num_jet vs. Num_lep; [ori/cut][01][01]
				 		histST[k][i][j]	->Fill(ST);
				 		//histHT[k][i][j]	->Fill(TotLepPT);
				 		//histLPT[k][i][j]->Fill(TotJetPT);
				 		//histMET[k][i][j]->Fill(MET);
					}
					//cut = DecideCut(std::make_pair(i,j),cut3[k],cut4[k]);
					//if(cut){
					//	//Spectrum for HT vs. MET
				 	//	if(k==0) histNumJet[k][i][j]->Fill(it->Jets.size());
				 	//	if(k==1) histNumJet[k][i][j]->Fill(Num_jet);
					//}
				}
			}
		}//end of k
		*/

	}//end of event iterator

	fout->Write();
	fout->Close();
	return 1;
}
