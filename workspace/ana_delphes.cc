// PID code/particle
// 1/d,2/u,3/s,4/c,5/b,6/t,7/b',8/t'
// 11/e-,12/nu_e,13/mu-,14/nu_mu,15/tau-,16/nu_tau
// 21/g, 22/gamma, 23/z, 24/W+, 25/h, 211/pi+, 2212/proton
// status==1 particles are those that should be processed directly by a detetor simulation
// status==2 particles are unstable, such as a pi0 that will later decay to gamma
// status==3 particles label the beam and the hard process
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
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TPad.h>
#include <TRandom3.h>

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"

using namespace std;
using namespace TMVA;

TRandom3 rnd(39);
const int num_bg = 8;
const int NUM = num_bg+1;
const int NUM_hist = 13;
const int NUM_CUT  = 9;
const double L=100.;
const double pb2fb=1000.;
const double L_MASS=60.;//mass lower bound for bonson candidates
const double U_MASS=110.;//mass upper bound for bonson candidates

double X[NUM], Err_X[NUM];//to be determined by SIG (the Xsections of all considered processes)
//double X_sig[2] = {0.0324, 0.0001134}; //1TeV, 2TeV
//double Err_X_sig[2] = {4.8e-05, 1.6e-07}; //1TeV, 2TeV
//double X_bg[8] = {472.1, 130.3, 1.047, 0.3969, 0.01295, 0.009666, 0.002586, 0.0001468}; //pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
//double Err_X_bg[8] = {11, 1.3, 0.018, 0.046, 0.0017, 3.1e-05, 2.2e-05, 3.2e-07}; //pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
double X_sig[2] = {0.051, 0.000177}; //1TeV, 2TeV
double Err_X_sig[2] = {0.001, 4e-06}; //1TeV, 2TeV
double X_bg[8] = {815.96, 192.4, 1.77, 0.621, 0.0219, 0.0164, 0.00405, 0.000249}; //pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
double Err_X_bg[8] = {45.51, 4.1, 0.10, 0.073, 0.0013, 9e-04, 9e-05, 3.4e-05}; //pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
//double X_bg[8] = {624.4, 130.3, 1.047, 0.3969, 0.0342, 0.009666, 0.002586, 0.000148}; //pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt

char *histName[13]={"Num_lep","Num_jet","PT_lep","PT_jet","tot_Lep_PT","HT","MET","ST","Num_boson","M_boson","PT_chosenJet","Eta_chosenJet","Phi_chosenJet"};
char *processes[NUM] = {"fireball", "pptt", "ppvv", "ppvtt", "ppvvv", "ppvvtt", "pptttt", "ppvvvv", "ppvvvtt"};
char *units[NUM_hist] = {"# of lep", "# of jet", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "# of boson", "GeV", "GeV", "eta", "rad."};
//int color[NUM] = {kRed, kOrange+3, kOrange+1, kBlue, kCyan, kGreen+2, kGreen, kViolet+2, kMagenta};
int color[NUM] = {kRed, kOrange+3, kOrange+1, kBlue, kCyan, kGreen+2, kGreen, kMagenta+2, kMagenta};
int timesStyle[4] = {20,22,24,26};//20 21 22 23:(solid) circle square triangle anti-triangle; 24 25 26 27 (hollow)~~~~
int timesColor[5] = {kRed, kOrange+1, kGreen+1, kBlue-2, kMagenta+2};
//int timesColor[5] = {kRed, kOrange+1, kBlue, kGreen+2, kViolet+2};
int CutLineColor = kBlack;

string METHOD;//record the inpurt method TMVA/MyGA/NoGA
int SIG;//record the input value
char *savingPath;
char *suffix;

//for finding segmentation error
void sighandler(int sig){
	while(1);
}
//signal(SIGSEGV,sighandler);


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

int COUNTER = 0;
TH1F *hist = new TH1F("hist",";significance;entries",20,0,2);
//#############################################
//### Define Fitness for TMVA GA
//#############################################
class MyFitness : public IFitterTarget{
public:
	MyFitness() : IFitterTarget() {
	}

	Double_t EstimatorFunction( std::vector<Double_t> & factors ){
		double CUT_Num_lep   = factors.at(0);
		double CUT_Num_jet   = factors.at(1);
		double CUT_PT_lep    = factors.at(2);
		double CUT_PT_jet    = factors.at(3);
		double CUT_SPT_lep   = factors.at(4);
		double CUT_SPT_jet   = factors.at(5);
		double CUT_TOT_MET   = factors.at(6);
		double CUT_ST 		 = factors.at(7);
		double CUT_Num_boson = factors.at(8);
		
		double eff_sig, eff_bg[num_bg], yield_sig, yield_bg;
		double S = 0., B[num_bg];
		for(int i=0; i<num_bg; i++)	B[i]=0.;
		
		std::vector<MyEvent> vec;
		//loop over signal & bg
		for(int k=0; k<NUM; k++){
			if(k==0) vec = vec_sig;
			else 	 vec = vec_bg[k-1];

			for( std::vector<MyEvent>::iterator it=vec.begin(); it!=vec.end(); it++){
				
				if((it->Electrons.size()+it->Muons.size()) < CUT_Num_lep) continue;
				if((it->Jets.size()) < CUT_Num_jet) continue;

				double Num_lep = 0, Num_jet = 0, SPT_lep = 0, SPT_jet = 0, TOT_met = 0, ST = 0;
				
				//boson
				if(it->H_Bosons.size() < CUT_Num_boson) continue;
				//lep
				for(int i=0; i<it->Electrons.size(); i++){
					if(it->Electrons[i].pt > CUT_PT_lep) Num_lep += 1;
					SPT_lep += it->Electrons[i].pt;
				}
				for(int i=0; i<it->Muons.size(); i++){
					if(it->Muons[i].pt > CUT_PT_lep) Num_lep += 1;
					SPT_lep += it->Muons[i].pt;
				}
				if(Num_lep < CUT_Num_lep) continue;
				if(SPT_lep < CUT_SPT_lep) continue;
				//met
				for(int i=0; i<it->METs.size(); i++) TOT_met += it->METs[i].pt;
				if(TOT_met < CUT_TOT_MET) continue;
				//jet
				for(int i=0; i<it->Jets.size(); i++){
					if(it->Jets[i].pt > CUT_PT_jet) Num_jet += 1;
					SPT_jet += it->Jets[i].pt;
				}
				if(Num_jet < CUT_Num_jet) continue;
				if(SPT_jet < CUT_SPT_jet) continue;
				//ST
				ST = SPT_lep + SPT_jet + TOT_met;
				if(ST < CUT_ST) continue;

				if(k==0) S +=1. ;
				else 	 B[k-1] += 1.;
			}//end of event loop

			if(k==0){
				eff_sig = S / (double)vec_sig.size(); 
				yield_sig = L*pb2fb*X_sig[SIG-1]*eff_sig; 
			} else{
				eff_bg[k-1] = B[k-1] / (double)vec_bg[k-1].size(); 
				yield_bg += L*pb2fb*X_bg[k-1]*eff_bg[k-1]; 
			}
		}//end of k

		//double significance = yield_sig / sqrt( yield_sig + yield_bg );
		//if((yield_sig+yield_bg)!=0)	return -significance;
		double significance;
		if(yield_bg>1e-20) significance = sqrt( 2*( (yield_sig + yield_bg)*TMath::Log(1+(yield_sig/yield_bg)) - yield_sig ) );
		else significance = -1;
		hist->Fill(significance);

		if(METHOD=="TMVA"){
			COUNTER += 1;
			std::cout<<"("<<COUNTER<<")"<<"\tsig = "<<significance<<std::endl;
		}

		//if((yield_bg)!=0){
		//	hist->Fill(significance);
		//	return -significance;
		//}else{//denominator==0
		//	hist->Fill(-1);
		//	return 1;
		//}

		return -significance;

	}//end of EstimatorFunction
};

//#############################################
//### Miscellaneous Functions
//#############################################
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
void PrintSelectionCuts(std::vector<double> factors){
	std::cout<<"factors are:";
	for( std::vector<Double_t>::iterator it = factors.begin(); it!=factors.end(); it++)
		std::cout<<' '<<*it; 
	std::cout<<'\n';
}
void ListSelectionCuts(std::vector<double> factors){
	int n = 0;
    for( std::vector<Double_t>::iterator it = factors.begin(); it<factors.end(); it++ ){
        std::cout << setw(11) <<histName[n] << ": " << (*it) << std::endl;
        n++;
    }
}

//#############################################
//### Define the MyGA Functions
//#############################################
//const int NCHROMOS = 500;
//const int NGEN = 100;
//const int Ntimes = 5;
const int NCHROMOS = 10;
const int NGEN = 5;
const int Ntimes = 3;
const int SpreadFactor = 3;//1 = no mutation effect
const int SamplingPeriod = 5;
//std::vector<Double_t> chromo[NCHROMOS*2];
std::vector< std::vector<double> > chromo;
double MySig[NCHROMOS*2];

//void PrepareChromosomes(int NCHROMOS, std::vector<double> &chromo[NCHROMOS*2]){
void PrepareChromosomes(){
	for(int i=0;i<NCHROMOS;i++){
		std::vector<double> DNAsegment;
		DNAsegment.push_back( (double)(int)rnd.Uniform(2,10) );//Num_lep
		DNAsegment.push_back( (double)(int)rnd.Uniform(0,15) );//Num_jet
		DNAsegment.push_back( (double)(int)(rnd.Uniform(20,80)/5)*5 );//PT_lep
		DNAsegment.push_back( (double)(int)(rnd.Uniform(30,80)/5)*5 );//PT_jet
		DNAsegment.push_back( (double)(int)(rnd.Uniform(0,500)/20)*20  );//SPT_lep
		DNAsegment.push_back( (double)(int)(rnd.Uniform(0,1500)/20)*20 );//SPT_jet
		DNAsegment.push_back( (double)(int)(rnd.Uniform(20,800)/20)*20 );//MET
		DNAsegment.push_back( (double)(int)(rnd.Uniform(0,2500)/50)*50 );//ST
		DNAsegment.push_back( (double)(int)rnd.Uniform(0,10) );//Num_boson
		chromo.push_back(DNAsegment);
	}
}
//void ApplyCrossOver(int NCHROMOS, std::vector<double> &chromo[NCHROMOS*2]){
void ApplyCrossOver(){
	for(int i=NCHROMOS;i<NCHROMOS*2;i++){
		int s1=(int)rnd.Uniform(0.,NCHROMOS);
		int s2=(int)rnd.Uniform(0.,NCHROMOS);
		std::vector<double> DNAsegment;
		for(int j=0;j<9;j++){
			if(rnd.Uniform()>0.5) DNAsegment.push_back(chromo.at(s1).at(j));
			else DNAsegment.push_back(chromo.at(s2).at(j));
		} 
		chromo.push_back(DNAsegment);
	}
}
//void ApplyMutation(int NCHROMOS, std::vector<double> &chromo[NCHROMOS*2]){
void ApplyMutation(){
	for(int i=NCHROMOS;i<NCHROMOS*2;i++){
		if(rnd.Uniform()<0.60){//fluctuate 
			chromo[i].at(0)+=(double)((int)rnd.Gaus(0,1*SpreadFactor));
			chromo[i].at(1)+=(double)((int)rnd.Gaus(0,1*SpreadFactor));
			chromo[i].at(2)+=(double)((int)(rnd.Gaus(0,5*SpreadFactor)/5)*5);
			chromo[i].at(3)+=(double)((int)(rnd.Gaus(0,5*SpreadFactor)/5)*5);
			chromo[i].at(4)+=(double)((int)(rnd.Gaus(0,20*SpreadFactor)/20)*20);
			chromo[i].at(5)+=(double)((int)(rnd.Gaus(0,20*SpreadFactor)/20)*20);
			chromo[i].at(6)+=(double)((int)(rnd.Gaus(0,20*SpreadFactor)/20)*20);
			chromo[i].at(7)+=(double)((int)(rnd.Gaus(0,50*SpreadFactor)/50)*50);
			chromo[i].at(8)+=(double)((int)rnd.Gaus(0,2*SpreadFactor));
			for(int j=0;j<9;j++){
				//if(j==0 || j==1 || j==8) chromo[i].at(j)=(double)(int)chromo[i].at(j);
				if(j==0 || j==2 || j==3 || j==6){
					if(j==0 && chromo[i].at(j)<2) chromo[i].at(j)=2;//Num_lep
					if(j==2 && chromo[i].at(j)<20) chromo[i].at(j)=20;//PT_lep
					if(j==3 && chromo[i].at(j)<30) chromo[i].at(j)=30;//PT_jet
					if(j==6 && chromo[i].at(j)<20) chromo[i].at(j)=20;//MET
				}
				else{ if(chromo[i].at(j)<0) chromo[i].at(j)=0;}
			}
		} 
	}
}
void MyGA(TNtuple *&SigGen, TH1D *&SigOpt){
	//###Initiation
	MyFitness *Fitness;
	PrepareChromosomes();
	for(int i=0;i<NCHROMOS;i++)
		MySig[i] = -(Fitness->MyFitness::EstimatorFunction(chromo[i]));
		//MySig[i] = fabs(Fitness->MyFitness::EstimatorFunction(chromo[i]));
	
	//###Evolution loop
	for(int gen=0;gen<NGEN;gen++){
		printf("generation=%d\n",gen+1);
		ApplyCrossOver();
		ApplyMutation();
		for(int i=NCHROMOS;i<NCHROMOS*2;i++)
			MySig[i] = -(Fitness->MyFitness::EstimatorFunction(chromo[i]));
			//MySig[i] = fabs(Fitness->MyFitness::EstimatorFunction(chromo[i]));

		//sorting
		for(int i=0;i<NCHROMOS*2;i++){
			for(int j=i+1;j<NCHROMOS*2;j++){
				if(MySig[j]>MySig[i]){
					double s=MySig[j];
					MySig[j]=MySig[i];
					MySig[i]=s;
					for(int k=0;k<9;k++){
						double t=chromo[j].at(k);
						chromo[j].at(k)=chromo[i].at(k);
						chromo[i].at(k)=t;
					}
				}
			}
		}
		//calculate the ave. and rms.of MySig
		int n=0;
		double s_ave=0, s_rms=0;
		for(int i=0;i<NCHROMOS;i++){
			if(MySig[i]!=-1){
				s_ave+=MySig[i];
				n+=1;
			}
		}
		s_ave/=n;
		SigGen->Fill(gen+1,s_ave,MySig[0]);//generation, ave, best

		for(int i=0;i<NCHROMOS;i++)
			if(MySig[i]!=-1) s_rms+=pow((MySig[i]-s_ave),2);
		s_rms/=n;
		s_rms=pow(s_rms,0.5);

		std::cout<<"nChromos = "<< n <<"\ts_ave = "<< s_ave <<"\ts_rms = "<< s_rms << std::endl;
		for(int i=0; i<5; i++){
			std::cout<<"order:"<<i<<"\tsig = "<<MySig[i]<<std::endl;
			PrintSelectionCuts(chromo[i]);
			std::cout<<std::endl;
		}
		if(s_rms<0.001*s_ave) break;
	}//end of evolution loop

	SigOpt->Fill(MySig[0]);

}

//#############################################
//### THE MAIN FUNCTION
//#############################################
int main(int argc, char* argv[]){
	std::cout << "The chosen bp mass is " << argv[2] << " TeV" << std::endl;
	METHOD = argv[1];// TMVA or MyGA or NoGA
	SIG = atoi(argv[2]);
	suffix = argv[3];
	//if(SIG==1) savingPath = Form("/afs/cern.ch/user/y/ykao/work/fireball/03output/1TeV%s",suffix); 
	//if(SIG==2) savingPath = Form("/afs/cern.ch/user/y/ykao/work/fireball/03output/2TeV%s",suffix); 
	if(SIG==1) savingPath = "/afs/cern.ch/user/y/ykao/work/fireball/03output/1TeV"; 
	if(SIG==2) savingPath = "/afs/cern.ch/user/y/ykao/work/fireball/03output/2TeV"; 

	//pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
	if(SIG==1) importEvents(vec_sig,"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_fireball_bp_1TeV.root");
	if(SIG==2) importEvents(vec_sig,"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_fireball_bp_2TeV.root");
	importEvents(vec_bg[0],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_pptt_jet_matching.root");
	importEvents(vec_bg[1],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvv_jet_matching.root");
	importEvents(vec_bg[2],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvtt_jet_matching.root");
	importEvents(vec_bg[3],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvv_jet_matching.root");
	importEvents(vec_bg[4],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvtt_jet_matching.root");
	importEvents(vec_bg[5],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_pptttt.root");
	importEvents(vec_bg[6],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvvv_jet_matching.root");
	importEvents(vec_bg[7],"/afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvvtt.root");

	//set cross section list
	for(int i=0; i<NUM; i++){
		if(i==0) {X[i] = X_sig[SIG-1]; Err_X[i] = Err_X_sig[SIG-1];}
		else {X[i] = X_bg[i-1]; Err_X[i] = Err_X_bg[i-1];}
	}
	
	std::vector<double> factors;//to store the optimized selection cuts
	//#############################################
 	//### Optimize the selection cuts with TMVA GA
	//#############################################
	//if(!METHOD.compare("TMVA")){};
	if(METHOD=="TMVA"){
		std::vector<Interval*> ranges;
		ranges.push_back( new Interval(2.,10.,9) );//di-lepton
		ranges.push_back( new Interval(0.,15.,16) );
		ranges.push_back( new Interval(20.,80.,13) );
		ranges.push_back( new Interval(30.,100.,15) );
		ranges.push_back( new Interval(0.,500.,26) );
		ranges.push_back( new Interval(0.,1500.,31) );
		ranges.push_back( new Interval(20.,800.,40) );
		ranges.push_back( new Interval(0.,2500.,51) );
		ranges.push_back( new Interval(0.,10.,11) );

		std::cout<<"\nInitial ranges:"<<std::endl;
		for(std::vector<Interval*>::iterator it=ranges.begin(); it != ranges.end(); it++){
			std::cout<<"range:"<< (*it)->GetMin() << "   " << (*it)->GetMax()<<std::endl;
		}
		
		IFitterTarget *myFitness = new MyFitness();
		GeneticAlgorithm mg( *myFitness, 5, ranges);
    	
		#define CONVSTEPS 5 
		#define CONVCRIT 0.0001
		#define SCSTEPS 10 
		#define SCRATE 5 
		#define SCFACTOR 0.8

    	do {
    	    mg.Init();
    	    mg.CalculateFitness();
    	    mg.GetGeneticPopulation().Print(0);
    	    mg.GetGeneticPopulation().TrimPopulation();
    	    mg.SpreadControl( SCSTEPS, SCRATE, SCFACTOR );
    	} while (!mg.HasConverged( CONVSTEPS, CONVCRIT ));
    	
    	GeneticGenes* genes = mg.GetGeneticPopulation().GetGenes( 0 );
		factors = genes->GetFactors();
    	int n = 0;
    	std::cout << "\nBest factors:" << std::endl;
		ListSelectionCuts(factors);//print out the Best selectionCuts!
    	std::cout << "\nCounter = " << COUNTER << std::endl;
	}//end of TMVA

	//#############################################
 	//### Optimize the selection cuts with MyGA
	//#############################################
	if(METHOD=="MyGA"){
		TCanvas *canMyGA = new TCanvas("canMyGA","",800,600);
		TLegend *legend_best = new TLegend(0.72,0.5,0.81,0.9);
		TLegend *legend_ave  = new TLegend(0.81,0.5,0.9,0.9);
		//legend -> SetNColumns(2);
		TNtuple *SigGen[Ntimes], *SigGenCopy[Ntimes];// new TNtuple("SigGen","SigGen","gen:ave:best");
		TH2D *SigGenHist = new TH2D("SigGenHist","Evolution;gen;significance",NGEN+10,0,(double)(NGEN+10),12,0,1.2);
		if(SIG==1) SigGenHist -> SetBins(NGEN+10,0,(double)(NGEN+10),18,0,18);
		TH1D *SigOpt; 
		if(SIG==1) SigOpt = new TH1D("SigOpt","SigOpt",200,0,20);
		if(SIG==2) SigOpt = new TH1D("SigOpt","SigOpt",120,0,1.2);
		
		//Loop over several times for validation of GA & Draw the BEST significance in evolution
		for(int i_time=0; i_time<Ntimes; i_time++){
			SigGen[i_time] = new TNtuple(Form("SigGen_%d",i_time),"SigGen","gen:ave:best");

			MyGA(SigGen[i_time],SigOpt);
			SigGenCopy[i_time] = (TNtuple*) SigGen[i_time]->Clone();

			SigGen[i_time]->SetMarkerSize(0.8);
			SigGen[i_time]->SetMarkerColor(timesColor[i_time % SamplingPeriod]);
			SigGen[i_time]->SetMarkerStyle(timesStyle[(int)(i_time / SamplingPeriod)]);
			SigGen[i_time]->SetLineColor(timesColor[i_time % SamplingPeriod]);
			SigGen[i_time]->SetLineStyle(1);
			SigGen[i_time]->SetLineWidth(2);
			legend_best->AddEntry(SigGen[i_time],Form("Best_{trial %d}",i_time+1),"lp");//lepf
			if(i_time==0) SigGenHist->Draw();
			SigGen[i_time]->Draw("best:gen","","lp,same");
			if(i_time==Ntimes-1) legend_best->Draw("same");

			SigGenCopy[i_time]->SetMarkerSize(0.8);
			SigGenCopy[i_time]->SetMarkerColor(timesColor[i_time % SamplingPeriod]);
			SigGenCopy[i_time]->SetMarkerStyle(timesStyle[(int)(i_time / SamplingPeriod)+2]);
			SigGenCopy[i_time]->SetLineColor(timesColor[i_time % SamplingPeriod]);
			SigGenCopy[i_time]->SetLineStyle(2);
			SigGenCopy[i_time]->SetLineWidth(1);
			legend_ave->AddEntry(SigGenCopy[i_time],Form("Ave_{trial %d}",i_time+1),"lp");//lepf
			SigGenCopy[i_time]->Draw("ave:gen","","lp,same");
			if(i_time==Ntimes-1) legend_ave->Draw("same");
		}
		gStyle->SetOptStat(0);
		canMyGA->SaveAs(Form("%s%s/rootfile/SigGen.root",savingPath,suffix));
		canMyGA->SaveAs(Form("%s%s/png/SigGen.png",savingPath,suffix));

		gStyle->SetOptStat(1);
		SigOpt->SetLineWidth(2);
		SigOpt->Draw();
		canMyGA->SaveAs(Form("%s%s/rootfile/SigOpt.root",savingPath,suffix));
		canMyGA->SaveAs(Form("%s%s/png/SigOpt.png",savingPath,suffix));
		
		factors = chromo[0];
    	std::cout << "\nBest factors:" << std::endl;
		ListSelectionCuts(factors);//print out the best seclectionCuts!
	}//end of MyGA

	if(METHOD=="NoGA"){//For ploting certain set of SelectionCuts
		double list01[9]={2,6,30,40,20,1200,50,1400,0};
		double list02[9]={2,6,30,40,20,1200,50,1400,0};
		std::vector<double>::iterator it;//to store the optimized selection cuts
		if(SIG==1) factors.insert(factors.begin(),list01,list01+9);
		if(SIG==2) factors.insert(factors.begin(),list02,list02+9);
    	std::cout << "\nChosen factors:" << std::endl;
		ListSelectionCuts(factors);//print out the best seclectionCuts!
	}
	
	//#############################################
	//### Ploting with the optimized selection cuts
	//#############################################
	TFile *fout = new TFile(Form("%s/result.root",savingPath),"recreate");
	TCanvas *can = new TCanvas("can","",800,600);
	SetStyle();
	TH1D *hist_Ori[NUM_hist][NUM];//before cuts 
	TH1D *hist_Cut[NUM_hist][NUM];//after cuts 
	TH1D *hist_Nm1_Plot[NUM_hist][NUM];//N-1 plots 
	TLine *CutLine[NUM_CUT];//the cut lines
	THStack *hist_Stack_Ori[NUM_CUT];//before cuts 
	THStack *hist_Stack_Cut[NUM_CUT];//after cuts 
	THStack *hist_Stack_Copy[NUM_CUT];//before cuts 

	int    n_bin_ori[2][NUM_hist]={{10,40,50,50,50,50,50,50,12,U_MASS-L_MASS,60,16,16},
								   {10,40,50,50,50,35,50,40,12,U_MASS-L_MASS,100,16,16}};
	double l_bin_ori[2][NUM_hist]={{0,0,0,0,0,0,0,0,0,L_MASS,0,-8,-8},
								   {0,0,0,0,0,0,0,0,0,L_MASS,0,-8,-8}};
	double h_bin_ori[2][NUM_hist]={{10,40, 800,1400, 600,5000,1200,7500,12,U_MASS,1200,8,8},
								   {10,40,1000,2000,1000,7000,1500,8000,12,U_MASS,1600,8,8}};
	int    n_bin_cut[2][NUM_hist]={{ 8,40,50,50,50,50,50,50,12,U_MASS-L_MASS,60,16,16},
								   {12,40,50,50,50,50,50,50,12,U_MASS-L_MASS,100,16,16}};
	double l_bin_cut[2][NUM_hist]={{0,0,0,0,0,0,0,0,0,L_MASS,0,-8,-8},
								   {0,0,0,0,0,0,0,0,0,L_MASS,0,-8,-8}};
	double h_bin_cut[2][NUM_hist]={{ 8,40, 800,1400, 600,5000,1200,7500,12,U_MASS,1200,8,8},
								   {12,40,1000,2000,1000,5500,1500,8000,12,U_MASS,1600,8,8}};
	double scaleBosonMass, scaleChosenJetEta, scaleChosenJetPhi, scaleLepSPT;//hist->SetMaximum(scale)
	if(SIG==1)	{scaleBosonMass = 20000, scaleChosenJetEta = 20000, scaleChosenJetPhi = 200000, scaleLepSPT = 1000;}
	if(SIG==2)  {scaleBosonMass = 500, scaleChosenJetEta = 10000, scaleChosenJetPhi = 5000, scaleLepSPT = 1000;}

	//printf("4.0: %s\n",savingPath);
	for(int k=0; k<NUM_hist; k++){//loop over physical quantities
		//printf("4.0.%d: %s\n",k,savingPath);
		for(int i=0; i<NUM; i++){//loop over processes
			hist_Ori[k][i]       = new TH1D(Form("%s_Ori_%s"    ,processes[i],histName[k]),"",n_bin_ori[SIG-1][k] ,l_bin_ori[SIG-1][k] ,h_bin_ori[SIG-1][k]);
			hist_Cut[k][i]       = new TH1D(Form("%s_Cut_%s"    ,processes[i],histName[k]),"",n_bin_ori[SIG-1][k] ,l_bin_ori[SIG-1][k] ,h_bin_ori[SIG-1][k]);
			hist_Nm1_Plot[k][i]  = new TH1D(Form("%s_N-1Plot_%s",processes[i],histName[k]),"",n_bin_ori[SIG-1][k] ,l_bin_ori[SIG-1][k] ,h_bin_ori[SIG-1][k]);
		}
		if(k<NUM_CUT){
			if(k==2 || k==3 || k==4 || k==5 || k==6 || k==7){
			double binWidth_ori = (h_bin_ori[SIG-1][k]-l_bin_ori[SIG-1][k])/n_bin_ori[SIG-1][k];
			double binWidth_cut = (h_bin_cut[SIG-1][k]-l_bin_cut[SIG-1][k])/n_bin_cut[SIG-1][k];
			hist_Stack_Ori[k] =new THStack(Form("Ori_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %.0f %s",binWidth_ori,units[k])));
			hist_Stack_Copy[k]=new THStack(Form("Copy_Stack_%s",histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %.0f %s",binWidth_ori,units[k])));
			hist_Stack_Cut[k] =new THStack(Form("Cut_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %.0f %s",binWidth_cut,units[k])));
			}
			else{
			hist_Stack_Ori[k] = new THStack(Form("Ori_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %s",units[k])));
			hist_Stack_Copy[k]= new THStack(Form("Copy_Stack_%s",histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %s",units[k])));
			hist_Stack_Cut[k] = new THStack(Form("Cut_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %s",units[k])));
			}
		}
	}
	
	//printf("4.1: %s\n",savingPath);
	
	//### Loop over "Full Cut" & "N-1 Cuts"
	double CUT[NUM_CUT];
	for(int N=0; N<NUM_CUT+1; N++){
		if(N==0) for(int k=0; k<NUM_CUT; k++) CUT[k] = factors.at(k);//for full cut
		else{//for N-1 Cuts
			for(int k=0; k<NUM_CUT; k++){
				if(k==N-1) CUT[k] = 0.;
				else CUT[k] = factors.at(k);
			}
		}
		double CUT_Num_lep   = CUT[0];
		double CUT_Num_jet   = CUT[1];
		double CUT_PT_lep    = CUT[2];
		double CUT_PT_jet    = CUT[3];
		double CUT_SPT_lep   = CUT[4];
		double CUT_SPT_jet   = CUT[5];
		double CUT_TOT_MET   = CUT[6];
		double CUT_ST		 = CUT[7];
		double CUT_Num_boson = CUT[8];
	
		double eff[NUM], yield[NUM], Err_eff[NUM], Err_yield[NUM];
		double count[NUM]={0.};
		int TagDieOut[NUM]={0};
		//for(int i=0; i<NUM; i++)	count[i]=0.;
		
		//### Loop over different processes
		std::vector<MyEvent> vec;
		for(int k=0; k<NUM; k++){
			if(k==0) vec = vec_sig; 
			else	 vec = vec_bg[k-1]; 
			
			//### Storing particle info & Apply selection cuts
			double weighting = L*pb2fb*X[k]/(double)vec.size();
			for( std::vector<MyEvent>::iterator it=vec.begin(); it!=vec.end(); it++){
				double Num_lep = 0, Num_jet = 0, SPT_lep = 0, SPT_jet = 0, TOT_met = 0, ST = 0;
				std::vector<double> PT_lep;        PT_lep.clear();
				std::vector<double> PT_jet;        PT_jet.clear();
				std::vector<double> MASS_boson;    MASS_boson.clear();
				std::vector<double> PT_chosenJet;  PT_chosenJet.clear();
				std::vector<double> Eta_chosenJet; Eta_chosenJet.clear();
				std::vector<double> Phi_chosenJet; Phi_chosenJet.clear();

				//lep
				for(int i=0; i < it->Electrons.size(); i++){
					if(it->Electrons[i].pt > CUT_PT_lep) Num_lep += 1;
					SPT_lep += it->Electrons[i].pt;
					PT_lep.push_back(it->Electrons[i].pt);
				}
				for(int i=0; i < it->Muons.size(); i++){
					if(it->Muons[i].pt > CUT_PT_lep) Num_lep += 1;
					SPT_lep += it->Muons[i].pt;
					PT_lep.push_back(it->Muons[i].pt);
				}
				//met
				for(int i=0; i < it->METs.size(); i++) TOT_met += it->METs[i].pt;
				//jet
				for(int i=0; i < it->Jets.size(); i++){
					if(it->Jets[i].pt > CUT_PT_jet) Num_jet += 1;
					SPT_jet += it->Jets[i].pt;
					PT_jet.push_back(it->Jets[i].pt);
				}
				//H-boson
				for(int i=0; i < it->H_Bosons.size(); i++){
					MASS_boson.push_back(it->H_Bosons[i].mass);
					PT_chosenJet.push_back(it->Jets[it->H_Bosons[i].jetIds.first].pt);
					PT_chosenJet.push_back(it->Jets[it->H_Bosons[i].jetIds.second].pt);
					Eta_chosenJet.push_back(it->Jets[it->H_Bosons[i].jetIds.first].eta);
					Eta_chosenJet.push_back(it->Jets[it->H_Bosons[i].jetIds.second].eta);
					Phi_chosenJet.push_back(it->Jets[it->H_Bosons[i].jetIds.first].phi);
					Phi_chosenJet.push_back(it->Jets[it->H_Bosons[i].jetIds.second].phi);
				}
				//ST
				ST = SPT_lep + SPT_jet + TOT_met;

				hist_Ori[0][k] -> Fill((it->Electrons.size()+it->Muons.size()),weighting);
    	    	hist_Ori[1][k] -> Fill(it->Jets.size(),weighting);
    	    	hist_Ori[4][k] -> Fill(SPT_lep,weighting);
    	    	hist_Ori[5][k] -> Fill(SPT_jet,weighting);
				hist_Ori[6][k] -> Fill(TOT_met,weighting);
				hist_Ori[7][k] -> Fill(ST,weighting);
				hist_Ori[8][k] -> Fill(it->H_Bosons.size(),weighting);
				for(int i=0; i<PT_lep.size(); i++) hist_Ori[2][k] -> Fill(PT_lep.at(i),weighting);
				for(int i=0; i<PT_jet.size(); i++) hist_Ori[3][k] -> Fill(PT_jet.at(i),weighting);
				for(int i=0; i<MASS_boson.size(); i++)    hist_Ori[9][k]  -> Fill(MASS_boson.at(i),weighting);
				for(int i=0; i<PT_chosenJet.size(); i++)  hist_Ori[10][k] -> Fill(PT_chosenJet.at(i),weighting);
				for(int i=0; i<Eta_chosenJet.size(); i++) hist_Ori[11][k] -> Fill(Eta_chosenJet.at(i),weighting);
				for(int i=0; i<Phi_chosenJet.size(); i++) hist_Ori[12][k] -> Fill(Phi_chosenJet.at(i),weighting);

				if(it->H_Bosons.size() < CUT_Num_boson) continue;
				if(Num_lep < CUT_Num_lep) continue;
				if(SPT_lep < CUT_SPT_lep) continue;
				if(TOT_met < CUT_TOT_MET) continue;
				if(Num_jet < CUT_Num_jet) continue;
				if(SPT_jet < CUT_SPT_jet) continue;
				if(ST < CUT_ST) continue;
				
				if(N==0){//for full cut
					hist_Cut[0][k] -> Fill(Num_lep,weighting);
    	    		hist_Cut[1][k] -> Fill(Num_jet,weighting);
    	    		hist_Cut[4][k] -> Fill(SPT_lep,weighting);
    	    		hist_Cut[5][k] -> Fill(SPT_jet,weighting);
					hist_Cut[6][k] -> Fill(TOT_met,weighting);
					hist_Cut[7][k] -> Fill(ST,weighting);
					hist_Cut[8][k] -> Fill(it->H_Bosons.size(),weighting);
					for(int i=0; i<PT_lep.size(); i++) if(PT_lep.at(i) > CUT_PT_lep) hist_Cut[2][k] -> Fill(PT_lep.at(i),weighting);
					for(int i=0; i<PT_jet.size(); i++) if(PT_jet.at(i) > CUT_PT_jet) hist_Cut[3][k] -> Fill(PT_jet.at(i),weighting);
					for(int i=0; i<MASS_boson.size(); i++)    hist_Cut[9][k]  -> Fill(MASS_boson.at(i),weighting);
					for(int i=0; i<PT_chosenJet.size(); i++)  hist_Cut[10][k] -> Fill(PT_chosenJet.at(i),weighting);
					for(int i=0; i<Eta_chosenJet.size(); i++) hist_Cut[11][k] -> Fill(Eta_chosenJet.at(i),weighting);
					for(int i=0; i<Phi_chosenJet.size(); i++) hist_Cut[12][k] -> Fill(Phi_chosenJet.at(i),weighting);
					count[k] += 1.;
				} else{//for N-1 cuts
					if(N-1==0) hist_Nm1_Plot[0][k] -> Fill(Num_lep,weighting);
					if(N-1==1) hist_Nm1_Plot[1][k] -> Fill(Num_jet,weighting);
					if(N-1==4) hist_Nm1_Plot[4][k] -> Fill(SPT_lep,weighting);
					if(N-1==5) hist_Nm1_Plot[5][k] -> Fill(SPT_jet,weighting);
					if(N-1==6) hist_Nm1_Plot[6][k] -> Fill(TOT_met,weighting);
					if(N-1==7) hist_Nm1_Plot[7][k] -> Fill(ST,weighting);
					if(N-1==8) hist_Nm1_Plot[8][k] -> Fill(it->H_Bosons.size(),weighting);
					if(N-1==2) for(int i=0; i<PT_lep.size(); i++) if(PT_lep.at(i) > CUT_PT_lep) hist_Nm1_Plot[2][k] -> Fill(PT_lep.at(i),weighting);
					if(N-1==3) for(int i=0; i<PT_jet.size(); i++) if(PT_jet.at(i) > CUT_PT_jet) hist_Nm1_Plot[3][k] -> Fill(PT_jet.at(i),weighting);
					if(N-1==9) for(int i=0; i<MASS_boson.size(); i++)    hist_Nm1_Plot[9][k]  -> Fill(MASS_boson.at(i),weighting);
					if(N-1==10) for(int i=0; i<PT_chosenJet.size(); i++)  hist_Nm1_Plot[10][k] -> Fill(PT_chosenJet.at(i),weighting);
					if(N-1==11) for(int i=0; i<Eta_chosenJet.size(); i++) hist_Nm1_Plot[11][k] -> Fill(Eta_chosenJet.at(i),weighting);
					if(N-1==12) for(int i=0; i<Phi_chosenJet.size(); i++) hist_Nm1_Plot[12][k] -> Fill(Phi_chosenJet.at(i),weighting);
				}
			}//end of iterator over events
			if(count[k]==0) {count[k]=1; TagDieOut[k]=1; } // to estimate the upper bound of yield for zero-count processes
			eff[k] = count[k] / (double)vec.size(); 
			Err_eff[k] = sqrt( eff[k]*(1-eff[k])/ (double)vec.size() );
			yield[k] = L*pb2fb*X[k]*eff[k]; 
			Err_yield[k] = sqrt( pow(Err_X[k]/X[k],2) + pow(Err_eff[k]/eff[k],2) )*yield[k];//relative_err times yield
		}//end of k over different processes

		if(N==0){//for full cut
			double yield_sig = yield[0];
			double yield_bg = 0.;
			for(int k=0; k<NUM; k++){
				if(TagDieOut[k]==0){
					std::cout<<setw(10)<<processes[k]<<"\tyield = "<<setw(10)<<yield[k]<<" +/- "<<setw(10)<<Err_yield[k]
							 <<"("<<setw(8)<<(Err_yield[k]/yield[k])*100<<"%)"
							 <<"\tSurvived events = "<<setw(5)<<count[k]<<" +/- "<<setw(10)<<sqrt(count[k]*(1-eff[k]))
							 <<"("<<setw(10)<<eff[k]*sqrt((1-eff[k])/count[k])*100<<"%)"
							 <<"\teff = "<<setw(10)<<eff[k]<<" +/- "<<setw(10)<<Err_eff[k]<<"("<<Err_eff[k]/eff[k]*100<<"%)"<<std::endl;
				} else if(TagDieOut[k]==1){
					std::cout<<setw(10)<<processes[k]<<"\tyield < "<<setw(10)<<yield[k]<<" +/- "<<setw(10)<<Err_yield[k]
							 <<"("<<setw(8)<<(Err_yield[k]/yield[k])*100<<"%)"
							 <<"\tSurvived events < "<<setw(5)<<count[k]<<" +/- "<<setw(10)<<sqrt(count[k]*(1-eff[k]))
							 <<"("<<setw(10)<<eff[k]*sqrt((1-eff[k])/count[k])*100<<"%)"
							 <<"\teff < "<<setw(10)<<eff[k]<<" +/- "<<setw(10)<<Err_eff[k]<<"("<<Err_eff[k]/eff[k]*100<<"%)"<<std::endl;
				}
				if(k!=0) yield_bg += yield[k];
			}
			//double significance = yield_sig / sqrt( yield_sig + yield_bg );
			printf("\nyield_sig = %f, yield_bg = %f\n",yield_sig,yield_bg);
			double significance = sqrt( 2*( (yield_sig + yield_bg)*TMath::Log(1+(yield_sig/yield_bg)) - yield_sig ) );
			std::cout<<"significance_worst = "<<significance<<std::endl;

			//calculate significance with yield = 0 processes
			for(int k=1; k<NUM; k++){
				printf("(%d) TagDieOut = %d\t",k,TagDieOut[k]);
				if(TagDieOut[k]==1) yield_bg -= yield[k];
			}
			printf("\nyield_sig = %f, yield_bg = %f\n",yield_sig,yield_bg);
			significance = sqrt( 2*( (yield_sig + yield_bg)*TMath::Log(1+(yield_sig/yield_bg)) - yield_sig ) );
			std::cout<<"significance_pure = "<<significance<<std::endl;

		} else{//for N-1 cuts
			CutLine[N-1] = new TLine(factors.at(N-1),0,factors.at(N-1),1000);
			CutLine[N-1]->SetLineWidth(3);
			CutLine[N-1]->SetLineColor(CutLineColor);
		}
	}//end of N


	//### Make Stack Plots!!!
	for(int i=0; i<NUM_CUT; i++){
		TPad *pad;
		if(i==0 || i==1) pad = new TPad("pad","",0.38,0.38,0.9,0.898,-1,-1,0);
		else			 pad = new TPad("pad","",0.32,0.32,0.9,0.898,-1,-1,0);
		TLegend *legend_ori = new TLegend(0.74,0.5,0.9,0.9);
		TLegend *legend_cut = new TLegend(0.74,0.5,0.9,0.9);

		//Num_lep, Num_jet, PT_lep, PT_jet, SPT_lep, SPT_jet, TOT_MET, Num_Boson, M_Boson, PT_cjet, Eta_cjet, Phi_cjet
		double binWidth;
		double scaleOri[2][NUM_CUT] = {{5e+4, 1e+4, 5e+3, 5e+4, 100, 8e+3, 100, 1e+4, 5e+4},
								 	   {3e+4, 5e+3, 100, 500, 10, 1e+4, 1e+4, 1e+4, 3e+4}};
		double scaleCut[2][NUM_CUT] = {{1e+4, 60, 40, 400, 100, 100, 100, 200, 200},
								 	   {100, 10, 10, 10, 10, 10, 10, 10, 10}};


		for(int i_proc=0; i_proc<NUM; i_proc++){
			DrawSetting(hist_Ori[i][i_proc], legend_ori, processes[i_proc], histName[i], units[i], "entries", 1001, color[i_proc]);
			DrawSetting(hist_Nm1_Plot[i][i_proc], legend_cut, processes[i_proc], histName[i], units[i], "entries", 1001, color[i_proc]);
			if(i_proc==0){
				hist_Stack_Ori[i]->Add(hist_Ori[i][i_proc]);
				hist_Stack_Copy[i]->Add(hist_Ori[i][i_proc]);
				hist_Stack_Cut[i]->Add(hist_Nm1_Plot[i][i_proc]);
			}
			else{
				hist_Stack_Ori[i]->Add(hist_Ori[i][NUM-i_proc]);
				hist_Stack_Copy[i]->Add(hist_Ori[i][NUM-i_proc]);
				hist_Stack_Cut[i]->Add(hist_Nm1_Plot[i][NUM-i_proc]);
			}
		}
		gPad->SetLogy();
		//if(SIG==2 && i==6) hist_Stack_Ori[i]->SetMinimum(0.0001);
		hist_Stack_Ori[i]->Draw();
		if(( (SIG==1 && !(i==2 || i==3 || i==4 || i==6)) || (SIG==2 && !(i==2 || i==3 || i==4 )) )){//zoom-in plot
			pad->Draw();
			pad->cd();
			gPad->SetLogy();
			hist_Stack_Copy[i]->SetTitle("");
			hist_Stack_Copy[i]->SetMaximum(scaleOri[SIG-1][i]);
			hist_Stack_Copy[i]->Draw();
			can->cd();
		}
		legend_ori->Draw("same");
		can->SaveAs(Form("%s%s/rootfile/OriStack_%s.root",savingPath,suffix,histName[i]));
		can->SaveAs(Form("%s%s/png/OriStack_%s.png",savingPath,suffix,histName[i]));

		//if(SIG==1) gPad->SetLogy(0);
		hist_Stack_Cut[i]->SetMaximum(scaleCut[SIG-1][i]);
		hist_Stack_Cut[i]->Draw();
		legend_cut->Draw("same");
		CutLine[i]->SetY2(scaleCut[SIG-1][i]);
		CutLine[i]->Draw("same");
		can->SaveAs(Form("%s%s/rootfile/CutStack_%s.root",savingPath,suffix,histName[i]));
		can->SaveAs(Form("%s%s/png/CutStack_%s.png",savingPath,suffix,histName[i]));
	}

	hist->Draw();
	can->SaveAs(Form("%s%s/png/AppearedSignificance.png",savingPath,suffix));
	can->SaveAs(Form("%s%s/png/EvolutionTMVA.png",savingPath,suffix));

	fout->Write();
	fout->Close();
	return 1;
}

