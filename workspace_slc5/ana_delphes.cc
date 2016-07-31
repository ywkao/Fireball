#include <stdio.h>
#include <iostream>
#include <iomanip> //for using setw(n)
#include <vector>
#include <math.h>
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ana_delphes.h"
#include "savingPath.h"
#include "cuts.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TMath.h>
using namespace std;

//const char *savingPath = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace/skimmed";
const char *PREFIX_fireball = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/PROC_4Gen_1/Events/run_";
const char *SUFFIX_tag1 = "tag_1_delphes_events.root";
const char *SUFFIX_tag2 = "tag_2_delphes_events.root";
const char *PREFIX = "/raid1/w/ykao/01source/simulation_delphes_";
const char *SUFFIX;
char *PROCESS;
const char *path;

int main(int argc, char *argv[]){
	//### Load data ###//
	PROCESS = argv[1];
	TAIL = atoi(argv[2]);
	PRESELECTION = atoi(argv[3]);
	if(TAIL==0) SUFFIX = ".root";
	if(TAIL==1) SUFFIX = "_jet_matching.root";
	
	TChain *chain = new TChain("Delphes"); message = Form("[INFO] Importing root files for: %s", PROCESS); INFOLOG(message);
	if((string)PROCESS=="pptt") 	  ChainingEvents("pptt",chain);
	else if((string)PROCESS=="ppvv")  ChainingEvents("ppvv",chain);
	else if((string)PROCESS=="ppvtt") ChainingEvents("ppvtt",chain);
	else if((string)PROCESS=="ppvvv") ChainingEvents("ppvvv",chain);
	else if((string)PROCESS=="fireball_1TeV"  ){path = Form("%s%s/%s",PREFIX_fireball,PROCESS,SUFFIX_tag2); chain->Add(path);}
	else if((string)PROCESS=="fireball_2TeV"  ){path = Form("%s%s/%s",PREFIX_fireball,PROCESS,SUFFIX_tag1); chain->Add(path);}
	else if((string)PROCESS=="fireball_1.1TeV"){path = Form("%s%s/%s",PREFIX_fireball,PROCESS,SUFFIX_tag1); chain->Add(path);}
	else if((string)PROCESS=="fireball_1.2TeV"){path = Form("%s%s/%s",PREFIX_fireball,PROCESS,SUFFIX_tag1); chain->Add(path);}
	else if((string)PROCESS=="fireball_1.3TeV"){path = Form("%s%s/%s",PREFIX_fireball,PROCESS,SUFFIX_tag2); chain->Add(path);}
	else if((string)PROCESS=="fireball_1.4TeV"){path = Form("%s%s/%s",PREFIX_fireball,PROCESS,SUFFIX_tag1); chain->Add(path);}
	else if((string)PROCESS=="fireball_1.5TeV"){path = Form("%s%s/%s",PREFIX_fireball,PROCESS,SUFFIX_tag1); chain->Add(path);}
	else if((string)PROCESS=="fireball_1.6TeV"){path = Form("%s%s/%s",PREFIX_fireball,PROCESS,SUFFIX_tag1); chain->Add(path);}
	else if((string)PROCESS=="fireball_1.7TeV"){path = Form("%s%s/%s",PREFIX_fireball,PROCESS,SUFFIX_tag1); chain->Add(path);}
	else if((string)PROCESS=="fireball_1.8TeV"){path = Form("%s%s/%s",PREFIX_fireball,PROCESS,SUFFIX_tag1); chain->Add(path);}
	else if((string)PROCESS=="fireball_1.9TeV"){path = Form("%s%s/%s",PREFIX_fireball,PROCESS,SUFFIX_tag1); chain->Add(path);}
	else {path = Form("%s%s%s",PREFIX,PROCESS,SUFFIX); chain->Add(path);}
	ExRootTreeReader *TreeReader = new ExRootTreeReader(chain);
	importEvents(vec,TreeReader);

	X 	    = Xsec.GetCrossSection(PROCESS); 
	Err_X   = Xsec.GetCrossSectionErr(PROCESS);
	message = Form("[INFO] Cross-section= %f +/- %f", X, Err_X); INFOLOG(message);

	//### Selection Cuts ###//
	vector<double> factors;
	//double list[9]={2,6,30,40,20,800,50,0,0};
	factors.insert(factors.begin(),list,list+9);

	//### Initialization ###//
	TFile *fout = new TFile(Form("%s/result_%s.root",savingPath,PROCESS),"recreate");
	for(int k=0; k<NUM_hist; k++){//loop over physical quantities
		hist_Ori[k] = new TH1D(Form("Ori_%s",histName[k]),"",n_bin[k],l_bin[k],h_bin[k]); hist_Ori[k] -> Sumw2();
		hist_Cut[k] = new TH1D(Form("Cut_%s",histName[k]),"",n_bin[k],l_bin[k],h_bin[k]); hist_Cut[k] -> Sumw2();
		hist_Nm1[k] = new TH1D(Form("N-1_%s",histName[k]),"",n_bin[k],l_bin[k],h_bin[k]); hist_Nm1[k] -> Sumw2();
		ntuple_Ori[k] = new TNtuple(Form("ntuple_Ori_%s",histName[k]),"","quantity");
		ntuple_Cut[k] = new TNtuple(Form("ntuple_Cut_%s",histName[k]),"","quantity");
		ntuple_Nm1[k] = new TNtuple(Form("ntuple_Nm1_%s",histName[k]),"","quantity");
		if(k==0) ntuple_Ori_HL = new TNtuple("ntuple_Ori_HL","","LPT:HT");
		if(k==0) ntuple_Cut_HL = new TNtuple("ntuple_Cut_HL","","LPT:HT");
	}

	//### Loop over "Full Cut" & "N-1 Cuts" ###//
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
	
		double eff, yield, Err_eff, Err_yield, weighting = 1;//double weighting = L*pb2fb*X[k]/(double)vec.size();
		double count={0.}; int TagDieOut={0};

		//### Storing particle info & Apply selection cuts
		for( std::vector<MyEvent>::iterator it=vec.begin(); it!=vec.end(); it++){

			if(PRESELECTION==1 && !it->PreSelectionPass) continue;

			double Num_lep = 0, Num_jet = 0, SPT_lep = 0, SPT_jet = 0, TOT_MET = 0, ST = 0, GenHT = 0;
			double SPT_lep_woPTcut = 0, SPT_jet_woPTcut = 0, ST_woPTcut = 0;
			std::vector<double> PT_lep;        PT_lep.clear();
			std::vector<double> PT_jet;        PT_jet.clear();
			std::vector<double> MASS_boson;    MASS_boson.clear();
			std::vector<double> PT_chosenJet;  PT_chosenJet.clear();
			std::vector<double> Eta_chosenJet; Eta_chosenJet.clear();
			std::vector<double> Phi_chosenJet; Phi_chosenJet.clear();

			//lep
			for(int i=0; i < it->Electrons.size(); i++){
				if((abs(it->Electrons[i].eta)< ETA) && (it->Electrons[i].pt > CUT_PT_lep)){ Num_lep += 1; SPT_lep += it->Electrons[i].pt; }
				SPT_lep_woPTcut += it->Electrons[i].pt;
				PT_lep.push_back(it->Electrons[i].pt);
			}
			for(int i=0; i < it->Muons.size(); i++){
				if(abs((it->Muons[i].eta) < ETA) && (it->Muons[i].pt > CUT_PT_lep)){ Num_lep += 1; SPT_lep += it->Muons[i].pt; }
				SPT_lep_woPTcut += it->Muons[i].pt;
				PT_lep.push_back(it->Muons[i].pt);
			}
			//met
			for(int i=0; i < it->METs.size(); i++) TOT_MET += it->METs[i].pt;
			//jet
			for(int i=0; i < it->Jets.size(); i++){
				if(abs((it->Jets[i].eta) < ETA) && (it->Jets[i].pt > CUT_PT_jet)){ Num_jet += 1; SPT_jet += it->Jets[i].pt; }
				SPT_jet_woPTcut += it->Jets[i].pt;
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
			ST = SPT_lep + SPT_jet + TOT_MET;
			ST_woPTcut = SPT_lep_woPTcut + SPT_jet_woPTcut + TOT_MET;
			//Gen_HT
			GenHT = it->GenHT;
			//for(int i=0; i < it->GenParticles.size(); i++){
			//	int PID = it->GenParticles[i].PID;
			//	if((abs(PID)==1 || abs(PID)==2 || abs(PID)==3 || abs(PID)==4 || abs(PID)==5 || abs(PID)==21)) GenHT += it->GenParticles[i].pt;
			//}

			if(N==0){//Store the original info once!
				hist_Ori[0] -> Fill((it->Electrons.size()+it->Muons.size()),weighting);
    			hist_Ori[1] -> Fill(it->Jets.size(),weighting);
    			hist_Ori[4] -> Fill(SPT_lep_woPTcut,weighting);
    			hist_Ori[5] -> Fill(SPT_jet_woPTcut,weighting);
				hist_Ori[6] -> Fill(TOT_MET,weighting);
				hist_Ori[7] -> Fill(ST_woPTcut,weighting);
				hist_Ori[8] -> Fill(it->H_Bosons.size(),weighting);
				hist_Ori[13]-> Fill(GenHT,weighting);
				for(int i=0; i<PT_lep.size(); i++) hist_Ori[2] -> Fill(PT_lep.at(i),weighting);
				for(int i=0; i<PT_jet.size(); i++) hist_Ori[3] -> Fill(PT_jet.at(i),weighting);
				for(int i=0; i<MASS_boson.size(); i++)    hist_Ori[9]  -> Fill(MASS_boson.at(i),weighting);
				for(int i=0; i<PT_chosenJet.size(); i++)  hist_Ori[10] -> Fill(PT_chosenJet.at(i),weighting);
				for(int i=0; i<Eta_chosenJet.size(); i++) hist_Ori[11] -> Fill(Eta_chosenJet.at(i),weighting);
				for(int i=0; i<Phi_chosenJet.size(); i++) hist_Ori[12] -> Fill(Phi_chosenJet.at(i),weighting);

				ntuple_Ori[0] -> Fill((double)(it->Electrons.size()+it->Muons.size()));
    			ntuple_Ori[1] -> Fill((double)(it->Jets.size()));
    			ntuple_Ori[4] -> Fill(SPT_lep_woPTcut);
    			ntuple_Ori[5] -> Fill(SPT_jet_woPTcut);
				ntuple_Ori[6] -> Fill(TOT_MET);
				ntuple_Ori[7] -> Fill(ST_woPTcut);
				ntuple_Ori[8] -> Fill((double)(it->H_Bosons.size()));
				ntuple_Ori[13]-> Fill(GenHT);
				ntuple_Ori_HL -> Fill(SPT_lep_woPTcut, SPT_jet_woPTcut);
				for(int i=0; i<PT_lep.size(); i++) ntuple_Ori[2] -> Fill(PT_lep.at(i));
				for(int i=0; i<PT_jet.size(); i++) ntuple_Ori[3] -> Fill(PT_jet.at(i));
				for(int i=0; i<MASS_boson.size(); i++)    ntuple_Ori[9]  -> Fill(MASS_boson.at(i));
				for(int i=0; i<PT_chosenJet.size(); i++)  ntuple_Ori[10] -> Fill(PT_chosenJet.at(i));
				for(int i=0; i<Eta_chosenJet.size(); i++) ntuple_Ori[11] -> Fill(Eta_chosenJet.at(i));
				for(int i=0; i<Phi_chosenJet.size(); i++) ntuple_Ori[12] -> Fill(Phi_chosenJet.at(i));
			}

			if(it->H_Bosons.size() < CUT_Num_boson) continue;
			if(Num_lep < CUT_Num_lep) continue;
			if(Num_jet < CUT_Num_jet) continue;
			if(SPT_lep < CUT_SPT_lep) continue;
			if(SPT_jet < CUT_SPT_jet) continue;
			if(TOT_MET < CUT_TOT_MET) continue;
			if(ST < CUT_ST) continue;
			
			if(N==0){//for full cut
				hist_Cut[0] -> Fill(Num_lep,weighting);
    			hist_Cut[1] -> Fill(Num_jet,weighting);
    			hist_Cut[4] -> Fill(SPT_lep,weighting);
    			hist_Cut[5] -> Fill(SPT_jet,weighting);
				hist_Cut[6] -> Fill(TOT_MET,weighting);
				hist_Cut[7] -> Fill(ST,weighting);
				hist_Cut[8] -> Fill(it->H_Bosons.size(),weighting);
				hist_Cut[13]-> Fill(GenHT,weighting);
				for(int i=0; i<PT_lep.size(); i++) if(PT_lep.at(i) > CUT_PT_lep) hist_Cut[2] -> Fill(PT_lep.at(i),weighting);
				for(int i=0; i<PT_jet.size(); i++) if(PT_jet.at(i) > CUT_PT_jet) hist_Cut[3] -> Fill(PT_jet.at(i),weighting);
				for(int i=0; i<MASS_boson.size(); i++)    hist_Cut[9]  -> Fill(MASS_boson.at(i),weighting);
				for(int i=0; i<PT_chosenJet.size(); i++)  hist_Cut[10] -> Fill(PT_chosenJet.at(i),weighting);
				for(int i=0; i<Eta_chosenJet.size(); i++) hist_Cut[11] -> Fill(Eta_chosenJet.at(i),weighting);
				for(int i=0; i<Phi_chosenJet.size(); i++) hist_Cut[12] -> Fill(Phi_chosenJet.at(i),weighting);

				ntuple_Cut[0] -> Fill(Num_lep);
    			ntuple_Cut[1] -> Fill(Num_jet);
    			ntuple_Cut[4] -> Fill(SPT_lep);
    			ntuple_Cut[5] -> Fill(SPT_jet);
				ntuple_Cut[6] -> Fill(TOT_MET);
				ntuple_Cut[7] -> Fill(ST);
				ntuple_Cut[8] -> Fill((double)it->H_Bosons.size());
				ntuple_Cut[13]-> Fill(GenHT);
				ntuple_Cut_HL -> Fill(SPT_lep, SPT_jet);
				for(int i=0; i<PT_lep.size(); i++) if(PT_lep.at(i) > CUT_PT_lep) ntuple_Cut[2] -> Fill(PT_lep.at(i));
				for(int i=0; i<PT_jet.size(); i++) if(PT_jet.at(i) > CUT_PT_jet) ntuple_Cut[3] -> Fill(PT_jet.at(i));
				for(int i=0; i<MASS_boson.size(); i++)    ntuple_Cut[9]  -> Fill(MASS_boson.at(i));
				for(int i=0; i<PT_chosenJet.size(); i++)  ntuple_Cut[10] -> Fill(PT_chosenJet.at(i));
				for(int i=0; i<Eta_chosenJet.size(); i++) ntuple_Cut[11] -> Fill(Eta_chosenJet.at(i));
				for(int i=0; i<Phi_chosenJet.size(); i++) ntuple_Cut[12] -> Fill(Phi_chosenJet.at(i));

				count += 1.;
			} else{//for N-1 cuts
				if(N-1==0) hist_Nm1[0] -> Fill(Num_lep,weighting);
				if(N-1==1) hist_Nm1[1] -> Fill(Num_jet,weighting);
				if(N-1==4) hist_Nm1[4] -> Fill(SPT_lep,weighting);
				if(N-1==5) hist_Nm1[5] -> Fill(SPT_jet,weighting);
				if(N-1==6) hist_Nm1[6] -> Fill(TOT_MET,weighting);
				if(N-1==7) hist_Nm1[7] -> Fill(ST,weighting);
				if(N-1==8) hist_Nm1[8] -> Fill(it->H_Bosons.size(),weighting);
				if(N-1==2) for(int i=0; i<PT_lep.size(); i++) if(PT_lep.at(i) > CUT_PT_lep) hist_Nm1[2] -> Fill(PT_lep.at(i),weighting);
				if(N-1==3) for(int i=0; i<PT_jet.size(); i++) if(PT_jet.at(i) > CUT_PT_jet) hist_Nm1[3] -> Fill(PT_jet.at(i),weighting);
				if(N-1==9) for(int i=0; i<MASS_boson.size(); i++)    hist_Nm1[9]  -> Fill(MASS_boson.at(i),weighting);
				if(N-1==10) for(int i=0; i<PT_chosenJet.size(); i++)  hist_Nm1[10] -> Fill(PT_chosenJet.at(i),weighting);
				if(N-1==11) for(int i=0; i<Eta_chosenJet.size(); i++) hist_Nm1[11] -> Fill(Eta_chosenJet.at(i),weighting);
				if(N-1==12) for(int i=0; i<Phi_chosenJet.size(); i++) hist_Nm1[12] -> Fill(Phi_chosenJet.at(i),weighting);

				if(N-1==0) ntuple_Nm1[0] -> Fill(Num_lep);
				if(N-1==1) ntuple_Nm1[1] -> Fill(Num_jet);
				if(N-1==4) ntuple_Nm1[4] -> Fill(SPT_lep);
				if(N-1==5) ntuple_Nm1[5] -> Fill(SPT_jet);
				if(N-1==6) ntuple_Nm1[6] -> Fill(TOT_MET);
				if(N-1==7) ntuple_Nm1[7] -> Fill(ST);
				if(N-1==8) ntuple_Nm1[8] -> Fill((double)it->H_Bosons.size());
				if(N-1==2) for(int i=0; i<PT_lep.size(); i++) if(PT_lep.at(i) > CUT_PT_lep) ntuple_Nm1[2] -> Fill(PT_lep.at(i));
				if(N-1==3) for(int i=0; i<PT_jet.size(); i++) if(PT_jet.at(i) > CUT_PT_jet) ntuple_Nm1[3] -> Fill(PT_jet.at(i));
				if(N-1==9) for(int i=0; i<MASS_boson.size(); i++)    ntuple_Nm1[9]  -> Fill(MASS_boson.at(i));
				if(N-1==10) for(int i=0; i<PT_chosenJet.size(); i++)  ntuple_Nm1[10] -> Fill(PT_chosenJet.at(i));
				if(N-1==11) for(int i=0; i<Eta_chosenJet.size(); i++) ntuple_Nm1[11] -> Fill(Eta_chosenJet.at(i));
				if(N-1==12) for(int i=0; i<Phi_chosenJet.size(); i++) ntuple_Nm1[12] -> Fill(Phi_chosenJet.at(i));
			}
		}//end of iterator over events

		if(N==0){//for full cut
			if(count==0) {count=1; TagDieOut=1; } // to estimate the upper bound of yield for zero-count processes
			eff       = count / (double)vec.size(); 
			Err_eff   = sqrt( eff*(1-eff)/ (double)vec.size() );
			yield     = L*pb2fb*X*eff; 
			Err_yield = sqrt( pow(Err_X/X,2) + pow(Err_eff/eff,2) )*yield;//relative_err times yield

			message = Form("\nNEvents = %d; Factor = %f", vec.size(), L*pb2fb*X/(double)vec.size());
			INFOLOG(message);
			message = Form("Events = %7.0f +/- %10.4f (%4.2f%)", count,sqrt(count*(1-eff)),eff*sqrt((1-eff)/count)*100);
			INFOLOG(message);
			message = Form("Yield = %10.4f +/- %10.4f (%4.2f%)",yield,Err_yield,(Err_yield/yield)*100);
			INFOLOG(message);
			message = Form("Eff = %10.4f +/- %10.4f (%4.2f%)", eff,Err_eff,Err_eff/eff*100);
			INFOLOG(message);
			message = Form("TagDieOut = %d\n",TagDieOut);
			INFOLOG(message);
		}
	}//end of N

	fout->Write();
	fout->Close();
	return 1;
}
