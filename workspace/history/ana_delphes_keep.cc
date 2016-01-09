#include <stdio.h>
#include <iostream>
#include <iomanip> //for using setw(n)
#include <vector>
#include <math.h>
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ana_delphes.h"

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

int main(int argc, char *argv[]){
	LOGFile = fopen("INFOLOG","w+"); fclose(LOGFile);// write/update file, to empty the previous content
	message = Form("The chosen bp mass is %d TeV",atoi(argv[1])); INFOLOG(message);
	//std::cout << "The chosen bp mass is " << argv[1] << " TeV" << std::endl;
	SIG = atoi(argv[1]);
	if(SIG==1) savingPath = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace"; 
	if(SIG==2) savingPath = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace"; 

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

	//### X-sec ###//
	if(SIG==1){ X[0] = Xsec.GetCrossSection("fireball_bp_1TeV"); Err_X[0] = Xsec.GetCrossSectionErr("fireball_bp_1TeV");}
	if(SIG==2){ X[0] = Xsec.GetCrossSection("fireball_bp_2TeV"); Err_X[0] = Xsec.GetCrossSectionErr("fireball_bp_2TeV");}
	X[1] = Xsec.GetCrossSection("pptt")   ; Err_X[1] = Xsec.GetCrossSectionErr("pptt")	 ;
	X[2] = Xsec.GetCrossSection("ppvv")   ; Err_X[2] = Xsec.GetCrossSectionErr("ppvv")   ;
	X[3] = Xsec.GetCrossSection("ppttv")  ; Err_X[3] = Xsec.GetCrossSectionErr("ppttv")  ;
	X[4] = Xsec.GetCrossSection("ppvvv")  ; Err_X[4] = Xsec.GetCrossSectionErr("ppvvv")  ;
	X[5] = Xsec.GetCrossSection("ppttvv") ; Err_X[5] = Xsec.GetCrossSectionErr("ppttvv") ;
	X[6] = Xsec.GetCrossSection("pptttt") ; Err_X[6] = Xsec.GetCrossSectionErr("pptttt") ;
	X[7] = Xsec.GetCrossSection("ppvvvv") ; Err_X[7] = Xsec.GetCrossSectionErr("ppvvvv") ;
	X[8] = Xsec.GetCrossSection("ppttvvv"); Err_X[8] = Xsec.GetCrossSectionErr("ppttvvv");

	//### Selection Cuts ###//
	vector<double> factors;
	double list01[9]={2,6,30,40,20,1200,50,1400,0};
	double list02[9]={2,6,30,40,20,1200,50,1400,0};
	std::vector<double>::iterator it;//to store the optimized selection cuts
	if(SIG==1) factors.insert(factors.begin(),list01,list01+9);
	if(SIG==2) factors.insert(factors.begin(),list02,list02+9);
    std::cout << "\nChosen factors:" << std::endl;
	ListSelectionCuts(factors);//print out the best seclectionCuts!

	//### Initialization ###//
	TFile *fout = new TFile(Form("%s/result_weighting_1.root",savingPath),"recreate");
	for(int k=0; k<NUM_hist; k++){//loop over physical quantities
		for(int i=0; i<NUM; i++){//loop over processes
			hist_Ori[k][i] = new TH1D(Form("%s_Ori_%s",processes[i],histName[k]),"",n_bin[SIG-1][k],l_bin[SIG-1][k],h_bin[SIG-1][k]); hist_Ori[k][i] -> Sumw2();
			hist_Cut[k][i] = new TH1D(Form("%s_Cut_%s",processes[i],histName[k]),"",n_bin[SIG-1][k],l_bin[SIG-1][k],h_bin[SIG-1][k]); hist_Cut[k][i] -> Sumw2();
			hist_Nm1[k][i] = new TH1D(Form("%s_N-1_%s",processes[i],histName[k]),"",n_bin[SIG-1][k],l_bin[SIG-1][k],h_bin[SIG-1][k]); hist_Nm1[k][i] -> Sumw2();
		}
		//if(k<NUM_CUT){
		//	if(k==2 || k==3 || k==4 || k==5 || k==6 || k==7){
		//	double binWidth_ori = (h_bin_ori[SIG-1][k]-l_bin_ori[SIG-1][k])/n_bin_ori[SIG-1][k];
		//	double binWidth_cut = (h_bin_cut[SIG-1][k]-l_bin_cut[SIG-1][k])/n_bin_cut[SIG-1][k];
		//	hist_Stack_Ori[k] =new THStack(Form("Ori_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %.0f %s",binWidth_ori,units[k])));
		//	hist_Stack_Copy[k]=new THStack(Form("Copy_Stack_%s",histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %.0f %s",binWidth_ori,units[k])));
		//	hist_Stack_Cut[k] =new THStack(Form("Cut_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %.0f %s",binWidth_cut,units[k])));
		//	}
		//	else{
		//	hist_Stack_Ori[k] = new THStack(Form("Ori_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %s",units[k])));
		//	hist_Stack_Copy[k]= new THStack(Form("Copy_Stack_%s",histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %s",units[k])));
		//	hist_Stack_Cut[k] = new THStack(Form("Cut_Stack_%s" ,histName[k]),Form("%s;%s;%s",histName[k],units[k],Form("entries / %s",units[k])));
		//	}
		//}
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
	
		double eff[NUM], yield[NUM], Err_eff[NUM], Err_yield[NUM];
		double count[NUM]={0.};
		int TagDieOut[NUM]={0};
		
		//### Loop over different processes
		std::vector<MyEvent> vec;
		for(int k=0; k<NUM; k++){
			if(k==0) vec = vec_sig; 
			else	 vec = vec_bg[k-1]; 
			
			//### CHECK PT_CUT ON SPT!!!!!!!!!!!!!!!!!!!!!!!
			//### Storing particle info & Apply selection cuts
			//double weighting = L*pb2fb*X[k]/(double)vec.size();
			double weighting = 1;
			for( std::vector<MyEvent>::iterator it=vec.begin(); it!=vec.end(); it++){
				double Num_lep = 0, Num_jet = 0, SPT_lep = 0, SPT_jet = 0, TOT_MET = 0, ST = 0;
				double SPT_lep_woPTcut = 0, SPT_jet_woPTcut = 0, ST_woPTcut = 0;
				std::vector<double> PT_lep;        PT_lep.clear();
				std::vector<double> PT_jet;        PT_jet.clear();
				std::vector<double> MASS_boson;    MASS_boson.clear();
				std::vector<double> PT_chosenJet;  PT_chosenJet.clear();
				std::vector<double> Eta_chosenJet; Eta_chosenJet.clear();
				std::vector<double> Phi_chosenJet; Phi_chosenJet.clear();

				//lep
				for(int i=0; i < it->Electrons.size(); i++){
					if(it->Electrons[i].pt > CUT_PT_lep){ Num_lep += 1; SPT_lep += it->Electrons[i].pt; }
					SPT_lep_woPTcut += it->Electrons[i].pt;
					PT_lep.push_back(it->Electrons[i].pt);
				}
				for(int i=0; i < it->Muons.size(); i++){
					if(it->Muons[i].pt > CUT_PT_lep){ Num_lep += 1; SPT_lep += it->Muons[i].pt; }
					SPT_lep_woPTcut += it->Muons[i].pt;
					PT_lep.push_back(it->Muons[i].pt);
				}
				//met
				for(int i=0; i < it->METs.size(); i++) TOT_MET += it->METs[i].pt;
				//jet
				for(int i=0; i < it->Jets.size(); i++){
					if(it->Jets[i].pt > CUT_PT_jet){ Num_jet += 1; SPT_jet += it->Jets[i].pt; }
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

				hist_Ori[0][k] -> Fill((it->Electrons.size()+it->Muons.size()),weighting);
    	    	hist_Ori[1][k] -> Fill(it->Jets.size(),weighting);
    	    	hist_Ori[4][k] -> Fill(SPT_lep_woPTcut,weighting);
    	    	hist_Ori[5][k] -> Fill(SPT_jet_woPTcut,weighting);
				hist_Ori[6][k] -> Fill(TOT_MET,weighting);
				hist_Ori[7][k] -> Fill(ST_woPTcut,weighting);
				hist_Ori[8][k] -> Fill(it->H_Bosons.size(),weighting);
				for(int i=0; i<PT_lep.size(); i++) hist_Ori[2][k] -> Fill(PT_lep.at(i),weighting);
				for(int i=0; i<PT_jet.size(); i++) hist_Ori[3][k] -> Fill(PT_jet.at(i),weighting);
				for(int i=0; i<MASS_boson.size(); i++)    hist_Ori[9][k]  -> Fill(MASS_boson.at(i),weighting);
				for(int i=0; i<PT_chosenJet.size(); i++)  hist_Ori[10][k] -> Fill(PT_chosenJet.at(i),weighting);
				for(int i=0; i<Eta_chosenJet.size(); i++) hist_Ori[11][k] -> Fill(Eta_chosenJet.at(i),weighting);
				for(int i=0; i<Phi_chosenJet.size(); i++) hist_Ori[12][k] -> Fill(Phi_chosenJet.at(i),weighting);

				if(it->H_Bosons.size() < CUT_Num_boson) continue;
				if(Num_lep < CUT_Num_lep) continue;
				if(Num_jet < CUT_Num_jet) continue;
				if(SPT_lep < CUT_SPT_lep) continue;
				if(SPT_jet < CUT_SPT_jet) continue;
				if(TOT_MET < CUT_TOT_MET) continue;
				if(ST < CUT_ST) continue;
				
				if(N==0){//for full cut
					hist_Cut[0][k] -> Fill(Num_lep,weighting);
    	    		hist_Cut[1][k] -> Fill(Num_jet,weighting);
    	    		hist_Cut[4][k] -> Fill(SPT_lep,weighting);
    	    		hist_Cut[5][k] -> Fill(SPT_jet,weighting);
					hist_Cut[6][k] -> Fill(TOT_MET,weighting);
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
					if(N-1==0) hist_Nm1[0][k] -> Fill(Num_lep,weighting);
					if(N-1==1) hist_Nm1[1][k] -> Fill(Num_jet,weighting);
					if(N-1==4) hist_Nm1[4][k] -> Fill(SPT_lep,weighting);
					if(N-1==5) hist_Nm1[5][k] -> Fill(SPT_jet,weighting);
					if(N-1==6) hist_Nm1[6][k] -> Fill(TOT_MET,weighting);
					if(N-1==7) hist_Nm1[7][k] -> Fill(ST,weighting);
					if(N-1==8) hist_Nm1[8][k] -> Fill(it->H_Bosons.size(),weighting);
					if(N-1==2) for(int i=0; i<PT_lep.size(); i++) if(PT_lep.at(i) > CUT_PT_lep) hist_Nm1[2][k] -> Fill(PT_lep.at(i),weighting);
					if(N-1==3) for(int i=0; i<PT_jet.size(); i++) if(PT_jet.at(i) > CUT_PT_jet) hist_Nm1[3][k] -> Fill(PT_jet.at(i),weighting);
					if(N-1==9) for(int i=0; i<MASS_boson.size(); i++)    hist_Nm1[9][k]  -> Fill(MASS_boson.at(i),weighting);
					if(N-1==10) for(int i=0; i<PT_chosenJet.size(); i++)  hist_Nm1[10][k] -> Fill(PT_chosenJet.at(i),weighting);
					if(N-1==11) for(int i=0; i<Eta_chosenJet.size(); i++) hist_Nm1[11][k] -> Fill(Eta_chosenJet.at(i),weighting);
					if(N-1==12) for(int i=0; i<Phi_chosenJet.size(); i++) hist_Nm1[12][k] -> Fill(Phi_chosenJet.at(i),weighting);
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
				message = Form("\n[%s]\tNEvents = %d",processes[k], vec.size());
				INFOLOG(message);
				message = Form("Events = %7.0f +/- %10.4f (%4.2f%)", count[k],sqrt(count[k]*(1-eff[k])),eff[k]*sqrt((1-eff[k])/count[k])*100);
				INFOLOG(message);
				message = Form("Yield = %10.4f +/- %10.4f (%4.2f%)",yield[k],Err_yield[k],(Err_yield[k]/yield[k])*100);
				INFOLOG(message);
				message = Form("Eff = %10.4f +/- %10.4f (%4.2f%)", eff[k],Err_eff[k],Err_eff[k]/eff[k]*100);
				INFOLOG(message);
				if(k!=0) yield_bg += yield[k];
			}
			message = Form("\nyield_sig = %f, yield_bg = %f",yield_sig,yield_bg); INFOLOG(message);
			double significance = sqrt( 2*( (yield_sig + yield_bg)*TMath::Log(1+(yield_sig/yield_bg)) - yield_sig ) );
			message = Form("significance_worst = %f", significance); INFOLOG(message);

			//calculate significance with yield = 0 processes
			for(int k=1; k<NUM; k++){
				printf("%s TagDieOut = %d\n",processes[k],TagDieOut[k]);
				if(TagDieOut[k]==1) yield_bg -= yield[k];
			}
			message = Form("\nyield_sig = %f, yield_bg = %f",yield_sig,yield_bg); INFOLOG(message);
			significance = sqrt( 2*( (yield_sig + yield_bg)*TMath::Log(1+(yield_sig/yield_bg)) - yield_sig ) );
			message = Form("significance_pure = %f", significance); INFOLOG(message);

		} else{//for N-1 cuts
			CutLine[N-1] = new TLine(factors.at(N-1),0,factors.at(N-1),1000);
			CutLine[N-1]->SetLineWidth(3);
			CutLine[N-1]->SetLineColor(CutLineColor);
		}
	}//end of N

	fout->Write();
	fout->Close();
	return 1;
}
