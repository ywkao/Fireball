// PID code/particle
// 1/d,2/u,3/s,4/c,5/b,6/t,7/b',8/t'
// 11/e-,12/nu_e,13/mu-,14/nu_mu,15/tau-,16/nu_tau
// 21/g, 22/gamma, 23/z, 24/W+, 25/h, 211/pi+, 2212/proton
// status==1 particles are those that should be processed directly by a detetor simulation
// status==2 particles are unstable, such as a pi0 that will later decay to gamma
// status==3 particles label the beam and the hard process
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
#include <TNtuple.h>
using namespace std;

int main(int argc, char *argv[]){
	//const char* path = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/pptt_test/Events/run_ht_min_1000/tag_1_delphes_events.root";
	const char* input  = argv[1];
	const char* output = argv[2];
	TChain *chain = new TChain("Delphes");
	chain->Add(input);
	ExRootTreeReader *TreeReader = new ExRootTreeReader(chain);
	TClonesArray *branchGenParticle = TreeReader->UseBranch("Particle");
	
	TFile   *file   = new TFile(output,"RECREATE");
	TH1D    *hist   = new TH1D("hist_gen_HT",";HT [GeV];Entries",100,0,3000);
	TNtuple *ntuple = new TNtuple("ntuple_gen_HT",";HT [GeV];Entries","quantity");
	//TCanvas *can  = new TCanvas("can","",800,600);

	for(int entry=0; entry<TreeReader->GetEntries(); entry++){
		TreeReader->ReadEntry(entry);
		double HT=0;
		for(int i=0; i<branchGenParticle->GetEntries(); i++){
			GenParticle *particle = (GenParticle*)branchGenParticle->At(i);
			if((abs(particle->PID)==1 || abs(particle->PID)==2 || abs(particle->PID)==3 || abs(particle->PID)==4 || abs(particle->PID)==5 || abs(particle->PID)==21) && particle->Status==3){
				HT += particle->PT;
			}
		}
		hist->Fill(HT);
		ntuple->Fill(HT);
	}

	//can->SaveAs();
	file->Write();
	file->Close();
}
