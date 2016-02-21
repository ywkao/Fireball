#define ana_pythia_1000_cxx
#include "ana_pythia_1000.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

const int STATUS=3;
void ana_pythia_1000::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L ana_pythia_1000.C
//      Root > ana_pythia_1000 t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *file = new TFile(OUTPUT,"RECREATE");
   TH1D *hist = new TH1D("HT",";HT [GeV]; Entries",280,0,7000);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //printf("Event %d \n", Event_Number[0]);
      //printf("Particle: %d \n", Particle_);

      int PID, Status; double HT=0;
      for(int i=0; i<Particle_size; i++){
         int PID = Particle_PID[i], Status = Particle_Status[i];
         if((abs(PID)==1 || abs(PID)==2 || abs(PID)==3 || abs(PID)==4 || abs(PID)==5 || abs(PID)==21) && Status==STATUS) HT+=Particle_PT[i];
         //if((abs(PID)==1 || abs(PID)==2 || abs(PID)==3 || abs(PID)==4 || abs(PID)==5 || abs(PID)==21) && abs(Status)==SetStatus) HT+=Particle_PT[i];
      }
	  hist->Fill(HT);
   }

   file->Write();
   file->Close();
}
