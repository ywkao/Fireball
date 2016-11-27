#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "header/ana_delphes.h"
using namespace std;
#define STATUS 3 

int main()
{
    //=== Input Info ===//
    TChain *chain = new TChain("Delphes");
    //chain->Add("/home/xiaokao/Desktop/test/pool/tag_2_delphes_events.root");
    //chain->Add("/home/xiaokao/Desktop/MG5_aMC_v2_2_3/simulation/fireball_testing/Events/run_fireball_2.0TeV/tag_1_delphes_events.root"); //pop out error message
    chain->Add("/home/xiaokao/Desktop/MG5_aMC_v2_2_3/Delphes-3.2.0/workspace/check/tag_1_delphes_events.root");

    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

    TClonesArray *branchParticle  = treeReader->UseBranch("Particle");
    TClonesArray *branchJet       = treeReader->UseBranch("Jet");
    TClonesArray *branchElectron  = treeReader->UseBranch("Electron");
    TClonesArray *branchMuon      = treeReader->UseBranch("Muon");
    TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
    
    //=== To Store ===//
    TFile *fout = new TFile("result.root","recreate");
    TTree *mytree =  new TTree("mytree","mytree");
    MyGenParticle VB(mytree,"VB"), W(mytree,"W"), Z(mytree,"Z"), GenJet(mytree,"GenJet"), GenLep(mytree,"GenLep");//Including Branch Setting!
    //MyParticle    Jet(mytree,"Jet"), Lep(mytree,"Lep");

    int n_entries = treeReader->GetEntries();
    for(int entry = 0; entry < n_entries; entry++) {
        treeReader->ReadEntry(entry);
        if (((entry+1) % 100)==0 || entry+1==n_entries)
            cout<<"processing "<<entry+1<<"/"<<n_entries<<'.'<<'\r'<<flush;
        //=== Gen Particle ===//
        VB.InitCounting(); W.InitCounting(); Z.InitCounting(); GenJet.InitCounting(); GenLep.InitCounting();
        for (int idx=0; idx<branchParticle->GetEntries(); idx++) {
            VB.Init(); W.Init(); Z.Init(); GenJet.Init(); GenLep.Init();
            GenParticle *gen = (GenParticle*) branchParticle->At(idx);
            if (gen->Status!=STATUS) continue;
            if (!(abs(gen->PID)==1  || abs(gen->PID)==2  || abs(gen->PID)==3  || abs(gen->PID)==4  || abs(gen->PID)==5  || abs(gen->PID)==6  ||
                  abs(gen->PID)==11 || abs(gen->PID)==13 || abs(gen->PID)==15 || abs(gen->PID)==21 || abs(gen->PID)==23 || abs(gen->PID)==24 )) continue;
            //=== GenLep GenJet ===/
            if ( abs(gen->PID)==11 || abs(gen->PID)==13 || abs(gen->PID)==15 )
                RegisterParticleInfo(GenLep, gen);
            if ( abs(gen->PID)==1  || abs(gen->PID)==2  || abs(gen->PID)==3  || abs(gen->PID)==4  || abs(gen->PID)==5  || abs(gen->PID)==6 || abs(gen->PID)==21)
                RegisterParticleInfo(GenJet, gen);
            //=== VB Z W bonson ===/
            if(abs(gen->PID)==23 || abs(gen->PID)==24) RegisterParticleInfo(VB, gen);
            if(abs(gen->PID)==23) RegisterParticleInfo(Z, gen);
            if(abs(gen->PID)==24) RegisterParticleInfo(W, gen);
            mytree->Fill();
        }
        RegisterMultiplicity(GenLep);
        RegisterMultiplicity(GenJet);
        RegisterMultiplicity(VB);
        RegisterMultiplicity(Z);
        RegisterMultiplicity(W);
        mytree->Fill();
    }

    fout->Write();
    fout->Close();
}

    //MyGenParticle VB, W, Z, Jet, Lep;
    //mytree->Branch("VB_Multiplicity"  , &VB.Multiplicity  , "VB.Multiplicity/D");
    //mytree->Branch("VB_Energy"        , &VB.Energy        , "VB.Energy      /D");
    //mytree->Branch("VB_Momentum"      , &VB.Momentum      , "VB.Momentum    /D");
    //mytree->Branch("VB_PT"            , &VB.PT            , "VB.PT          /D");
    //mytree->Branch("VB_Eta"           , &VB.Eta           , "VB.Eta         /D");
    //mytree->Branch("VB_Phi"           , &VB.Phi           , "VB.Phi         /D");
    //mytree->Branch("Z_Multiplicity"   , &Z.Multiplicity   , "Z.Multiplicity/D");
    //mytree->Branch("Z_Energy"         , &Z.Energy         , "Z.Energy      /D");
    //mytree->Branch("Z_Momentum"       , &Z.Momentum       , "Z.Momentum    /D");
    //mytree->Branch("Z_PT"             , &Z.PT             , "Z.PT          /D");
    //mytree->Branch("Z_Eta"            , &Z.Eta            , "Z.Eta         /D");
    //mytree->Branch("Z_Phi"            , &Z.Phi            , "Z.Phi         /D");
    //mytree->Branch("W_Multiplicity"   , &W.Multiplicity   , "W.Multiplicity/D");
    //mytree->Branch("W_Energy"         , &W.Energy         , "W.Energy      /D");
    //mytree->Branch("W_Momentum"       , &W.Momentum       , "W.Momentum    /D");
    //mytree->Branch("W_PT"             , &W.PT             , "W.PT          /D");
    //mytree->Branch("W_Eta"            , &W.Eta            , "W.Eta         /D");
    //mytree->Branch("W_Phi"            , &W.Phi            , "W.Phi         /D");
    //mytree->Branch("Lep_Multiplicity" , &Lep.Multiplicity , "Lep.Multiplicity/D");
    //mytree->Branch("Lep_Energy"       , &Lep.Energy       , "Lep.Energy      /D");
    //mytree->Branch("Lep_Momentum"     , &Lep.Momentum     , "Lep.Momentum    /D");
    //mytree->Branch("Lep_PT"           , &Lep.PT           , "Lep.PT          /D");
    //mytree->Branch("Lep_Eta"          , &Lep.Eta          , "Lep.Eta         /D");
    //mytree->Branch("Lep_Phi"          , &Lep.Phi          , "Lep.Phi         /D");
    //mytree->Branch("Jet_Multiplicity" , &Jet.Multiplicity , "Jet.Multiplicity/D");
    //mytree->Branch("Jet_Energy"       , &Jet.Energy       , "Jet.Energy      /D");
    //mytree->Branch("Jet_Momentum"     , &Jet.Momentum     , "Jet.Momentum    /D");
    //mytree->Branch("Jet_PT"           , &Jet.PT           , "Jet.PT          /D");
    //mytree->Branch("Jet_Eta"          , &Jet.Eta          , "Jet.Eta         /D");
    //mytree->Branch("Jet_Phi"          , &Jet.Phi          , "Jet.Phi         /D");
