#ifndef _LEAFDRAW_C_
#define _LEAFDRAW_C_
#include "../header/LeafDraw.h"
const char* DirOutput = "/home/xiaokao/Desktop/MG5_aMC_v2_2_3/Delphes-3.2.0/workspace/tmp";
void LeafDraw(){
    //TFile  *file = new TFile("/home/xiaokao/Desktop/MG5_aMC_v2_2_3/Delphes-3.2.0/workspace/result.root");
    //TTree  *tree = (TTree*) file->Get("mytree");
    
    TChain *chain = new TChain("mytree");
    chain->Add("/home/xiaokao/Desktop/MG5_aMC_v2_2_3/Delphes-3.2.0/workspace/result.root");

    MyGenParticle VB, W, Z, Jet, Lep;
    chain->SetBranchAddress("VB_Multiplicity", &VB.Multiplicity );
    chain->SetBranchAddress("VB_Momentum"    , &VB.Momentum     );
    chain->SetBranchAddress("VB_Energy"      , &VB.Energy       );
    chain->SetBranchAddress("VB_Eta"         , &VB.Eta          );
    chain->SetBranchAddress("VB_Phi"         , &VB.Phi          );
    chain->SetBranchAddress("Z_Multiplicity" , &Z.Multiplicity  );
    chain->SetBranchAddress("Z_Momentum"     , &Z.Momentum      );
    chain->SetBranchAddress("Z_Energy"       , &Z.Energy        );
    chain->SetBranchAddress("Z_Eta"          , &Z.Eta           );
    chain->SetBranchAddress("Z_Phi"          , &Z.Phi           );
    chain->SetBranchAddress("W_Multiplicity" , &W.Multiplicity  );
    chain->SetBranchAddress("W_Momentum"     , &W.Momentum      );
    chain->SetBranchAddress("W_Energy"       , &W.Energy        );
    chain->SetBranchAddress("W_Eta"          , &W.Eta           );
    chain->SetBranchAddress("W_Phi"          , &W.Phi           );
    chain->SetBranchAddress("Lep_Multiplicity" , &Lep.Multiplicity  );
    chain->SetBranchAddress("Lep_Momentum"     , &Lep.Momentum      );
    chain->SetBranchAddress("Lep_Energy"       , &Lep.Energy        );
    chain->SetBranchAddress("Lep_Eta"          , &Lep.Eta           );
    chain->SetBranchAddress("Lep_Phi"          , &Lep.Phi           );
    chain->SetBranchAddress("Jet_Multiplicity" , &Jet.Multiplicity  );
    chain->SetBranchAddress("Jet_Momentum"     , &Jet.Momentum      );
    chain->SetBranchAddress("Jet_Energy"       , &Jet.Energy        );
    chain->SetBranchAddress("Jet_Eta"          , &Jet.Eta           );
    chain->SetBranchAddress("Jet_Phi"          , &Jet.Phi           );

    TH1D *hist_VB_Multiplicity  = new TH1D("hist_VB_Multiplicity"  , "Multiplicity" , 50 , 0.   , 50.   );
    TH1D *hist_VB_Momentum      = new TH1D("hist_VB_Momentum"      , "Momentum"     , 50 , 0.   , 1000. );
    TH1D *hist_VB_Energy        = new TH1D("hist_VB_Energy"        , "Energy"       , 50 , 0.   , 1000. );
    TH1D *hist_VB_Eta           = new TH1D("hist_VB_Eta"           , "Eta"          , 50 , -10. , 10.   );
    TH1D *hist_VB_Phi           = new TH1D("hist_VB_Phi"           , "Phi"          , 50 , -10. , 10.   );
    TH1D *hist_Z_Multiplicity   = new TH1D("hist_Z_Multiplicity"   , "Multiplicity" , 50 , 0.   , 50.   );
    TH1D *hist_Z_Momentum       = new TH1D("hist_Z_Momentum"       , "Momentum"     , 50 , 0.   , 1000. );
    TH1D *hist_Z_Energy         = new TH1D("hist_Z_Energy"         , "Energy"       , 50 , 0.   , 1000. );
    TH1D *hist_Z_Eta            = new TH1D("hist_Z_Eta"            , "Eta"          , 50 , -10. , 10.   );
    TH1D *hist_Z_Phi            = new TH1D("hist_Z_Phi"            , "Phi"          , 50 , -10. , 10.   );
    TH1D *hist_W_Multiplicity   = new TH1D("hist_W_Multiplicity"   , "Multiplicity" , 50 , 0.   , 50.   );
    TH1D *hist_W_Momentum       = new TH1D("hist_W_Momentum"       , "Momentum"     , 50 , 0.   , 1000. );
    TH1D *hist_W_Energy         = new TH1D("hist_W_Energy"         , "Energy"       , 50 , 0.   , 1000. );
    TH1D *hist_W_Eta            = new TH1D("hist_W_Eta"            , "Eta"          , 50 , -10. , 10.   );
    TH1D *hist_W_Phi            = new TH1D("hist_W_Phi"            , "Phi"          , 50 , -10. , 10.   );
    TH1D *hist_Lep_Multiplicity = new TH1D("hist_Lep_Multiplicity" , "Multiplicity" , 50 , 0.   , 50.   );
    TH1D *hist_Lep_Momentum     = new TH1D("hist_Lep_Momentum"     , "Momentum"     , 50 , 0.   , 1000. );
    TH1D *hist_Lep_Energy       = new TH1D("hist_Lep_Energy"       , "Energy"       , 50 , 0.   , 1000. );
    TH1D *hist_Lep_Eta          = new TH1D("hist_Lep_Eta"          , "Eta"          , 50 , -10. , 10.   );
    TH1D *hist_Lep_Phi          = new TH1D("hist_Lep_Phi"          , "Phi"          , 50 , -10. , 10.   );
    TH1D *hist_Jet_Multiplicity = new TH1D("hist_Jet_Multiplicity" , "Multiplicity" , 50 , 0.   , 50.   );
    TH1D *hist_Jet_Momentum     = new TH1D("hist_Jet_Momentum"     , "Momentum"     , 50 , 0.   , 1000. );
    TH1D *hist_Jet_Energy       = new TH1D("hist_Jet_Energy"       , "Energy"       , 50 , 0.   , 1000. );
    TH1D *hist_Jet_Eta          = new TH1D("hist_Jet_Eta"          , "Eta"          , 50 , -10. , 10.   );
    TH1D *hist_Jet_Phi          = new TH1D("hist_Jet_Phi"          , "Phi"          , 50 , -10. , 10.   );

    //tree->Draw("VB_Multiplicity>>hist_VB_Multiplicity","");
    //tree->Draw("VB_Momentum>>hist_VB_Momentum","");
    
    for (int i=0; i<chain->GetEntries(); i++) {
        chain->GetEntry(i);
        FillHist( hist_VB_Multiplicity  , VB.Multiplicity );
        FillHist( hist_VB_Momentum      , VB.Momentum     );
        FillHist( hist_VB_Energy        , VB.Energy       );
        FillHist( hist_VB_Eta           , VB.Eta          );
        FillHist( hist_VB_Phi           , VB.Phi          );
        FillHist( hist_Z_Multiplicity   , Z.Multiplicity );
        FillHist( hist_Z_Momentum       , Z.Momentum     );
        FillHist( hist_Z_Energy         , Z.Energy       );
        FillHist( hist_Z_Eta            , Z.Eta          );
        FillHist( hist_Z_Phi            , Z.Phi          );
        FillHist( hist_W_Multiplicity   , W.Multiplicity );
        FillHist( hist_W_Momentum       , W.Momentum     );
        FillHist( hist_W_Energy         , W.Energy       );
        FillHist( hist_W_Eta            , W.Eta          );
        FillHist( hist_W_Phi            , W.Phi          );
        FillHist( hist_Lep_Multiplicity , Lep.Multiplicity );
        FillHist( hist_Lep_Momentum     , Lep.Momentum     );
        FillHist( hist_Lep_Energy       , Lep.Energy       );
        FillHist( hist_Lep_Eta          , Lep.Eta          );
        FillHist( hist_Lep_Phi          , Lep.Phi          );
        FillHist( hist_Jet_Multiplicity , Jet.Multiplicity );
        FillHist( hist_Jet_Momentum     , Jet.Momentum     );
        FillHist( hist_Jet_Energy       , Jet.Energy       );
        FillHist( hist_Jet_Eta          , Jet.Eta          );
        FillHist( hist_Jet_Phi          , Jet.Phi          );
    }

    DrawHist( hist_VB_Multiplicity );
    DrawHist( hist_VB_Momentum     );
    DrawHist( hist_VB_Energy       );
    DrawHist( hist_VB_Eta          );
    DrawHist( hist_VB_Phi          );
    DrawHist( hist_Z_Multiplicity );
    DrawHist( hist_Z_Momentum     );
    DrawHist( hist_Z_Energy       );
    DrawHist( hist_Z_Eta          );
    DrawHist( hist_Z_Phi          );
    DrawHist( hist_W_Multiplicity );
    DrawHist( hist_W_Momentum     );
    DrawHist( hist_W_Energy       );
    DrawHist( hist_W_Eta          );
    DrawHist( hist_W_Phi          );
    DrawHist( hist_Lep_Multiplicity );
    DrawHist( hist_Lep_Momentum     );
    DrawHist( hist_Lep_Energy       );
    DrawHist( hist_Lep_Eta          );
    DrawHist( hist_Lep_Phi          );
    DrawHist( hist_Jet_Multiplicity );
    DrawHist( hist_Jet_Momentum     );
    DrawHist( hist_Jet_Energy       );
    DrawHist( hist_Jet_Eta          );
    DrawHist( hist_Jet_Phi          );

}
void FillHist(TH1D* hist, double value){
    if(value==-999){}
    else hist->Fill(value);
}
void DrawHist(TH1D* hist){
    hist->Draw();
    c1->SaveAs(Form("%s/%s.png", DirOutput, hist->GetName()));
}
#endif
