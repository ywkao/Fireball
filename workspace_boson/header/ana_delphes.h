#ifndef _ANA_DELPHES_H_
#define _ANA_DELPHES_H_
#include "TTree.h"
#include "classes/DelphesClasses.h"

class MyGenParticle{
public:
    int    Num;
    double Multiplicity;
    double Momentum    ;
    double Energy      ;
    double Pt          ;
    double Eta         ;
    double Phi         ;
    //MyGenParticle(TTree *&mytree, const char* Type){
    //    mytree->Branch(Form("%s_Multiplicity", Type), &Multiplicity, Form("%s.Multiplicity/D", Type));
    //    mytree->Branch(Form("%s_Momentum"    , Type), &Momentum    , Form("%s.Momentum    /D", Type));
    //    mytree->Branch(Form("%s_Energy"      , Type), &Energy      , Form("%s.Energy      /D", Type));
    //    mytree->Branch(Form("%s_Eta"         , Type), &Eta         , Form("%s.Eta         /D", Type));
    //    mytree->Branch(Form("%s_Phi"         , Type), &Phi         , Form("%s.Phi         /D", Type));
    //}
    void InitCounting(){ Num = 0; }
    void Init(){
        Multiplicity = -999;
        Momentum     = -999;
        Energy       = -999;
        Pt           = -999;
        Eta          = -999;
        Phi          = -999;
    }
    friend void RegisterParticleInfo(MyGenParticle &object, GenParticle *gen){
        object.Num  += 1;
        object.Momentum = sqrt(pow(gen->Px,2) + pow(gen->Py,2) + pow(gen->Pz,2));
        object.Energy   = gen->E;
        object.Pt       = gen->PT;
        object.Eta      = gen->Eta;
        object.Phi      = gen->Phi;
    }
    friend void RegisterMultiplicity(MyGenParticle &object){
        object.Init();
        object.Multiplicity = (double)object.Num;
    }
};

#endif
