#ifndef _LEAFDRAW_H_
#define _LEAFDRAW_H_

#include "TH1D.h"
class MyGenParticle{
public:
    int    Num;
    double Multiplicity;
    double Energy      ;
    double Momentum    ;
    double Momentum2   ;
    double PT          ;
    double Eta         ;
    double Phi         ;

    TH1D *hist_Multiplicity;
    TH1D *hist_Energy      ;
    TH1D *hist_Momentum    ;
    TH1D *hist_Momentum2   ;
    TH1D *hist_PT          ;
    TH1D *hist_Eta         ;
    TH1D *hist_Phi         ;
    
    MyGenParticle(const char* Type){
        hist_Multiplicity  = new TH1D( Form("hist_%s_Multiplicity",Type) , "Multiplicity" , 50 , 0.   , 50.   ); hist_Multiplicity -> Sumw2(); 
        hist_Energy        = new TH1D( Form("hist_%s_Energy"      ,Type) , "Energy"       , 50 , 0.   , 1000. ); hist_Energy       -> Sumw2(); 
        hist_Momentum      = new TH1D( Form("hist_%s_Momentum"    ,Type) , "Momentum"     , 50 , 0.   , 1000. ); hist_Momentum     -> Sumw2(); 
        hist_Momentum2     = new TH1D( Form("hist_%s_Momentum2"   ,Type) , "Momentum2"    , 50 , 0.   , 1000. ); hist_Momentum2    -> Sumw2(); 
        hist_PT            = new TH1D( Form("hist_%s_PT"          ,Type) , "PT"           , 50 , 0.   , 1000. ); hist_PT           -> Sumw2(); 
        hist_Eta           = new TH1D( Form("hist_%s_Eta"         ,Type) , "Eta"          , 50 , -10. , 10.   ); hist_Eta          -> Sumw2(); 
        hist_Phi           = new TH1D( Form("hist_%s_Phi"         ,Type) , "Phi"          , 50 , -10. , 10.   ); hist_Phi          -> Sumw2(); 
    }

};

#endif
