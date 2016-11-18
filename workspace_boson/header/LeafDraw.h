#ifndef _LEAFDRAW_H_
#define _LEAFDRAW_H_

class MyGenParticle{
public:
    int    Num;
    double Multiplicity;
    double Momentum    ;
    double Energy      ;
    double Pt          ;
    double Eta         ;
    double Phi         ;
    void InitCounting(){ Num = 0; }
    void Init(){
        Multiplicity = -999;
        Momentum     = -999;
        Energy       = -999;
        Pt           = -999;
        Eta          = -999;
        Phi          = -999;
    }
};

#endif
