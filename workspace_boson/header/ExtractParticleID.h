#ifndef _EXTRACTPARTICLEID_H_
#define _EXTRACTPARTICLEID_H_
#include <string>

string ExtractParticleID(const char* histName){
    string particleID;
    if((string)histName=="hist_VB_Multiplicity")  particleID ="Vector boson";
    if((string)histName=="hist_VB_Energy")        particleID ="Vector boson";
    if((string)histName=="hist_VB_Momentum")      particleID ="Vector boson";
    if((string)histName=="hist_VB_Momentum2")     particleID ="Vector boson";
    if((string)histName=="hist_VB_PT")            particleID ="Vector boson";
    if((string)histName=="hist_VB_Eta")           particleID ="Vector boson";
    if((string)histName=="hist_VB_Phi")           particleID ="Vector boson";
    if((string)histName=="hist_Z_Multiplicity")   particleID ="Z boson";
    if((string)histName=="hist_Z_Energy")         particleID ="Z boson";
    if((string)histName=="hist_Z_Momentum")       particleID ="Z boson";
    if((string)histName=="hist_Z_Momentum2")      particleID ="Z boson";
    if((string)histName=="hist_Z_PT")             particleID ="Z boson";
    if((string)histName=="hist_Z_Eta")            particleID ="Z boson";
    if((string)histName=="hist_Z_Phi")            particleID ="Z boson";
    if((string)histName=="hist_W_Multiplicity")   particleID ="W boson";
    if((string)histName=="hist_W_Energy")         particleID ="W boson";
    if((string)histName=="hist_W_Momentum")       particleID ="W boson";
    if((string)histName=="hist_W_Momentum2")      particleID ="W boson";
    if((string)histName=="hist_W_PT")             particleID ="W boson";
    if((string)histName=="hist_W_Eta")            particleID ="W boson";
    if((string)histName=="hist_W_Phi")            particleID ="W boson";
    if((string)histName=="hist_Lep_Multiplicity") particleID ="Lepton";
    if((string)histName=="hist_Lep_Energy")       particleID ="Lepton";
    if((string)histName=="hist_Lep_Momentum")     particleID ="Lepton";
    if((string)histName=="hist_Lep_Momentum2")    particleID ="Lepton";
    if((string)histName=="hist_Lep_PT")           particleID ="Lepton";
    if((string)histName=="hist_Lep_Eta")          particleID ="Lepton";
    if((string)histName=="hist_Lep_Phi")          particleID ="Lepton";
    if((string)histName=="hist_Jet_Multiplicity") particleID ="Jet";
    if((string)histName=="hist_Jet_Energy")       particleID ="Jet";
    if((string)histName=="hist_Jet_Momentum")     particleID ="Jet";
    if((string)histName=="hist_Jet_Momentum2")    particleID ="Jet";
    if((string)histName=="hist_Jet_PT")           particleID ="Jet";
    if((string)histName=="hist_Jet_Eta")          particleID ="Jet";
    if((string)histName=="hist_Jet_Phi")          particleID ="Jet";
    return particleID;
}
#endif
