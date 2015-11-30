
#include <iostream>
#include <vector>
#include <math.h>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"

using namespace std;
using namespace TMVA;

class MyEvent {
public:
    float mass,pt,eta;
};

std::vector<MyEvent> vec_signal, vec_data;

void importEvents(std::vector<MyEvent>& vec, TString filename)
{
    std::cout << "Loading " << filename << std::endl;
    
    TFile *fin = new TFile(filename);
    TTree *tree = (TTree*)fin->Get("tree");
    
    float _mass,_pt,_eta;
    tree->SetBranchAddress("mass",&_mass);
    tree->SetBranchAddress("pt",&_pt);
    tree->SetBranchAddress("eta",&_eta);
    
    for (int evt=0; evt<tree->GetEntries(); evt++) {
        tree->GetEntry(evt);
        
        if (fabs(_mass-10.)<0.1*6.5 && _pt>8. && fabs(_eta)<2.5) {
            MyEvent event;
            event.mass = _mass;
            event.pt   = _pt;
            event.eta  = _eta;
            
            vec.push_back(event);
        }
    }
    
    std::cout << vec.size() << " events stored." << std::endl;
}

class MyFitness : public IFitterTarget {
public:
    MyFitness() : IFitterTarget() {
    }
    
    Double_t EstimatorFunction( std::vector<Double_t> & factors ){
        
        double CUT_PT  = factors.at(0);
        double CUT_ETA = factors.at(1);
        
        double S = 0., N = 0.;
        
        // loop over signal sample
        for( std::vector<MyEvent>::iterator it=vec_signal.begin(); it!=vec_signal.end(); it++) {
            
            if (fabs(it->mass-10.)>=0.1*3) continue; // 3 sigma mass cut
            if (it->pt<=CUT_PT) continue;
            if (fabs(it->eta)>=CUT_ETA) continue;
            
            S += 1.;
        }
        S *= 0.25; // assuming signal should be scaled by a factor of 0.25 from MC
        
        // loop over data
        for( std::vector<MyEvent>::iterator it=vec_data.begin(); it!=vec_data.end(); it++) {
            
            if (fabs(it->mass-10.)<0.1*3.5 || fabs(it->mass-10.)>=0.1*6.5) continue; // 3.5-6.5 sigma mass cut
            if (it->pt<=CUT_PT) continue;
            if (fabs(it->eta)>=CUT_ETA) continue;
            
            N += 1.;
        }

        double significance = S/sqrt(S+N);
        
        return -significance;
    }
};

int main()
{
    importEvents(vec_signal,"signal.root");
    importEvents(vec_data,"data.root");
    
    std::vector<Interval*> ranges;
    ranges.push_back( new Interval(8.,100.,93) );
    ranges.push_back( new Interval(0.1,2.5,25) );
    
    std::cout << "\nInitial ranges:" << std::endl;
    for( std::vector<Interval*>::iterator it = ranges.begin(); it != ranges.end(); it++ ){
        std::cout << " range: " << (*it)->GetMin() << "   " << (*it)->GetMax() << std::endl;
    }
    
    IFitterTarget* myFitness = new MyFitness();
    
    GeneticAlgorithm mg( *myFitness, 100, ranges );
    
#define CONVSTEPS 20
#define CONVCRIT 0.0001
#define SCSTEPS 10
#define SCRATE 5
#define SCFACTOR 0.95
    
    do {
        mg.Init();
        
        mg.CalculateFitness();
        
        mg.GetGeneticPopulation().Print(0);

        mg.GetGeneticPopulation().TrimPopulation();
        
        mg.SpreadControl( SCSTEPS, SCRATE, SCFACTOR );
        
    } while (!mg.HasConverged( CONVSTEPS, CONVCRIT ));
    
    std::cout << "\nBest factors:" << std::endl;
    
    GeneticGenes* genes = mg.GetGeneticPopulation().GetGenes( 0 );
    std::vector<Double_t> gvec;
    gvec = genes->GetFactors();
    int n = 0;
    for( std::vector<Double_t>::iterator it = gvec.begin(); it<gvec.end(); it++ ){
        std::cout << "FACTOR " << n << " : " << (*it) << std::endl;
        n++;
    }
}

