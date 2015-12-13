#include <iostream>
#include <stdlib.h> //system
#include <TCanvas.h>
#include <TCollection.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TStyle.h>
#include "ana_draw_ABCD.h"
using namespace std;

//char *processes[NUM] = {"fireball", "pptt", "ppvv", "ppvtt", "ppvvv", "ppvvtt", "pptttt", "ppvvvv", "ppvvvtt"};
//int color[NUM] = {kRed, kOrange+3, kOrange+1, kBlue, kCyan, kGreen+2, kGreen, kMagenta+2, kMagenta};
//int color[8] = {kOrange+3, kOrange+1, kBlue, kCyan, kGreen+2, kGreen, kMagenta+2, kMagenta};
//importEvents(vec,"/home/xiaokao/Desktop/work_ABCD/skimmed/simulation_delphes_pptt_skimmed.root");

void ana_draw_ABCD(){
	std::vector<MyProcess> vec;//to store spectrum for each bg process
	importEvents(vec,"ppttvvv");
	importEvents(vec,"ppvvvv");
	importEvents(vec,"pptttt");
	importEvents(vec,"ppttvv");
	importEvents(vec,"ppvvv");
	importEvents(vec,"ppttv");
	importEvents(vec,"ppvv");
	importEvents(vec,"pptt");

	quickCheck(vec.at(0).HT); //system("PAUSE");
	//quickCheck(vec.at(1).HT); //system("PAUSE");
	//quickCheck(vec.at(2).HT); //system("PAUSE");
	//quickCheck(vec.at(3).HT); //system("PAUSE");
	//quickCheck(vec.at(4).HT); //system("PAUSE");
	//quickCheck(vec.at(5).HT); //system("PAUSE");
	//quickCheck(vec.at(6).HT); //system("PAUSE");
	//quickCheck(vec.at(7).HT);

	MyProcessCombined sum;
	HistCombined(vec,sum);
	
	//TCanvas *c5 = new TCanvas("c5","ST" ,0,0,650,350);
	//TCanvas *c6 = new TCanvas("c6","HT" ,650,0,650,350);
	//TCanvas *c7 = new TCanvas("c7","LPT",0,380,650,350);
	//TCanvas *c8 = new TCanvas("c8","MET",650,380,650,350);
	//TCanvas *can = new TCanvas("can","",800,600);
	
	//DrawStackSpectrum(c5,"ST"  , vec, sum);
	//DrawStackSpectrum(c6,"HT"  , vec, sum);
	//DrawStackSpectrum(c7,"LPT" , vec, sum);
	//DrawStackSpectrum(c8,"MET" , vec, sum);
	//DrawStackSpectrum(can,"NumJet", vec, sum);

	//DrawHistSpectrum(c5,"ST" , sum);
	//DrawHistSpectrum(c6,"HT" , sum);
	//DrawHistSpectrum(c7,"LPT", sum);
	//DrawHistSpectrum(c8,"MET", sum);
	////DrawHistSpectrum(can,"NumJet", sum);

	//DrawMap(can, sum);
	
	std::cout<<"Completed!"<<std::endl;

}
