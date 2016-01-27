#include <iostream>
#include <stdlib.h> //system
#include <TCanvas.h>
#include <TCollection.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <THStack.h>
#include <TStyle.h>
#include "ana_draw_ABCD.h"
using namespace std;

//char *processes[NUM] = {"fireball", "pptt", "ppvv", "ppvtt", "ppvvv", "ppvvtt", "pptttt", "ppvvvv", "ppvvvtt"};
//int color[NUM] = {kRed, kOrange+3, kOrange+1, kBlue, kCyan, kGreen+2, kGreen, kMagenta+2, kMagenta};
//int color[8] = {kOrange+3, kOrange+1, kBlue, kCyan, kGreen+2, kGreen, kMagenta+2, kMagenta};
//importEvents(vec,"/home/xiaokao/Desktop/work_ABCD/skimmed/simulation_delphes_pptt_skimmed.root");

void ana_draw_ABCD(){
	LOGFile = fopen("INFOLOG","w+"); fclose(LOGFile);// write/update file, to empty the previous content
	std::vector<MyProcess> vec;//to store spectrum for each bg process
	importEvents(vec,"ppttvvv");//2
	importEvents(vec,"ppvvvv");
	importEvents(vec,"pptttt");//3
	importEvents(vec,"ppttvv");//1
	importEvents(vec,"ppvvv");
	importEvents(vec,"ppttv");
	importEvents(vec,"ppvv");
	importEvents(vec,"pptt");

	//quickCheck(vec.at(0).HT); //system("PAUSE");
	//quickCheck(vec.at(1).HT); //system("PAUSE");
	//quickCheck(vec.at(2).HT); //system("PAUSE");
	//quickCheck(vec.at(3).HT); //system("PAUSE");
	//quickCheck(vec.at(4).HT); //system("PAUSE");
	//quickCheck(vec.at(5).HT); //system("PAUSE");
	//quickCheck(vec.at(6).HT); //system("PAUSE");
	//quickCheck(vec.at(7).HT);

	MyProcessCombined sum;
	HistCombined(vec,sum);
	
	//TCanvas *can = new TCanvas("can","",800,600);
	//TCanvas *c1 = new TCanvas("c1","ST" ,0,0,650,350);
	//TCanvas *c2 = new TCanvas("c2","HT" ,650,0,650,350);
	//TCanvas *c3 = new TCanvas("c3","LPT",0,380,650,350);
	//TCanvas *c4 = new TCanvas("c4","MET",650,380,650,350);
	TLegend *legend = new TLegend(0.74,0.5,0.9,0.9);

	//DrawHistSpectrum(can,"Expected","ST", vec, sum, legend);
	TCanvas *c1 = new TCanvas("c1","ST" ,0,0,650,350);
	DrawHistSpectrum(c1,"Expected","ST" , vec, sum, legend); legend->Clear();
	TCanvas *c2 = new TCanvas("c2","HT" ,650,0,650,350);
	DrawHistSpectrum(c2,"Expected","HT" , vec, sum, legend); legend->Clear();
	TCanvas *c3 = new TCanvas("c3","LPT",0,380,650,350);
	DrawHistSpectrum(c3,"Expected","LPT", vec, sum, legend); legend->Clear();
	TCanvas *c4 = new TCanvas("c4","MET",650,380,650,350);
	DrawHistSpectrum(c4,"Expected","MET", vec, sum, legend); legend->Clear();
	TCanvas *can = new TCanvas("can","NumJet",800,600);
	DrawHistSpectrum(can,"Expected","NumJet", vec, sum, legend); legend->Clear();

	//DrawStackSpectrum(can,"ST"  , vec, sum);
	//DrawStackSpectrum(c1,"ST"  , vec, sum);
	//DrawStackSpectrum(c2,"HT"  , vec, sum);
	//DrawStackSpectrum(c3,"LPT" , vec, sum);
	//DrawStackSpectrum(c4,"MET" , vec, sum);
	//DrawStackSpectrum(can,"NumJet", vec, sum);

	//DrawHistSpectrum(c1,"ABCD","ST" , sum, legend);
	//DrawHistSpectrum(c2,"ABCD","HT" , sum);
	//DrawHistSpectrum(c3,"ABCD","LPT", sum);
	//DrawHistSpectrum(c4,"ABCD","MET", sum);
	//DrawHistSpectrum(can,"ABCD","NumJet", sum);

	//gStyle->SetOptStat(0);
	//DrawMap(can, sum);
	
	std::cout<<"Completed!"<<std::endl;

}
