#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>

const int Nbins=20;
const char *savingPath = "output";

const int NUM_BG = 8;
const double L=100.;
const double pb2fb=1000.;
double X_bg[8] = {815.96, 192.4, 1.77, 0.621, 0.0219, 0.0164, 0.00405, 0.000249}; //pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
double Err_X_bg[8] = {45.51, 4.1, 0.10, 0.073, 0.0013, 9e-04, 9e-05, 3.4e-05}; //pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
char *CutTag[2] = {"woPTcut","PTcut"};
//double X_sig[2] = {0.051, 0.000177}; //1TeV, 2TeV
//double Err_X_sig[2] = {0.001, 4e-06}; //1TeV, 2TeV
//double X[NUM], Err_X[NUM];//to be determined by SIG (the Xsections of all considered processes)
//char *processes[NUM] = {"fireball", "pptt", "ppvv", "ppvtt", "ppvvv", "ppvvtt", "pptttt", "ppvvvv", "ppvvvtt"};
using namespace std;


//#############################################
//### ImportEvents & calculate Baysian Errors
//#############################################
TH2D *histRaw2DMap[NUM_BG][2][2];//processes,PTcut, NjvsNl & HTvsMET
TH1D *histRaw[NUM_BG][2][5][2][2];//processes,PTcut,(ST,HT,MET,LPT), 2*2
TH1D *histYield[NUM_BG][2][5][2][2];//processes,PTcut,(ST,HT,MET,LPT), 2*2
TGraphAsymmErrors *graph[NUM_BG][2][5][2][2];//processes,PTcut,(ST,HT,MET,LPT), 2*2
double *YieldErrors[NUM_BG][2][5][2][2][Nbins][2];//processes,PTcut,(ST,HT,MET,LPT), 2*2, bins, high&low errors
void calculateBayesErrors(int BgID){
	TH1D *hist;
	for(int k=0; k<2; k++){
		//###sum over A B C D regions
		for(int i=0; i<2; i++){
			for(int j=0; j<2; j++){
				if(i==0 && j==0) hist = (TH1D*)histRaw[BgID-1][k][0][i][j]->Clone();
				else hist->Add(Form("Hist.hist%d%d",i,j),1);
			}
		}
		double scale_factor = L*pb2fb*X_bg[BgID-1];
		//###calculate Bayesian Errors
		for(int i=0; i<2; i++){
			for(int j=0; j<2; j++){
				for(int s=0; s<5; s++){//spectrum
					graph[BgID-1][k][s][i][j] = new TGraphAsymmErrors(Nbins);
					graph[BgID-1][k][s][i][j] -> BayesDivide(histRaw[BgID-1][k][s][i][j],hist);
					for(int bin=0; bin<graph[BgID-1][k][s][i][j]->GetN(); bin++){
						double x,y,err;
						graph[BgID-1][k][s][i][j]->GetPoint(bin,x,y);
						histYield[BgID-1][k][s][i][j]->SetBinContent(bin,y*scale_factor);
						err = sqrt( pow(Err_X_bg[BgID-1]/X_bgi[BgID-1],2) + pow(graph[BgID-1][k][s][i][j]->GetErrorYlow(bin)/y,2) );
						YieldErrors[BgID-1][k][s][i][j][bin][0] = err*(y*scale_factor);
						err = sqrt( pow(Err_X_bg[BgID-1]/X_bgi[BgID-1],2) + pow(graph[BgID-1][k][s][i][j]->GetErrorYhigh(bin)/y,2) );
						YieldErrors[BgID-1][k][s][i][j][bin][1] = err*(y*scale_factor);
					}
				}
			}
		}
	}//end of k
}

TH2D *hist2DMap[2][2];//PTcut, NjvsNl & HTvsMET
TH1D *histSpectrum[2][5][2][2];//PTcut,(ST,HT,MET,LPT), 2*2
void histCombined(){
	//###hist2DMap
	for(int k=0; k<2; k++){
		for(int iMap=0; iMap<2; iMap++){
			for(int BgID=0; BgID<NUM_BG; BgID++){
				if(BgID==0)hist2DMap[k][iMap]= (TH2D*)histRaw2DMap[BgID][k][iMap]->Clone();
				hist2DMap[k][iMap]->Add(histRaw2DMap[BgID][k][iMap],1);
			}
		}
	}
	//###histSpectrum
	for(int k=0; k<2; k++){
		for(int i=0; i<2; i++){
			for(int j=0; j<2; j++){
				for(int s=0; s<5; s++){//spectrum
					for(int BgID=0; BgID<NUM_BG; BgID++){
						if(BgID==0) histSpectrum[k][s][i][j] = (TH1D*)histRaw[BgID-1][k][0][i][j]->Clone();
						histSpectrum[k][s][i][j]->Add(histYield[BgID-1][k][s][i][j],1);//PTcut,(ST,HT,MET,LPT), 2*2
					}
				}
			}		
		}
	}
}

void importEvents(int BgID, char *inputFile){
	TFile *fin = new TFile(inputFile,"READ");
	for(int k=0; k<2; k++){
		for(int i=0; i<2; i++){
			for(int j=0; j<2; j++){
				histRaw[BgID-1][k][0][i][j] = (TH1D*)fin->Get(Form("histST_%s%d%d",CutTag[k],i,j));
				histRaw[BgID-1][k][1][i][j] = (TH1D*)fin->Get(Form("histHT_%s%d%d",CutTag[k],i,j));
				histRaw[BgID-1][k][2][i][j] = (TH1D*)fin->Get(Form("histLPT_%s%d%d",CutTag[k],i,j));
				histRaw[BgID-1][k][3][i][j] = (TH1D*)fin->Get(Form("histMET_%s%d%d",CutTag[k],i,j));
				histRaw[BgID-1][k][4][i][j] = (TH1D*)fin->Get(Form("histNumJet_%s%d%d",CutTag[k],i,j));
			}//end of j
		}//end of i
		histRaw2DMap[BgID-1][k][0] = (TH2D*)fin->Get(Form("hist_NjvsNl_%s",CutTag[k]));
		histRaw2DMap[BgID-1][k][1] = (TH2D*)fin->Get(Form("hist_HTvsMET_%s",CutTag[k]));
	}//end of k
	calculateBayesErrors(BgID);
}

//#############################################
//### Miscellaneous Functions
//#############################################
int DecideLocation(int k){//k=2i+j
	if(k==0) return 3; //B(0,0)
	if(k==1) return 1; //A(0,1)
	if(k==2) return 4; //C(1,0)
	if(k==3) return 2; //D(1,1)
}

//class MyHistCollection {
//public:	
//	TH1F *hist00, *hist01, *hist10, *hist11;
//};
//class MyGraphCollection {
//public:	
//	TGraphAsymmErrors *graph00, *graph01, *graph10, *graph11;
//};

//#############################################
//### THE MAIN FUNCTION
//#############################################
void ana_draw_macro(){
	//###import MC events
	//pptt, ppvv, ppvtt, ppvvv, ppvvtt, pptttt, ppvvvv, ppvvvtt
	importEvents(1,"/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/output/skimmed/simulation_delphes_pptt_skimmed.root");
	importEvents(2,"/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/output/skimmed/simulation_delphes_ppvv_skimmed.root");
	importEvents(3,"/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/output/skimmed/simulation_delphes_ppvtt_skimmed.root");
	importEvents(4,"/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/output/skimmed/simulation_delphes_ppvvv_skimmed.root");
	importEvents(5,"/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/output/skimmed/simulation_delphes_ppvvtt_skimmed.root");
	importEvents(6,"/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/output/skimmed/simulation_delphes_pptttt_skimmed.root");
	importEvents(7,"/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/output/skimmed/simulation_delphes_ppvvvv_skimmed.root");
	importEvents(8,"/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/output/skimmed/simulation_delphes_ppvvvtt_skimmed.root");
	std::cout<<"All events have been imported!!!"<<std::endl;
	histCombined();

	TCanvas *can = new TCanvas("can","",1200,600);
	//can->SetBorderSize(1.5);

	//###2d-EntriesPlot: Num_jet vs. Num_lep
	gStyle->SetOptStat(1111); //including over/under flow
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);

	TLine *lineH, *lineV;
	char *MapName[2] = {"NjvsNl", "HTvsMET"};
	double cutH[2]={6,1200}, cutV[2]={2,50}, boundH[2]={7,1200}, boundV[2]={25,5000};
	for(int k=0; k<2; k++){
		for(int iMap=0; iMap<2; iMap++){
			lineH = new TLine(0,cutH[iMap],boundH[iMap],cutH[iMap]);
			lineV = new TLine(cutV[iMap],0,cutV[iMap],boundV[iMap]);
			hist2DMap[k][iMap] -> Draw("colz");
			lineH -> SetLineWidth(2); lineH -> SetLineColor(1); lineH -> Draw("same");
			lineV -> SetLineWidth(2); lineV -> SetLineColor(1); lineV -> Draw("same");
			can -> SaveAs(Form("%s/result_%s_%s.png",savingPath,MapName[iMap],CutTag[k]));
		}
	} can -> Clear();
	
	/*
	//###Spectrum for ABCD regions: ST, HT, MET, LPT; k=0,1,2,3
	//###Calculate Estimated Spectrum
	char *ParaName[4] = {"ST", "HT", "MET", "LPT"};
	double Limit[4] = {2000, 2000, 800, 800};
	TH1D *histSpectrum[4][2][2], *histEstimated[4], *hist;
	//const double constRatio=0.688736;
	for(int k=0; k<4; k++){
		for(int i=0; i<2; i++){
			for(int j=0; j<2; j++) histSpectrum[k][i][j] = (TH1D*)file->Get(Form("hist_%s_%d%d",ParaName[k],i,j));
		}
		histEstimated[k] = new TH1D(Form("hist_%s_estimated",ParaName[k]),"",20,0,Limit[k]); histEstimated[k]->Sumw2(); 
		histEstimated[k]->Divide(histSpectrum[k][0][1],histSpectrum[k][0][0],1,1);
		histEstimated[k]->Multiply(histEstimated[k],histSpectrum[k][1][0],1,1);
		
		//histEstimated[k] = (TH1D*)histSpectrum[k][1][0]->Clone();
		//for(int bin=1; bin<MaxBin+1; bin++)	histEstimated[k]->SetBinContent(bin,constRatio*histEstimated[k]->GetBinContent(bin));
	}

	//###Draw Spectrum
	can->Divide(2,2);
	gStyle->SetOptStat(111111); //including over/under flow

	for(int k=0; k<4; k++){
		for(int i=0; i<2; i++){
			for(int j=0; j<2; j++){
				can->cd(DecideLocation(2*i+j));
				histSpectrum[k][i][j]->SetFillColor(kBlue-9);
				histSpectrum[k][i][j]->Draw("E2");
				//hist = (TH1D*)histEstimated[k]->DrawClone("hist same");
				hist = (TH1D*)histSpectrum[k][i][j]->Clone();
				hist->Draw("hist,same");//h==histogram+errorBar, hist==histogram only
				hist->SetFillStyle(0);
				if(i+j==2){
					int max_bin = histEstimated[k]->GetMaximumBin();
					double max = histEstimated[k]->GetBinContent(max_bin)+histEstimated[k]->GetBinError(max_bin);
					histSpectrum[k][i][j]->SetMaximum(max);

					histEstimated[k]->SetFillColor(kRed-9);
					histEstimated[k]->Draw("E2,same");
					//hist = (TH1D*)histEstimated[k]->DrawClone("hist same");
					hist = (TH1D*)histEstimated[k]->Clone();
					hist->Draw("hist,same");
					hist->SetLineColor(2);
					hist->SetFillStyle(0);
				}
			}//end of j
		}//end of i
		can -> SaveAs(Form("%s/result_%s.png",savingPath,ParaName[k]));
	}//end of k; ST, HT, MET, LPT
	//file->Close();
	*/
}

/*
//-----
	int SIG=1;
	std::vector<double> factors;//to store the optimized selection cuts
	double list[9]={2,6,30,40,20,1200,50,1400,0};
	factors.insert(factors.begin(),list,list+9);
	double CUT_Num_lep   = factors.at(0);
	double CUT_Num_jet   = factors.at(1);
	double CUT_PT_lep    = factors.at(2);
	double CUT_PT_jet    = factors.at(3);
	double CUT_TotLepPT  = factors.at(4);
	double CUT_TotJetPT  = factors.at(5);
	double CUT_MET   	 = factors.at(6);
	double CUT_ST 		 = factors.at(7);
	double CUT_Num_boson = factors.at(8);
//-----
*/

