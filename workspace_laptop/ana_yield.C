#include <TCanvas.h>
#include <TH1D.h>
#include <TMath.h>
#include <TPad.h>
#include <Math/DistFuncMathCore.h>

#include <stdio.h> //FILE
#include <algorithm>
#include <iterator>
#include <vector>
#include <iomanip>
#include <iostream>

#include "ana_uncertainty.h"
#include "ana_yield.h"
#include "ana_table.h"
#include "FileGenerator.h"

using namespace std;
const int SHOWTABLE_bkg    = 0;
const int SHOWTABLE_sig    = 0;
const int SHOWTABLE_nm1    = 0;
const int SHOWTEST         = 1;
const int SHOWSIGNIFICANCE = 1;
const int MCESTIMATION     = 0;
const int CREATEFILE       = 0;

void ana_yield(){
	importHist(vec_bkg,"pptt");
	importHist(vec_bkg,"ppvv");
	importHist(vec_bkg,"ppvtt");
	importHist(vec_bkg,"ppvvv");
	importHist(vec_bkg,"ppvvtt");
	importHist(vec_bkg,"pptttt");
	importHist(vec_bkg,"ppvvvv");
	importHist(vec_bkg,"ppvvvtt");

	importHist(vec,"fireball_1TeV");
	importHist(vec,"fireball_1.1TeV");
	importHist(vec,"fireball_1.2TeV");
	importHist(vec,"fireball_1.3TeV");
	importHist(vec,"fireball_1.4TeV");
	importHist(vec,"fireball_1.5TeV");
	importHist(vec,"fireball_1.6TeV");
	importHist(vec,"fireball_1.7TeV");
	importHist(vec,"fireball_1.8TeV");
	importHist(vec,"fireball_1.9TeV");
	importHist(vec,"fireball_2TeV");

    if(SHOWTABLE_bkg==1) MakeTable_bkg(vec_bkg);
    if(SHOWTABLE_sig==1) MakeTable_sig(vec);
    if(SHOWTABLE_nm1==1) MakeTable_nm1(vec);

    CombineBackgroundProcess(vec_bkg);

    //#################################//
    //### MC method for uncertainty ###//
    //#################################//
    if(MCESTIMATION==1){
        TCanvas *can = new TCanvas("can","",800,600); //can->SetLogx();
        TH1D *hist = new TH1D("hist","",100,-0.5,0.8);
        f.mean_lumi = L*pb2fb;
        f.rms_lumi  = f.mean_lumi*Err_L;
	    for(unsigned int i=0; i<vec_bkg.size(); i++){
            f.N_b[i]       = vec_bkg.at(i).Entries;
            f.eff_b[i]     = vec_bkg.at(i).Eff[2][0];
            f.mean_xsec[i] = vec_bkg.at(i).X;
            f.rms_xsec[i]  = vec_bkg.at(i).Err_X;
        }
	    for(std::vector<MyProcess>::iterator it=vec.begin(); it!=vec.end(); it++){
            f.mean_N = TotalYield;
            f.N_s    = it->Yield[2][0];
            f.eff_s  = it->Eff[2][0];
            std::pair<double, double> Answer = ana_uncertainty(f,hist);
        }
        hist->Draw();
    }

    if(SHOWTEST==1){
        //conversion between p-value and significance Z
        for(int i=0; i<6; i++){
            double A = (double)i;
            double B = ROOT::Math::normal_cdf_c(A);
            //printf("%.0f, ", A);
            //printf("%5.4e, ", B);
            printf("Z = %d, p-value = %5.4e\n", i, B);
        } printf("\n");
        
        //significance Z with manual observed events number N
        for(int i=0; i<11; i++){
            N = 150. + (double)i*15.;
            NB = 148.;
	        significance = sqrt( 2*( N*TMath::Log(N/NB) - (N-NB) ));
            p_value = ROOT::Math::normal_cdf_c(significance,1.0,0); //x, sigma, x0
            printf("N = %.0f, NB = %0.f, significance = %6.4f, p-value = %5.4f\n", N, NB, significance, p_value);
        } printf("\n");
    }

    //################################//
	//### List the statistics info ###//
    //################################//
    int LoopCounter=0;
	for(std::vector<MyProcess>::iterator it=vec.begin(); it!=vec.end(); it++){
        LoopCounter += 1;
        //=========================//
        //=== X-sec Measurement ===//
        //=========================//
        TotalYield                = it->Yield[2][0] + yield_bg;
        Err_ExpectedSigYield      = sqrt(TotalYield+pow(Err_yield_bg,2));
        MeasuredCrossSection      = it->Yield[2][0]/(it->Eff[2][0]*L*pb2fb);

        //### Intuitive scale up/down ###//
        MeasuredCrossSection_up   = (it->Yield[2][0]+Err_ExpectedSigYield)/((it->Eff[2][0]-it->Err_Eff[2][0])*(L-L*Err_L)*pb2fb);//!?
        MeasuredCrossSection_down = (it->Yield[2][0]-Err_ExpectedSigYield)/((it->Eff[2][0]+it->Err_Eff[2][0])*(L+L*Err_L)*pb2fb);//!?
        Err_MeasuredCrossSection  =pow( Err_ExpectedSigYield / it->Yield[2][0],2) + pow(it->Err_Eff[2][0]/it->Eff[2][0],2) + pow(Err_L,2);//relative error^2
        Err_MeasuredCrossSection  = sqrt(Err_MeasuredCrossSection)*MeasuredCrossSection;//error of measured X-sec

        //### \sigma = N/(L*eff) - (N_B/L)/eff_sig ###//
        First      = TotalYield / (L*pb2fb*it->Eff[2][0]);
        Err_First  = First*sqrt((1/TotalYield) + pow(Err_L,2) + pow(it->Err_Eff[2][0]/it->Eff[2][0],2));
        Second     = yield_bg/(L*pb2fb*it->Eff[2][0]);
        Err_Second = Second*sqrt(pow(Err_yield_bg/yield_bg,2) - pow(Err_L,2) + pow(it->Err_Eff[2][0]/it->Eff[2][0],2));

        //### Summary ###//
        //printf("%8d & $%1.f\\%%$ &    -     & $%.2f\\%%$ & $%3.2e \\pm %2.1e + %2.1e - %2.1e$ & $%3.2e \\pm %2.1e$\\\\\n", 
        //        LoopCounter*100+900, Err_L*100, it->Err_Eff[2][0]/it->Eff[2][0]*100, MeasuredCrossSection, Err_MeasuredCrossSection,
        //        MeasuredCrossSection_up-MeasuredCrossSection, MeasuredCrossSection-MeasuredCrossSection_down,
        //        First-Second, sqrt(pow(Err_First,2)+pow(Err_Second,2)) );
        //printf("N = %.0f \u00B1 %.0f, N_B = %.0f \u00B1 %.0f, lumi = %.0f \u00B1 %.0f, eff = %.3f \u00B1 %2.1e, xsec = %3.2e \u00B1 %3.2e\n", 
        //printf("%8d: %3.2e \u00B1 %3.2e\n", 
        //        LoopCounter*100+900, it->Yield[2][0], Err_ExpectedSigYield);
        //printf("%8d & $%.0f \\pm %.0f$ & $%.0f \\pm %.0f$ & $%.0f \\pm %.0f$ & $%.4f \\pm %.4f$ & $%3.2e \\pm %3.2e$\\\\\n", 
        //        LoopCounter*100+900,
        //        TotalYield, Err_ExpectedSigYield, yield_bg, Err_yield_bg, L, L*Err_L, 
        //        it->Eff[2][0], it->Err_Eff[2][0], MeasuredCrossSection, sqrt(pow(Err_First,2)+pow(Err_Second,2)));
        ////systematic uncertainty!!
        //printf("%3.2e, ",sqrt(pow(Err_First,2)+pow(Err_Second,2)));
        
        
                
        //=====================//
        //=== full cut info ===//
        //=====================//
        if(SHOWSIGNIFICANCE==1){
            yield_sig = it->Yield[2][0];    
            N = yield_sig + yield_bg;
	        significance = sqrt( 2*( N*TMath::Log(N/yield_bg) - (N-yield_bg) ));
            p_value = ROOT::Math::normal_cdf_c(significance,1.0,0); //x, sigma, x0
            //printf("%3.2e, ", significance);
            //printf("%3.2e, ", p_value);
            printf("%16s: significance = %6.4f, p-value = %5.4f (yield_sig = %.3f, yield_bg = %.3f)\n", it->Name, significance, p_value, yield_sig, yield_bg);
        }
        
        
        //============================//
        //=== Files for Limit Plot ===//
        //============================//
        if(CREATEFILE==1){
            Begin(it->Name);
            INFOLOG(it->Name,Form("rate           %5.2f     %5.2f     %5.2f    %5.2f     %5.2f",it->Yield[2][0], ybg1, ybg2, ybg3, ybg4));
            INFOLOG(it->Name,Form("------------"));
            INFOLOG(it->Name,Form("luminosity     lnN %3.2f  %3.2f  %3.2f  %3.2f  %3.2f ", Err_L+1, Err_L+1, Err_L+1, Err_L+1, Err_L+1));
            INFOLOG(it->Name,Form("X_pptt         lnN -  %5.4f -        -        -         ",1+vec_bkg.at(0).Err_X/vec_bkg.at(0).X));
            INFOLOG(it->Name,Form("X_ppvv         lnN -        -  %5.4f -        -         ",1+vec_bkg.at(1).Err_X/vec_bkg.at(1).X));
            INFOLOG(it->Name,Form("X_ppvtt        lnN -        -        -  %5.4f -         ",1+vec_combine.at(0).Err_X/vec_combine.at(0).X));
            INFOLOG(it->Name,Form("X_ppvvv        lnN -        -        -        -  %5.4f  ",1+vec_combine.at(1).Err_X/vec_combine.at(1).X));
            INFOLOG(it->Name,Form("Eff_fireball   lnN %6.4f    -        -        -       - ",1+vec.at(0).Err_Yield[2][0]/vec.at(0).Yield[2][0]));
            INFOLOG(it->Name,Form("Eff_pptt       lnN -  %5.4f -        -        -         ",1+vec_bkg.at(0).Err_Eff[2][0]/vec_bkg.at(0).Eff[2][0]));
            INFOLOG(it->Name,Form("Eff_ppvv       lnN -        -  %5.4f -        -         ",1+vec_bkg.at(1).Err_Eff[2][0]/vec_bkg.at(1).Eff[2][0]));
            INFOLOG(it->Name,Form("Eff_ppvtt      lnN -        -        -  %5.4f -         ",1+vec_combine.at(0).Err_Eff[2][0]/vec_combine.at(0).Eff[2][0]));
            INFOLOG(it->Name,Form("Eff_ppvvv      lnN -        -        -        -  %5.4f  ",1+vec_combine.at(1).Err_Eff[2][0]/vec_combine.at(1).Eff[2][0]));
            INFOLOG(it->Name,"------------\n");
        }
    }
}
