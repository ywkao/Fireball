#include <TCanvas.h>
#include <THStack.h>
#include <TMath.h>
#include <TPad.h>

#include <algorithm>
#include <iterator>
#include <vector>
#include <iomanip>

#include "ana_draw.h"
//#include "ana_data.h"

using namespace std;
const double CUT[NUM] = {2,6,30,40,20,1200,50,1400,0};
const int RebinFactor = 8;

void ana_draw(){

    SetLegend(legend_MC_sig, 1, 43, 16, 0 ,0, 0);
    SetLegend(legend_MC_bkg, 2, 43, 16, 0 ,0, 0);
	//if(SIG==1) importHist("fireball_bp_1TeV");
	//if(SIG==2) importHist("fireball_bp_2TeV");
	importHist("fireball_2TeV");
	importHist("pptt");
	importHist("ppvv");
	importHist("ppvtt");
	importHist("ppvvv");
	importHist("ppvvtt");
	importHist("pptttt");
	importHist("ppvvvv");
	importHist("ppvvvtt");
	importHist("fireball_1TeV");

	int Nbin; double binWidth, binMin, binMax;
	for(int k=0; k<NUM; k++){
		if(!(k==0 || k==1 || k==8)){
            vec.at(0).hist_Ori[k]->Rebin(RebinFactor);
            vec.at(0).hist_Nm1[k]->Rebin(RebinFactor);
        }
		Nbin   = vec.at(0).hist_Ori[k]->GetXaxis()->GetNbins();
		binMin = vec.at(0).hist_Ori[k]->GetXaxis()->GetXmin();
		binMax = vec.at(0).hist_Ori[k]->GetXaxis()->GetXmax();
		binWidth = (binMax-binMin)/(double)Nbin;
		if(k==2 || k==3 || k==4 || k==5 || k==6 || k==7){
		stack_ori[k] =new THStack(Form("Ori_Stack_%s" ,histName[k]),Form(";%s [%s];%s",histName[k],units[k],Form("Entries / %.0f %s",binWidth,units[k])));
		stack_nm1[k] =new THStack(Form("Nm1_Stack_%s" ,histName[k]),Form(";%s [%s];%s",histName[k],units[k],Form("Entries / %.0f %s",binWidth,units[k])));
		} else{
		stack_ori[k] =new THStack(Form("Ori_Stack_%s" ,histName[k]),Form(";%s;%s",units[k],"Entries"));
		stack_nm1[k] =new THStack(Form("Nm1_Stack_%s" ,histName[k]),Form(";%s;%s",units[k],"Entries"));
		}
	}


	//### combine process & stack background! ###//
	for(int k=0; k<NUM; k++){
		for(unsigned int i=0; i<vec.size(); i++){
            if(i==0) continue;//fireball process has been rebinned
			if(!(k==0 || k==1 || k==8)){
			    vec.at(i).hist_Ori[k]->Rebin(RebinFactor);
                vec.at(i).hist_Nm1[k]->Rebin(RebinFactor);
            }

            //for plotting bkg errors
            if(i==1){
                hist_err_ori[k] = (TH1D*)vec.at(i).hist_Ori[k]->Clone();
                hist_err_nm1[k] = (TH1D*)vec.at(i).hist_Nm1[k]->Clone();
                hist_err_ori[k]->SetFillStyle(3001);
                hist_err_nm1[k]->SetFillStyle(3001);
                hist_err_ori[k]->SetFillColor(kGray);
                hist_err_nm1[k]->SetFillColor(kGray);
                hist_err_ori[k]->SetLineColor(kGray);
                hist_err_nm1[k]->SetLineColor(kGray);
                //if(k==0) legend_MC_bkg->AddEntry(hist_err_nm1[0],"Bkg. Errors","f");
            }else if(i!=vec.size()-1){
                hist_err_ori[k] -> Add(vec.at(i).hist_Ori[k]);
                hist_err_nm1[k] -> Add(vec.at(i).hist_Nm1[k]);
            }else{
                continue;
            }
            

            //combine bkg processes
            if(i==1||i==2) continue;//keep pp>tt, pp>vv not mix with other processes
            else if(i==4){//pp>vvv
                hist_ori_ppvvv[k] = (TH1D*)vec.at(i).hist_Ori[k]->Clone();
                hist_nm1_ppvvv[k] = (TH1D*)vec.at(i).hist_Nm1[k]->Clone();
            }else if(i==7){//pp>vvvv
                hist_ori_ppvvv[k] -> Add(vec.at(i).hist_Ori[k]);
                hist_nm1_ppvvv[k] -> Add(vec.at(i).hist_Nm1[k]);
            }else if(i==3){//pp>ttv
                hist_ori_ppttv[k] = (TH1D*)vec.at(i).hist_Ori[k]->Clone();
                hist_nm1_ppttv[k] = (TH1D*)vec.at(i).hist_Nm1[k]->Clone();
            }else{//pp>ttvv, pp>ttvvv, pp>tttt
                hist_ori_ppttv[k] -> Add(vec.at(i).hist_Ori[k]);
                hist_nm1_ppttv[k] -> Add(vec.at(i).hist_Nm1[k]);
            }
        }

        //stack bkg processes
		for(int i=0; i<5; i++){
			if(i+1==5) continue;//this corresponds to signal which should be separated from bkg
            if(i==0){
			    stack_ori[k]->Add(hist_ori_ppvvv[k]);
			    stack_nm1[k]->Add(hist_nm1_ppvvv[k]);
            }else if(i==1){
			    stack_ori[k]->Add(vec.at(2).hist_Ori[k]);
			    stack_nm1[k]->Add(vec.at(2).hist_Nm1[k]);
            }else if(i==2){
			    stack_ori[k]->Add(hist_ori_ppttv[k]);
			    stack_nm1[k]->Add(hist_nm1_ppttv[k]);
            }else{//pp>vv, pp>tt
			    stack_ori[k]->Add(vec.at(5-i-1).hist_Ori[k]);
			    stack_nm1[k]->Add(vec.at(5-i-1).hist_Nm1[k]);
            }
		}
	}


	//### List the statistics info ###//
	double yield_sig = 0, yield_bg=0, Err_yield_bg=0, significance;
    double sur_nm1[NUM]={0}, Err_sur_nm1[NUM]={0};//for each N-1 cut, store the tot survived events 
    double tot_nm1[NUM]={0}, Err_tot_nm1[NUM]={0};//for each N-1 cut, store the tot generated events
    double eff_nm1[NUM]={0}, Err_eff_nm1[NUM]={0};//store the cut eff for each N-1 cut
    double yield_nm1[NUM]={0}, Err_yield_nm1[NUM]={0};//store the cut yield for each N-1 cut
	for(std::vector<MyProcess>::iterator it=vec.begin(); it!=vec.end(); it++){
		if(it==vec.begin()||it+1==vec.end()) continue;//signals
        //cope with full cut info for each bkg process
		yield_bg     += it->Yield[2][0];
		Err_yield_bg += pow(it->Err_Yield[2][0],2);
        //printf("%8s: X_ori = %3.2e \u00B1 %3.2e, K = %3.2f, X = %3.2e \u00B1 %3.2e, Eff = %3.2e \u00B1 %3.2e, Yield = %3.2e \u00B1 %3.2e\n", 
        printf("%8s: & $%3.2e \\pm %3.2e$ & %3.2f & $%3.2e \\pm %3.2e$ & $%3.2e \\pm %3.2e$ & $%3.2e \\pm %3.2e$\\\\ \n", 
                it->Name, it->X_ori, it->Err_X_ori, it->KFactor, it->X, it->Err_X, it->Eff[2][0], it->Err_Eff[2][0], it->Yield[2][0], it->Err_Yield[2][0]);
        //cope with N-1 cut info for tot bkg process
        for(int i=0; i<NUM; i++){
            sur_nm1[i]   += it->Eff[1][i]*it->Entries;
            tot_nm1[i]   += it->Entries;
            yield_nm1[i] += it->Yield[1][i];
            //printf("%8s(%d): Eff = %8.7f \u00B1 %8.7f, Yield = %12.4f \u00B1 %12.4f\n",
            //        it->Name, i, it->Eff[1][i], it->Err_Eff[1][i], it->Yield[1][i], it->Err_Yield[1][i]);
        }
        if(it+2==vec.end()){
            printf("Total: Yield = $%3.2e \\pm %3.2e$\n", yield_bg, sqrt(Err_yield_bg));
            for(int i=0; i<NUM; i++){
                eff_nm1[i] = sur_nm1[i]/tot_nm1[i];
                Err_eff_nm1[i] = sqrt(eff_nm1[i]*(1-eff_nm1[i])/tot_nm1[i]);
                printf("%9s: %8.0f, Eff = %3.2e \u00B1 %3.2e, Yield = %3.2e\n", histName[i], tot_nm1[i], eff_nm1[i], Err_eff_nm1[i], yield_nm1[i]);
                //printf("%9s: %8.0f, Eff = %8.7f \u00B1 %8.7f, Yield = %12.4f\n", histName[i], tot_nm1[i], eff_nm1[i], Err_eff_nm1[i], yield_nm1[i]);
            }
        }
	}
	//calculate the combined bkg
    const char* Combined[2]={"Combined ppvvv","Combined ppttv"};
	double n[2]={0}, x[2]={0}, ex[2]={0}, y[2]={0}, ey[2]={0}, eff[2]={0}, err[2]={0};//0:ppvvv, 1:ppttv
	for(unsigned int i=0; i<vec.size(); i++){
        if(i==4||i==7)            {n[0]+=vec.at(i).Entries;x[0]+=vec.at(i).X;ex[0]+=pow(vec.at(i).Err_X,2);
                                   y[0]+=vec.at(i).Yield[2][0];ey[0]+=pow(vec.at(i).Err_Yield[2][0],2);}
        if(i==3||i==5||i==6||i==8){n[1]+=vec.at(i).Entries;x[1]+=vec.at(i).X;ex[1]+=pow(vec.at(i).Err_X,2);
                                   y[1]+=vec.at(i).Yield[2][0];ey[1]+=pow(vec.at(i).Err_Yield[2][0],2);}
	}
    for(int i=0; i<2; i++){//0:ppvvv, 1:ppttv
	    ex[i] = sqrt(ex[i]), ey[i] = sqrt(ey[i]);
        eff[i] = y[i]/(L*pb2fb*x[i]); err[i] = sqrt(eff[i]*(1-eff[i])/n[i]);
        printf("%16s: X-sec = %7.4f \u00B1 %8.6f (%6.4f), yield = %5.2f \u00B1 %5.2f (%6.4f), Eff = %8.7f \u00B1 %8.7f, test=%6.4f\n"
                ,Combined[i], x[i], ex[i], ex[i]/x[i], y[i], ey[i], ey[i]/y[i], eff[i], err[i], sqrt(pow(ex[i]/x[i],2)+pow(err[i]/eff[i],2))*L*pb2fb*x[i]*eff[i]);
    }
    //calculate the significance 
	yield_sig = vec.at(0).Yield[2][0];
	significance = sqrt( 2*( (yield_sig + yield_bg)*TMath::Log(1+(yield_sig/yield_bg)) - yield_sig ) );
    printf("significance = %5.4f (yield_sig = %.3f, yield_bg = %.3f)\n", significance, yield_sig, yield_bg);
	//cout<<"significance = "<<significance<<endl;
 

    //TCanvas * can = new TCanvas("can","",800,600);
	//can->SetLogy();
    //TPad    * pad = new TPad("pad","",0,0,1,1);
    //pad->Draw(); pad->cd();
    //vec.at(1).hist_Nm1[0]->Draw();
	//stack_nm1[5]->Draw();

	//### make plot ###//
	line->SetLineColor(kBlack); line->SetLineWidth(3);
	text->SetTextFont(43); text->SetTextSize(16);
	axis->SetLabelSize(0);
	TCanvas *can = new TCanvas("can","",800,600);
	for(int k=0; k<NUM; k++){
		can->SetLogy();
        stack_ori[k]->SetMaximum(MAX_ori);
        stack_ori[k]->SetMinimum(MIN);
		stack_ori[k]->Draw("hist");
        //hist_err_ori[k]->Draw("E2,same");
		vec.at(0).hist_Ori[k]->Draw("hist,same");//"hf,same"
		vec.at(9).hist_Ori[k]->Draw("hist,same");//"hf,same"
        legend_MC_sig->Draw("same");
		legend_MC_bkg->Draw("same");
		binMin = stack_ori[k]->GetXaxis()->GetXmin();
		binMax = stack_ori[k]->GetXaxis()->GetXmax();
		text->DrawLatex((binMax-binMin)*0.01,3e+07,"MC Simulation from MadGraph5");
		text->DrawLatex((binMax-binMin)*0.78,3e+07,"L = 100 fb^{-1} (13TeV)");
	    //axis->DrawAxis(0,MAX_ori,binMax,MAX_ori,0,binMax,510,"-");
	    //axis->DrawAxis(binMax,MIN,binMax,MAX_ori,MIN,MAX_ori,510,"+");
		can->SaveAs(Form("./output/stack_ori_%s.png",histName[k]));
		//can->SaveAs(Form("./output/stack_ori_%s.pdf",histName[k]));
		//can->SaveAs(Form("./output/stack_ori_%s.root",histName[k]));

        stack_nm1[k]->SetMaximum(MAX_nm1);
        stack_nm1[k]->SetMinimum(MIN);
		stack_nm1[k]->Draw("hist");
        //hist_err_nm1[k]->Draw("E2,same");
		vec.at(0).hist_Nm1[k]->Draw("hist,same");
		vec.at(9).hist_Nm1[k]->Draw("hist,same");//"hf,same"
        legend_MC_sig->Draw("same");
		legend_MC_bkg->Draw("same");
		line->DrawLine(CUT[k],0,CUT[k],2e+04);
		line->DrawLine(CUT[k],0,CUT[k],stack_nm1[k]->GetMaximum());
		binMin = stack_nm1[k]->GetXaxis()->GetXmin();
		binMax = stack_nm1[k]->GetXaxis()->GetXmax();
		text->DrawLatex((binMax-binMin)*0.01,2.3e+04,"MC Simulation from MadGraph5");
		text->DrawLatex((binMax-binMin)*0.78,2.3e+04,"L = 100 fb^{-1} (13TeV)");
	    //axis->DrawAxis(0,MAX_nm1,binMax,MAX_nm1,0,binMax,510,"-");
	    //axis->DrawAxis(binMax,MIN,binMax,MAX_nm1,MIN,MAX_nm1,510,"+");
		can->SaveAs(Form("./output/stack_nm1_%s.png",histName[k]));
		//can->SaveAs(Form("./output/stack_nm1_%s.pdf",histName[k]));
		//can->SaveAs(Form("./output/stack_nm1_%s.root",histName[k]));
	}
}
