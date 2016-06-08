
//=== full bkg info ===//
void MakeTable_bkg(std::vector<MyProcess> vec_bkg){
    printf(" Process &   lumi   &   xsec   &   K-factor   &   eff   &   Yield \n");
	for(std::vector<MyProcess>::iterator it=vec_bkg.begin(); it!=vec_bkg.end(); it++){
        printf("%8s: lumi = %.0f \u00B1 %.0f, xsec = %3.2e, K-factor = %3.2f \u00B1 %3.2f, eff = %3.2e \u00B1 %2.1e, Yield = %3.2e \u00B1 %3.2e \n", 
        //printf("%8s & $%.0f \\pm %.0f$ & $%3.2e$ & $%3.2f \\pm %3.2f$ & $%3.2e \\pm %3.2e$ & $%3.0f \\pm %3.0f$ \\\\\n", 
                it->Name, L, L*Err_L, it->X, it->KFactor, it->Err_KFactor, it->Eff[2][0], it->Err_Eff[2][0], it->Yield[2][0], it->Err_Yield[2][0]);
    }
    printf("\n");
}

//=== full cut info ===//
void MakeTable_sig(std::vector<MyProcess> vec){
    printf(" Process &   xsec_ori   &   xsec   &   eff   &   Yield \n"); 
	for(std::vector<MyProcess>::iterator it=vec.begin(); it!=vec.end(); it++){
        printf("%16s: & $%3.2e \\pm %3.2e$ & $%3.2e \\pm %3.2e$ & $%3.2e \\pm %3.2e$ & $%3.2e \\pm %3.2e$\\\\ \n", 
                it->Name, it->X_ori, it->Err_X_ori, it->X, it->Err_X, it->Eff[2][0], it->Err_Eff[2][0], it->Yield[2][0], it->Err_Yield[2][0]);
    }
    printf("\n");
}

//===  N-1 info  ===//
void MakeTable_nm1(std::vector<MyProcess> vec){
	for(std::vector<MyProcess>::iterator it=vec.begin(); it!=vec.end(); it++){
        if((it+1!=vec.end())) continue;
        for(int i=0; i<NUM; i++){
            if(i==2||i==3||i==8) continue;
            printf("%8s(%d): %8.7f $\\pm$ %8.7f,%9.4f $\\pm$ %9.4f\n", it->Name, i, it->Eff[1][i], it->Err_Eff[1][i], it->Yield[1][i], it->Err_Yield[1][i]);
        }
    }
    printf("\n");
}

