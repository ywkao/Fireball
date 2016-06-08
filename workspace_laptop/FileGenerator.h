#include <stdio.h> //FILE
FILE *file;
const int observation = 300;//bkg=147.84; 150, 165, 225, 300 
void INFOLOG(const char* name, const char *message){
    file = fopen(Form("%s.txt",name),"a");
    fputs(Form("%s\n",message),file);
    //std::cout<<message<<std::endl;
    fclose(file);
}
void Begin(const char* name){
    file = fopen(Form("%s.txt",name),"w");
    fputs("# Simple counting experiment, with one signal and a few background processes \n",file);
    fputs("imax 1  number of channels\n",file);
    fputs("jmax 4  number of backgrounds\n",file);
    fputs("kmax 10  number of nuisance parameters (sources of systematical uncertainties)\n",file);
    fputs("------------\n",file);
    fputs("# we have just one channel\n",file);
    fputs("bin 1\n",file);
    fputs(Form("observation %d \n",observation),file);
    fputs("------------\n",file);
    fputs("bin              1         1       1        1        1   \n",file);
    fputs("process       fireball    pptt    ppvv    ppvtt    ppvvv \n",file);
    fputs("process          0         1       2        3        4   \n",file);
    fclose(file);
}
