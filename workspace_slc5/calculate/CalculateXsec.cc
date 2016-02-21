#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]){
	const int datapoint = atoi(argv[1]);
	FILE *logfile = fopen("log_xsec","r");
	double data1, data2, xsec=0, err=0;
	for(int i=0; i<datapoint; i++){
		fscanf(logfile,"%lf %lf",&data1, &data2);
		xsec += data1;
		err  += pow(data2,2);
	}
	xsec /= (double)datapoint;
	err   = sqrt(err)/(double)datapoint;
	cout<<"Ave. xsec = "<<xsec<<" \u00B1 "<<err<<endl;

}
