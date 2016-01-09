#include <iostream>
#include "ana_delphes.h"
using namespace std;

const char* Name;
int main(int argc, char *argv[]){
	Name= argv[1];
	//X = Xsec.GetCrossSection(Name);
	cout<<Name<<endl;

	if((string)Name=="pptttt") cout<<"Shit!!!"<<endl;
	X = Xsec.GetCrossSection("pptttt");
	cout<< "X = " << X << endl;

	return 1;
}
