//void Processing(const char* NAME){
//	gROOT->ProcessLine(Form(".L %s.C++",NAME));
//	gSystem->Load(Form("%s_C.so",NAME));
//	ana_delphes t;	t.Loop();
//	gROOT->ProcessLine(".q");
//}

void Ana(const char* TARGET){
	if(TARGET=="pythia"){
		//gROOT->ProcessLine(".L ana_pythia.C++");
		//gSystem->Load("ana_pythia_C.so");
		//ana_pythia t;
		//t.Loop();
		gROOT->ProcessLine(".L ana_pythia_500.C++");
		gSystem->Load("ana_pythia_500_C.so");
		ana_pythia_500 t_500;	t_500.Loop();

		gROOT->ProcessLine(".L ana_pythia_1000.C++");
		gSystem->Load("ana_pythia_1000_C.so");
		ana_pythia_1000 t_1000;	t_1000.Loop();

		gROOT->ProcessLine(".L ana_pythia_1500.C++");
		gSystem->Load("ana_pythia_1500_C.so");
		ana_pythia_1500 t_1500;	t_1500.Loop();

		gROOT->ProcessLine(".L ana_pythia_2000.C++");
		gSystem->Load("ana_pythia_2000_C.so");
		ana_pythia_2000 t_2000;	t_2000.Loop();

		gROOT->ProcessLine(".L ana_pythia_2500.C++");
		gSystem->Load("ana_pythia_2500_C.so");
		ana_pythia_2500 t_2500;	t_2500.Loop();
	} else if(TARGET=="unweighted"){
		gROOT->ProcessLine(".L ana_unweighted.C++");
		gSystem->Load("ana_unweighted_C.so");
		ana_unweighted t;
		t.Loop();
	} else if(TARGET=="delphes"){
		gROOT->ProcessLine(".L ana_delphes_500.C++");
		gSystem->Load("ana_delphes_500_C.so");
		ana_delphes_500 t_500;	t_500.Loop();

		gROOT->ProcessLine(".L ana_delphes_1000.C++");
		gSystem->Load("ana_delphes_1000_C.so");
		ana_delphes_1000 t_1000;	t_1000.Loop();

		gROOT->ProcessLine(".L ana_delphes_1500.C++");
		gSystem->Load("ana_delphes_1500_C.so");
		ana_delphes_1500 t_1500;	t_1500.Loop();

		gROOT->ProcessLine(".L ana_delphes_2000.C++");
		gSystem->Load("ana_delphes_2000_C.so");
		ana_delphes_2000 t_2000;	t_2000.Loop();

		gROOT->ProcessLine(".L ana_delphes_2500.C++");
		gSystem->Load("ana_delphes_2500_C.so");
		ana_delphes_2500 t_2500;	t_2500.Loop();
	} else cout<<"[ERROR] The given TARGET is incorrect."<<endl;
}
