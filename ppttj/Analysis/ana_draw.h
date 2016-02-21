int XERR = 0;
void ScaleAndDrawSetting(TH1D *&hist, double &Yield, double &Err, double factor, double X, double Err_X, int lcolor, int lstyle, int lwidth, int fcolor, int style){
	double content, error, yield=0, err=0;
	for(int i=0; i<hist->GetNbinsX(); i++){
		content = hist->GetBinContent(i+1);
		error   = hist->GetBinError(i+1);
		yield  += content*factor;
		err    += pow(error*factor,2);

		hist->SetBinContent(i+1,content*factor);
		if(XERR==0) hist->SetBinError(i+1,error*factor);
		else		hist->SetBinError(i+1, (pow(Err_X/X,2)+ pow(error/content,2))*content*factor );
	}
	Yield = yield;
	Err   = sqrt(err);
	hist->SetLineColor(lcolor);
	hist->SetLineStyle(lstyle);
	hist->SetLineWidth(lwidth);
	hist->SetFillColor(fcolor);
	hist->SetFillStyle(style);
}

/*
void HistDivide(TH1D *&hist, TH1D *hist_base){
	double content, error;
	for(int i=0; i<hist->GetNbinsX(); i++){
		content = hist->GetBinContent(i+1);
		error   = hist->GetBinError(i+1);

		hist->SetBinContent(i+1,content/hist_base->GetBinContent(i+1));
		hist->SetBinError(i+1,error*factor);
		hist->SetBinError(i+1, (pow(Err_X/X,2)+ pow(error/content,2))*content*factor );
	}
	
}
*/
