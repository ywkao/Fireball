#include <TRandom3.h>
#include <TH1D.h>

const int BKG  = 8;
const double iteration_N = 100;
class Features{
public:
    double mean_lumi, rms_lumi;//gaus
    double mean_N;//poisson
    double N_s, eff_s;//binomial
    double N_b[BKG], eff_b[BKG];//binomial
    double mean_xsec[BKG], rms_xsec[BKG];//gaus
}; Features f;

std::pair<double, double> ana_uncertainty(Features f, TH1D *& hist){
    delete gRandom;
    gRandom = new TRandom3(100);

    double xsec_ToBeMeasured;
    double luminosity_variation;//
    double total_events_variation;//
    double total_sig_variation = 0;
    double total_bkg_variation = 0;
    double eff_sig_variation;//
    double eff_bkg_variation[BKG];//
    double xsec_bkg_variation[BKG];//

    for(int i=0; i<iteration_N; i++){
        luminosity_variation   = gRandom->Gaus(f.mean_lumi, f.rms_lumi);
        total_events_variation = gRandom->Poisson(f.mean_N);
        eff_sig_variation      = gRandom->Binomial(f.N_s,f.eff_s); eff_sig_variation /= f.N_s;
        total_bkg_variation = 0;
        for(int k=0; k<BKG; k++){
            eff_bkg_variation[k]  = gRandom->Binomial(f.N_b[k],f.eff_b[k]); eff_bkg_variation[k] /= f.N_b[k];
            xsec_bkg_variation[k] = gRandom->Gaus(f.mean_xsec[k], f.rms_xsec[k]);
            total_bkg_variation  += luminosity_variation*xsec_bkg_variation[k]*eff_bkg_variation[k];
        }
        total_sig_variation = total_events_variation - total_bkg_variation;
        xsec_ToBeMeasured   = total_sig_variation/(luminosity_variation*eff_sig_variation);
        hist->Fill(xsec_ToBeMeasured);
        //printf("xsec = %2.1e, total_bkg_variation = %2.1e\t",xsec_ToBeMeasured, total_bkg_variation);
    }

    double mean = hist->GetMean();
    double rms  = hist->GetRMS();
    printf("mean = %2.1e, rms = %2.1e\n", mean, rms);
    return make_pair(mean,rms);
}
