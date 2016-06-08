#include <TLine.h>
#include <TLatex.h>
TLatex *text = new TLatex(0,0,"");
TLine  *line = new TLine(0,0,0,0);
TLegend*legend = new TLegend(0.7,0.7,0.9,0.9);

const double UPPER_BOUND=1e+05, LOWER_BOUND=1e-14;
const double NUM_min = 0.4, NUM_max = 2.2;
//double Z[5]        = {1, 2, 3, 4, 5};
//double Z_p_value[5]= {1.5866e-01, 2.2750e-02, 1.3499e-03, 3.1671e-05, 2.8665e-07};
//const int Color[5] = {kBlack, kOrange+2, kBlue, kGreen+2, kRed};
//const int Color[5] = {kGray, kGray,kGray,kGray,kGray};

//const int NUM = 11;
//const double MASS[NUM]   = {1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};//TeV
//double significance[NUM] = {1.34e+01, 9.59e+00, 7.01e+00, 4.80e+00, 3.43e+00, 2.26e+00, 1.41e+00, 9.13e-01, 5.98e-01, 3.70e-01, 2.27e-01};
//double p_value[NUM]      = {3.12e-41, 4.28e-22, 1.22e-12, 7.83e-07, 3.01e-04, 1.18e-02, 7.90e-02, 1.81e-01, 2.75e-01, 3.56e-01, 4.10e-01};

const int NUM_Z = 7;
double Z[NUM_Z]        = {1, 2, 3, 4, 5, 6, 7};
double Z_p_value[NUM_Z]= {1.59e-01, 2.27e-02, 1.35e-03, 3.17e-05, 2.87e-07, 9.87e-10, 1.28e-12};
const int Color[NUM_Z] = {kRed, kRed, kRed, kRed, kRed, kRed, kRed};

const int NUM = 15;
const double MASS[NUM]   = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};//TeV
double significance[NUM] = {1.90e+01, 2.02e+01, 1.83e+01, 1.51e+01, 1.33e+01, 9.60e+00, 7.01e+00, 4.71e+00, 3.44e+00, 2.27e+00, 1.42e+00, 9.11e-01, 5.98e-01, 3.70e-01, 2.27e-01};
double p_value[NUM]      = {1.42e-80, 8.21e-91, 3.09e-75, 4.39e-52, 1.78e-40, 4.17e-22, 1.23e-12, 1.26e-06, 2.94e-04, 1.17e-02, 7.85e-02, 1.81e-01, 2.75e-01, 3.56e-01, 4.10e-01};
double significance_pkg[NUM] = {9.098, 9.563, 8.849, 7.612, 6.834, 5.191, 3.928, 2.726, 2.027, 1.359, 0.859, 0.557, 0.367, 0.228, 0.140};
double p_value_pkg[NUM]      = {4.58e-20, 5.73e-22, 4.41e-19, 1.35e-14, 4.11e-12, 1.04e-07, 4.28e-05, 3.21e-03, 2.13e-02, 8.70e-02, 1.95e-01, 2.89e-01, 3.57e-01, 4.10e-01, 4.44e-01};
