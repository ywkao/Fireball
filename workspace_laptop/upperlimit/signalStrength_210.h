const int NUM=15;
const double MASS_min = 0.4, MASS_max = 2.2;
const double K=1.69;

double MASS[NUM]          = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};//TeV
double Err_MASS[NUM]      = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double Xsec[NUM]          = {8.03e-01, 3.21e-01, 1.40e-01, 6.57e-02, 3.24e-02, 1.67e-02, 8.90e-03, 4.88e-03, 2.74e-03, 1.56e-03, 9.05e-04, 5.30e-04, 3.15e-04, 1.88e-04, 1.13e-04};
double Err_Xsec[NUM]      = {6.45e-04, 2.72e-04, 1.71e-04, 5.30e-05, 3.04e-05, 1.77e-05, 7.78e-06, 5.02e-06, 2.63e-06, 1.31e-06, 1.20e-06, 5.70e-07, 2.70e-07, 1.60e-07, 1.10e-07};
double Xmeasured[NUM]     = {1.36e+00, 5.42e-01, 2.37e-01, 1.11e-01, 5.48e-02, 2.82e-02, 1.50e-02, 8.24e-03, 4.63e-03, 2.64e-03, 1.53e-03, 8.96e-04, 5.32e-04, 3.18e-04, 1.91e-04};
double Err_Xmeasured[NUM] = {8.03e-02, 3.21e-02, 1.40e-02, 6.57e-03, 3.24e-03, 1.67e-03, 8.90e-04, 4.88e-04, 2.74e-04, 1.56e-04, 9.05e-05, 5.30e-05, 3.15e-05, 1.88e-05, 1.13e-05};
const double obs[NUM]   = {0.3781, 0.3394, 0.3756, 0.4662, 0.5424, 0.7770, 1.0928, 1.6739, 2.3278, 3.5825, 5.8001, 9.0735, 13.8795, 22.4988, 36.8458};
const double mean[NUM]  = {0.1646, 0.1499, 0.1665, 0.2061, 0.2412, 0.3447, 0.4863, 0.7441, 1.0352, 1.5859, 2.5703, 4.0469, 6.1719, 10.0312, 16.4375};
const double s1[2][NUM] = {{0.1160, 0.1064, 0.1192, 0.1481, 0.1727, 0.2481, 0.3483, 0.5329, 0.7413, 1.1480, 1.8468, 2.8980, 4.4198, 7.1836, 11.7712},
                           {0.2387, 0.2126, 0.2349, 0.2923, 0.3402, 0.4863, 0.6860, 1.0497, 1.4519, 2.2371, 3.6256, 5.6761, 8.6566, 14.0697, 23.0551}};
const double s2[2][NUM] = {{0.0865, 0.0799, 0.0894, 0.1115, 0.1296, 0.1852, 0.2612, 0.3997, 0.5560, 0.8580, 1.3906, 2.1736, 3.3150, 5.3879, 8.8287},
                           {0.3378, 0.2934, 0.3229, 0.3981, 0.4647, 0.6597, 0.9307, 1.4240, 1.9752, 3.0349, 4.9187, 7.7220, 11.7769, 19.1411, 31.1548}};
int Nobs = 210;
