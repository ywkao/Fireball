const int NUM=15;
const double MASS_min = 0.4, MASS_max = 2.2;
const double K=1.69;

double MASS[NUM]          = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};//TeV
double Err_MASS[NUM]      = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double Xsec[NUM]          = {8.03e-01, 3.21e-01, 1.40e-01, 6.57e-02, 3.24e-02, 1.67e-02, 8.90e-03, 4.88e-03, 2.74e-03, 1.56e-03, 9.05e-04, 5.30e-04, 3.15e-04, 1.88e-04, 1.13e-04};
double Err_Xsec[NUM]      = {6.45e-04, 2.72e-04, 1.71e-04, 5.30e-05, 3.04e-05, 1.77e-05, 7.78e-06, 5.02e-06, 2.63e-06, 1.31e-06, 1.20e-06, 5.70e-07, 2.70e-07, 1.60e-07, 1.10e-07};
double Xmeasured[NUM]     = {1.36e+00, 5.42e-01, 2.37e-01, 1.11e-01, 5.48e-02, 2.82e-02, 1.50e-02, 8.24e-03, 4.63e-03, 2.64e-03, 1.53e-03, 8.96e-04, 5.32e-04, 3.18e-04, 1.91e-04};
double Err_Xmeasured[NUM] = {8.03e-02, 3.21e-02, 1.40e-02, 6.57e-03, 3.24e-03, 1.67e-03, 8.90e-04, 4.88e-04, 2.74e-04, 1.56e-04, 9.05e-05, 5.30e-05, 3.15e-05, 1.88e-05, 1.13e-05};
const double obs[NUM]   = {0.0810, 0.0739, 0.0826, 0.1027, 0.1193, 0.1711, 0.2413, 0.3701, 0.5131, 0.7912, 1.2814, 2.0035, 3.0664, 4.9727, 8.1436};
const double mean[NUM]  = {0.1323, 0.1206, 0.1343, 0.1675, 0.1948, 0.2783, 0.3926, 0.6035, 0.8359, 1.2852, 2.0859, 3.2656, 4.9844, 8.0938, 13.1875};
const double s1[2][NUM] = {{0.0932, 0.0854, 0.0956, 0.1192, 0.1387, 0.1986, 0.2802, 0.4285, 0.5966, 0.9172, 1.4888, 2.3307, 3.5575, 5.7767, 9.4438},
                           {0.1930, 0.1730, 0.1915, 0.2376, 0.2764, 0.3948, 0.5569, 0.8513, 1.1858, 1.8230, 2.9423, 4.6064, 7.0705, 11.4167, 18.7069}};
const double s2[2][NUM] = {{0.0690, 0.0636, 0.0716, 0.0893, 0.1039, 0.1484, 0.2093, 0.3218, 0.4457, 0.6852, 1.1122, 1.7412, 2.6577, 4.3156, 7.0831},
                           {0.2756, 0.2404, 0.2635, 0.3257, 0.3788, 0.5412, 0.7634, 1.1703, 1.6149, 2.4827, 4.0183, 6.2909, 9.6290, 15.5918, 25.6434}};
int Nobs = 120;
