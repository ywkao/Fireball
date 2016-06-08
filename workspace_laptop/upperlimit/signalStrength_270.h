const int NUM=15;
const double MASS_min = 0.4, MASS_max = 2.2;
const double K=1.69;

double MASS[NUM]          = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};//TeV
double Err_MASS[NUM]      = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double Xsec[NUM]          = {8.03e-01, 3.21e-01, 1.40e-01, 6.57e-02, 3.24e-02, 1.67e-02, 8.90e-03, 4.88e-03, 2.74e-03, 1.56e-03, 9.05e-04, 5.30e-04, 3.15e-04, 1.88e-04, 1.13e-04};
double Err_Xsec[NUM]      = {6.45e-04, 2.72e-04, 1.71e-04, 5.30e-05, 3.04e-05, 1.77e-05, 7.78e-06, 5.02e-06, 2.63e-06, 1.31e-06, 1.20e-06, 5.70e-07, 2.70e-07, 1.60e-07, 1.10e-07};
double Xmeasured[NUM]     = {1.36e+00, 5.42e-01, 2.37e-01, 1.11e-01, 5.48e-02, 2.82e-02, 1.50e-02, 8.24e-03, 4.63e-03, 2.64e-03, 1.53e-03, 8.96e-04, 5.32e-04, 3.18e-04, 1.91e-04};
double Err_Xmeasured[NUM] = {8.03e-02, 3.21e-02, 1.40e-02, 6.57e-03, 3.24e-03, 1.67e-03, 8.90e-04, 4.88e-04, 2.74e-04, 1.56e-04, 9.05e-05, 5.30e-05, 3.15e-05, 1.88e-05, 1.13e-05};
const double obs[NUM]   = {0.6406, 0.5700, 0.6284, 0.7789, 0.9047, 1.2937, 1.8192, 2.7940, 3.8740, 5.9606, 9.6502, 15.0981, 23.0621, 37.4563, 61.3403};
const double mean[NUM]  = {0.1880, 0.1704, 0.1899, 0.2354, 0.2744, 0.3926, 0.5566, 0.8477, 1.1836, 1.8203, 2.9453, 4.6094, 7.0312, 11.4062, 18.6875};
const double s1[2][NUM] = {{0.1326, 0.1220, 0.1360, 0.1704, 0.1972, 0.2842, 0.4000, 0.6136, 0.8504, 1.3079, 2.1162, 3.3119, 5.0899, 8.2160, 13.4272},
                           {0.2712, 0.2417, 0.2679, 0.3320, 0.3871, 0.5538, 0.7785, 1.1957, 1.6554, 2.5532, 4.1311, 6.4651, 9.8620, 15.9983, 26.2109}};
const double s2[2][NUM] = {{0.0988, 0.0915, 0.1020, 0.1273, 0.1485, 0.2124, 0.3012, 0.4586, 0.6403, 0.9848, 1.5935, 2.4937, 3.8040, 6.1710, 10.1102},
                           {0.3827, 0.3335, 0.3659, 0.4504, 0.5251, 0.7463, 1.0535, 1.6113, 2.2400, 3.4501, 5.5824, 8.7364, 13.3267, 21.6188, 35.4193}};
int Nobs = 270;
