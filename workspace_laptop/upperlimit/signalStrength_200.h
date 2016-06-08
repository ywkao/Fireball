const int NUM=15;
const double MASS_min = 0.4, MASS_max = 2.2;
const double K=1.69;

double MASS[NUM]          = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};//TeV
double Err_MASS[NUM]      = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double Xsec[NUM]          = {8.03e-01, 3.21e-01, 1.40e-01, 6.57e-02, 3.24e-02, 1.67e-02, 8.90e-03, 4.88e-03, 2.74e-03, 1.56e-03, 9.05e-04, 5.30e-04, 3.15e-04, 1.88e-04, 1.13e-04};
double Err_Xsec[NUM]      = {6.45e-04, 2.72e-04, 1.71e-04, 5.30e-05, 3.04e-05, 1.77e-05, 7.78e-06, 5.02e-06, 2.63e-06, 1.31e-06, 1.20e-06, 5.70e-07, 2.70e-07, 1.60e-07, 1.10e-07};
double Xmeasured[NUM]     = {1.36e+00, 5.42e-01, 2.37e-01, 1.11e-01, 5.48e-02, 2.82e-02, 1.50e-02, 8.24e-03, 4.63e-03, 2.64e-03, 1.53e-03, 8.96e-04, 5.32e-04, 3.18e-04, 1.91e-04};
double Err_Xmeasured[NUM] = {8.03e-02, 3.21e-02, 1.40e-02, 6.57e-03, 3.24e-03, 1.67e-03, 8.90e-04, 4.88e-04, 2.74e-04, 1.56e-04, 9.05e-05, 5.30e-05, 3.15e-05, 1.88e-05, 1.13e-05};
const double obs[NUM]   = {0.3364, 0.3020, 0.3344, 0.4157, 0.4830, 0.6911, 0.9736, 1.4916, 2.0739, 3.1917, 5.1679, 8.0658, 12.3632, 20.0192, 32.8262};
const double mean[NUM]  = {0.1606, 0.1460, 0.1626, 0.2021, 0.2354, 0.3369, 0.4746, 0.7285, 1.0117, 1.5586, 2.5234, 3.9531, 6.0156, 9.7812, 16.0625};
const double s1[2][NUM] = {{0.1133, 0.1042, 0.1164, 0.1448, 0.1685, 0.2413, 0.3416, 0.5217, 0.7245, 1.1161, 1.8071, 2.8309, 4.3223, 7.0045, 11.5026},
                           {0.2330, 0.2083, 0.2307, 0.2851, 0.3320, 0.4752, 0.6695, 1.0276, 1.4271, 2.1861, 3.5393, 5.5446, 8.4854, 13.7191, 22.5291}};
const double s2[2][NUM] = {{0.0844, 0.0778, 0.0873, 0.1086, 0.1264, 0.1810, 0.2549, 0.3913, 0.5434, 0.8371, 1.3554, 2.1233, 3.2545, 5.2536, 8.6273},
                           {0.3298, 0.2865, 0.3141, 0.3894, 0.4534, 0.6447, 0.9082, 1.3941, 1.9361, 2.9740, 4.8151, 7.5431, 11.5119, 18.6641, 30.6496}};
int Nobs = 200;
