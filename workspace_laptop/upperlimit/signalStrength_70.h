const int NUM=15;
const double MASS_min = 0.4, MASS_max = 2.2;
const double K=1.69;

double MASS[NUM]          = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};//TeV
double Err_MASS[NUM]      = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double Xsec[NUM]          = {8.03e-01, 3.21e-01, 1.40e-01, 6.57e-02, 3.24e-02, 1.67e-02, 8.90e-03, 4.88e-03, 2.74e-03, 1.56e-03, 9.05e-04, 5.30e-04, 3.15e-04, 1.88e-04, 1.13e-04};
double Err_Xsec[NUM]      = {6.45e-04, 2.72e-04, 1.71e-04, 5.30e-05, 3.04e-05, 1.77e-05, 7.78e-06, 5.02e-06, 2.63e-06, 1.31e-06, 1.20e-06, 5.70e-07, 2.70e-07, 1.60e-07, 1.10e-07};
double Xmeasured[NUM]     = {1.36e+00, 5.42e-01, 2.37e-01, 1.11e-01, 5.48e-02, 2.82e-02, 1.50e-02, 8.24e-03, 4.63e-03, 2.64e-03, 1.53e-03, 8.96e-04, 5.32e-04, 3.18e-04, 1.91e-04};
double Err_Xmeasured[NUM] = {8.03e-02, 3.21e-02, 1.40e-02, 6.57e-03, 3.24e-03, 1.67e-03, 8.90e-04, 4.88e-04, 2.74e-04, 1.56e-04, 9.05e-05, 5.30e-05, 3.15e-05, 1.88e-05, 1.13e-05};
const double obs[NUM]   = {0.0335, 0.0305, 0.0343, 0.0426, 0.0493, 0.0708, 0.1001, 0.1529, 0.2126, 0.3270, 0.5304, 0.8296, 1.2697, 2.0539, 3.3684};
const double mean[NUM]  = {0.1177, 0.1069, 0.1187, 0.1479, 0.1724, 0.2471, 0.3467, 0.5332, 0.7402, 1.1367, 1.8359, 2.8828, 4.4062, 7.1562, 11.7188};
const double s1[2][NUM] = {{0.0819, 0.0757, 0.0840, 0.1052, 0.1220, 0.1757, 0.2474, 0.3793, 0.5265, 0.8113, 1.3103, 2.0506, 3.1368, 5.0904, 8.3359},
                           {0.1716, 0.1534, 0.1702, 0.2111, 0.2459, 0.3505, 0.4945, 0.7564, 1.0500, 1.6215, 2.6190, 4.0894, 6.2504, 10.1514, 16.6235}};
const double s2[2][NUM] = {{0.0611, 0.0564, 0.0626, 0.0783, 0.0912, 0.1308, 0.1849, 0.2822, 0.3918, 0.6061, 0.9789, 1.5259, 2.3494, 3.7878, 6.2027},
                           {0.2451, 0.2145, 0.2350, 0.2904, 0.3383, 0.4836, 0.6804, 1.0436, 1.4488, 2.2165, 3.6031, 5.6422, 8.6239, 13.9155, 22.7874}};
int Nobs = 70;
