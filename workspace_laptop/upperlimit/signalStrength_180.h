const int NUM=15;
const double MASS_min = 0.4, MASS_max = 2.2;
const double K=1.69;

double MASS[NUM]          = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};//TeV
double Err_MASS[NUM]      = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double Xsec[NUM]          = {8.03e-01, 3.21e-01, 1.40e-01, 6.57e-02, 3.24e-02, 1.67e-02, 8.90e-03, 4.88e-03, 2.74e-03, 1.56e-03, 9.05e-04, 5.30e-04, 3.15e-04, 1.88e-04, 1.13e-04};
double Err_Xsec[NUM]      = {6.45e-04, 2.72e-04, 1.71e-04, 5.30e-05, 3.04e-05, 1.77e-05, 7.78e-06, 5.02e-06, 2.63e-06, 1.31e-06, 1.20e-06, 5.70e-07, 2.70e-07, 1.60e-07, 1.10e-07};
double Xmeasured[NUM]     = {1.36e+00, 5.42e-01, 2.37e-01, 1.11e-01, 5.48e-02, 2.82e-02, 1.50e-02, 8.24e-03, 4.63e-03, 2.64e-03, 1.53e-03, 8.96e-04, 5.32e-04, 3.18e-04, 1.91e-04};
double Err_Xmeasured[NUM] = {8.03e-02, 3.21e-02, 1.40e-02, 6.57e-03, 3.24e-03, 1.67e-03, 8.90e-04, 4.88e-04, 2.74e-04, 1.56e-04, 9.05e-05, 5.30e-05, 3.15e-05, 1.88e-05, 1.13e-05};
const double obs[NUM]   = {0.2539, 0.2296, 0.2543, 0.3165, 0.3675, 0.5261, 0.7412, 1.1378, 1.5769, 2.4267, 3.9359, 6.1583, 9.4167, 15.2712, 25.0092};
const double mean[NUM]  = {0.1528, 0.1392, 0.1548, 0.1929, 0.2256, 0.3213, 0.4551, 0.6934, 0.9648, 1.4883, 2.3984, 3.7656, 5.7656, 9.3438, 15.3125};
const double s1[2][NUM] = {{0.1083, 0.0993, 0.1108, 0.1381, 0.1610, 0.2301, 0.3259, 0.4965, 0.6909, 1.0658, 1.7262, 2.6966, 4.1289, 6.6912, 10.9656},
                           {0.2229, 0.1985, 0.2196, 0.2736, 0.3182, 0.4532, 0.6383, 0.9780, 1.3610, 2.0993, 3.3831, 5.3116, 8.0868, 13.1799, 21.4772}};
const double s2[2][NUM] = {{0.0803, 0.0742, 0.0831, 0.1036, 0.1203, 0.1726, 0.2444, 0.3724, 0.5182, 0.7994, 1.2882, 2.0226, 3.0968, 5.0186, 8.2245},
                           {0.3145, 0.2749, 0.3010, 0.3726, 0.4346, 0.6189, 0.8684, 1.3357, 1.8464, 2.8481, 4.6203, 7.2061, 11.0017, 17.8808, 29.2185}};
int Nobs = 180;