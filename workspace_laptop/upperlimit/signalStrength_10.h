const int NUM=15;
const double MASS_min = 0.4, MASS_max = 2.2;
const double K=1.69;

double MASS[NUM]          = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};//TeV
double Err_MASS[NUM]      = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double Xsec[NUM]          = {8.03e-01, 3.21e-01, 1.40e-01, 6.57e-02, 3.24e-02, 1.67e-02, 8.90e-03, 4.88e-03, 2.74e-03, 1.56e-03, 9.05e-04, 5.30e-04, 3.15e-04, 1.88e-04, 1.13e-04};
double Err_Xsec[NUM]      = {6.45e-04, 2.72e-04, 1.71e-04, 5.30e-05, 3.04e-05, 1.77e-05, 7.78e-06, 5.02e-06, 2.63e-06, 1.31e-06, 1.20e-06, 5.70e-07, 2.70e-07, 1.60e-07, 1.10e-07};
double Xmeasured[NUM]     = {1.36e+00, 5.42e-01, 2.37e-01, 1.11e-01, 5.48e-02, 2.82e-02, 1.50e-02, 8.24e-03, 4.63e-03, 2.64e-03, 1.53e-03, 8.96e-04, 5.32e-04, 3.18e-04, 1.91e-04};
double Err_Xmeasured[NUM] = {8.03e-02, 3.21e-02, 1.40e-02, 6.57e-03, 3.24e-03, 1.67e-03, 8.90e-04, 4.88e-04, 2.74e-04, 1.56e-04, 9.05e-05, 5.30e-05, 3.15e-05, 1.88e-05, 1.13e-05};
const double obs[NUM]   = {0.0145, 102.0078, 0.0148, 0.0185, 0.0214, 0.0307, 0.0435, 0.0663, 0.0922, 0.1421, 0.2299, 0.3594, 0.5507, 0.8929, 1.4625};
const double mean[NUM]  = {0.1011, 0.0923, 0.1021, 0.1274, 0.1479, 0.2119, 0.2998, 0.4590, 0.6387, 0.9805, 1.5859, 2.4766, 3.7969, 6.1719, 10.0938};
const double s1[2][NUM] = {{0.0706, 0.0649, 0.0723, 0.0903, 0.1049, 0.1507, 0.2125, 0.3254, 0.4504, 0.6963, 1.1223, 1.7616, 2.6869, 4.3754, 7.1558},
                           {0.1490, 0.1331, 0.1472, 0.1828, 0.2134, 0.3040, 0.4277, 0.6584, 0.9111, 1.4065, 2.2750, 3.5526, 5.4465, 8.8042, 14.3988}};
const double s2[2][NUM] = {{0.0517, 0.0479, 0.0538, 0.0672, 0.0777, 0.1122, 0.1575, 0.2411, 0.3356, 0.5170, 0.8394, 1.3108, 2.0097, 3.2426, 5.3032},
                           {0.2140, 0.1867, 0.2052, 0.2540, 0.2956, 0.4223, 0.5922, 0.9090, 1.2615, 1.9417, 3.1408, 4.9046, 7.5193, 12.1903, 19.9365}};
int Nobs = 10;
