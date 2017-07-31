/*
 *  fitData.h
 *  FissionYields
 *
 *  Created by Ionel Stetcu on 2/3/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */



/*
 #define u236_minHFFmass 118
 #define u236_nmass 40
 
 //Djanchenko data
static double u236_avTKE[u236_nmass]={160.97,159.94,157.25,156.38,158.33,164.43,169.54,172.35,173.53,174.52,
                                      176.22,178.37,179.66,179.69,179.59,179.35,178.80,177.76,176.39,174.89, 
                                      173.66,172.6 ,171.4 ,170.17,169.02,168.13,167.13,166.29,165.34,164.14, 
                                      163.18,162.19,161.32,160.46,159.6 ,159.09,158.82,157.2 ,154.37,153.20 };
 
 static double u236_sTKE[u236_nmass]={10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
	10., 10., 10., 10., 10., 10., 10., 10., 10., 10. };
*/

// Hambsch data (from YATKE file)

#define u236_minHFFmass 118
#define u236_nmass 48

static double u236_avTKE[u236_nmass]={151.5000,166.1429,151.2667,156.5263,159.2500,164.1250,171.8857,173.1961,173.2660,173.5195,
    176.8648,177.2044,178.0114,178.3366,178.1227,178.0794,177.6892,176.1792,174.6586,173.2974,
    172.4086,171.9007,169.9520,169.0901,168.5960,167.5808,166.3030,165.7135,164.9405,163.6641,
    163.3955,161.7382,160.4763,159.4670,158.7371,157.7636,157.1034,155.9000,156.6667,152.0000,
    149.4737,150.1538,149.8571,149.6667,146.0000,148.0000,150.0000,149.0000};

static double u236_sTKE[u236_nmass]={14.1038,6.8542,15.0176,19.3455,15.2705,13.2769,9.9564,12.3796,12.6038,14.1170,
    10.4094,10.7006,11.1495,10.4034,10.3023,9.9951,10.0074,9.3691,9.3979,9.2273,
    8.8022,8.4325,8.1027,8.0632,7.7664,7.7215,7.4167,7.5380,7.1638,7.2199,
    7.3231,6.9727,6.9534,7.3366,6.8193,6.2040,9.1852,6.6038,5.5817,6.2941,
    7.5768,7.6846,6.4902,7.4714,7.0000,7.0000,7.0000,7.0000}; // change all the zero's by 7.00



#define u239_minHFFmass 120
#define u239_nmass 40

static double u239_avTKE[u239_nmass]={160.97,159.94,157.25,156.38,158.33,164.43,169.54,172.35,173.53,174.52,  // for now this is U236 data
	176.22,178.37,179.66,179.69,179.59,179.35,178.80,177.76,176.39,174.89, 
	173.66,172.6 ,171.4 ,170.17,169.02,168.13,167.13,166.29,165.34,164.14, 
	163.18,162.19,161.32,160.46,159.6 ,159.09,158.82,157.2 ,154.37,153.20 };

static double u239_sTKE[u239_nmass]={10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
	10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
	10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
	10., 10., 10., 10., 10., 10., 10., 10., 10., 10. };

#define pu240_minHFFmass 120
#define pu240_nmass 51

static double pu240_avTKE[pu240_nmass]=
  {1.6979E+02,1.7017E+02,1.7227E+02,1.7498E+02,1.7667E+02,1.7920E+02,1.8076E+02,1.8241E+02,1.8338E+02,1.8427E+02,
   1.8464E+02,1.8488E+02,1.8489E+02,1.8475E+02,1.8424E+02,1.8351E+02,1.8246E+02,1.8140E+02,1.8043E+02,1.7919E+02,
   1.7791E+02,1.7672E+02,1.7542E+02,1.7449E+02,1.7340E+02,1.7227E+02,1.7104E+02,1.6990E+02,1.6875E+02,1.6768E+02,
   1.6655E+02,1.6541E+02,1.6432E+02,1.6337E+02,1.6254E+02,1.6168E+02,1.6054E+02,1.5979E+02,1.5889E+02,1.5772E+02,
   1.5659E+02,1.5459E+02,1.5360E+02,1.5224E+02,1.5086E+02,1.4991E+02,1.4788E+02,1.4653E+02,1.4319E+02,1.4264E+02,
   1.3925E+02};



static double pu240_sTKE[pu240_nmass]={10.439,11.327,12.190,12.543,12.079,11.615,11.101,10.918,10.863,10.884,
				       10.829,10.774,10.769,10.816,10.812,10.680,10.369,10.187,10.029,9.8210,
                                       9.5109,9.2001,8.9663,8.7065,8.6260,8.4689,8.4643,8.3069,8.2005,8.1452,
                                       8.1918,8.1110,8.1829,8.1278,8.1489,8.2469,8.3195,8.4167,8.4636,8.3830,
                                       8.2772,8.6039,8.1146,8.4676,8.1059,7.7190,7.7190,7.7190,7.7190,7.7190,
                                       7.7190};  // data by M.Asghar, F.Caitucoli, P.Perrin, C.Wagemans: REFERENCE  (J,NP/A,311,205,197811)
                                                // the last 5 are copies of the sigma_TKE^2(A=165), to account for missing experimental data

