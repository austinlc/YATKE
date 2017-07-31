/*
 *  Yields.h
 *  FissionYields
 *
 *  Created by Ionel Stetcu on 7/17/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

class Yields {

public:

	Yields (int id);
	~Yields (void);

	void setYieldParameters(void);
	void setFissioningSystem(int id,double ex);
	double yieldA(int A);
	double yieldATKE(int A,double TKE);
	double yieldTKE(int A,double TKE);
	double sampleTKE(int A);
	double get_avTKE(int A);
	void setRescaleTKE(double *,int,int);

private:

	int numModes;
	int id_cn0;
	int Acn0,Zcn0;
	int a_cn;
	double n0[5],d[5],sig0[5],sig1[5],sig2[5],e1[5],e2[5];
	double w_e[5],sig_e[5],Abar[5];
	double gaussian(double x, double w, double av, double s);
	double avTKE[300];
	double sTKE[300];
	double sf_flag;
	double einc_ref;
	double tke_av_ref;
	void setParametersAverageTKE(void);
	double TKE_fit_param[5]; int n_poly; // average TKE parameterization
	double SystViola1966(void);
	double SystViola1966(int);
	double NewSyst(int);
	double averageTKE(double e); // calculates average TKE for the system as a function of incident energy
	double rescaleAverageTKE;
	double sampleGaussian(void);
	double einc_now;

};

