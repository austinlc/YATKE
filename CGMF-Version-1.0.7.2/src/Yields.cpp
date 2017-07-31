/*
 *  Yields.cpp
 *  FissionYields
 *
 *  Created by Ionel Stetcu on 7/17/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 */
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "Yields.h"
#include "constant.h"
#include "fitData.h"
#include "config-ff.h"
#include "mt19937ar.h"

#include "config.h"

using namespace std;

Yields::Yields (int id){
   id_cn0=id;
   Zcn0=id/1000;
   Acn0=id%1000;
//   cout << Acn0 << " " << Zcn0 << endl;
   
   setYieldParameters();
   // set here the reference energy
   setParametersAverageTKE();
   return;
}

Yields::~Yields (void){
   return;
}

void Yields::setYieldParameters(void){
   d[0]=0.;
   numModes=3;
   string gaussDataFit = DATADIR;
   gaussDataFit += "gaussfit.dat" ;
   ifstream gaussfit(gaussDataFit.c_str(),ios::in); // Parameters from PRC 85 (2012) 024608 or the FREYA input
   if (!gaussfit.is_open()) {
      cout << "The file " << gaussDataFit << " not found" <<endl;
      return;
   }
   for (int i=0; i<7; i++) {
      string line;
      getline(gaussfit, line);
   }
   bool not_found=true;
   while (!gaussfit.eof()) {
      int zz,aa;
      gaussfit >> zz >> aa >> n0[1] >> e1[1] >> e2[1] ;
      gaussfit >> n0[2] >> e1[2] >> e2[2] ;
      gaussfit >> d[1] >> d[2];
      for (int m=0; m<3; ++m) {
         gaussfit >> sig0[m] >> sig1[m] >> sig2[m] ;
      }
      if(id_cn0==zz*1000+aa){
         not_found=false;
         break;
      }
   }
   gaussfit.close();
   if(not_found){
      cout << "Compound system not found" << endl;
      return;
   }
   sf_flag = false;
  
  // TODO: Why do we need this here? --> Ionel
   if(id_cn0==92238 || id_cn0==98252 ){
      cout << "Only spontaneous fission available" << endl;
      sf_flag = true;
      cout << "revert to e_inc=0" << endl;
   }
   
   
   return;
   
}

void Yields::setFissioningSystem(int id,double einc){
   
   //determine first the equivalent incident energy
   int emitted_neutrons=id_cn0-id;
   /*
    double einc=ex-sepEn[idx];
    if(! sf_flag ){
    einc=0.;
    }
    */
   a_cn=id%1000;
   einc_now=einc;
  
   switch (id) {  // reset the first parameter in the Madland fit to the viola systematics, if not on the list
      case 94240:
			 // TKE_fit_param[0]=177.8; // original
			TKE_fit_param[0]=178.6;
	 break;
      case 92236:
				//original         TKE_fit_param[0]=171.45;
				TKE_fit_param[0]=171.52;
				break;
      case 92239:
				TKE_fit_param[0]=171.7;
				break;
      default:
				TKE_fit_param[0]=SystViola1966(id);
				break;
   }
   
   for (int m=1; m<numModes; m++) {
      sig_e[m]=sig0[m]+sig1[m]*einc+sig2[m]*einc*einc;
      w_e[m]=n0[m]/(1.+exp((einc-e1[m])/e2[m]));
      Abar[m]=a_cn/2.+d[m]+.5*emitted_neutrons;
   }
   Abar[0]=a_cn/2.;
   w_e[0]=2.*(1.-w_e[1]-w_e[2]);
   sig_e[0]=sig0[0];
   
}

void Yields::setRescaleTKE(double *y,int Amin,int Amax){

  double sum=0.;
  for(int a=Amin;a<=Amax;a++){
    sum+=avTKE[a]*y[a];
  }

  rescaleAverageTKE=averageTKE(einc_now)/sum;

}

double Yields::yieldA(int A){
   double a;
   if(A>=a_cn/2)
      a=(double) A;
   else {
      a=(double)(a_cn-A);
   }
   
   double sum=gaussian(a, w_e[0], Abar[0], sig_e[0]);
   for (int m=1; m<numModes; m++) {
      sum+=gaussian(a, w_e[m], Abar[m], sig_e[m]);
      sum+=gaussian(a, w_e[m], a_cn-Abar[m], sig_e[m]);
   }
   if (sum>1e-5) {
      return sum;
   }else{
      return 0.;
   }
}

double Yields::yieldATKE(int A,double TKE){
   int a;
   if(A>=a_cn/2)
      a= A;
   else {
      a=a_cn-A ;
   }
   double yatke=gaussian(TKE, yieldA(a), rescaleAverageTKE*avTKE[a], sTKE[a]);
   return yatke;
}

double Yields::yieldTKE(int A,double TKE){
   int a;
   if(A>=a_cn/2)
      a= A;
   else {
      a=a_cn-A ;
   }
   double ytke=gaussian(TKE, 1.0, rescaleAverageTKE*avTKE[a], sTKE[a]);
   if(ytke>1.e-7)
      return ytke;
   else
      return 0.;
}

double Yields::gaussian(double x, double w, double c, double s){
   double ct=w/(sqrt(PI2)*s);
   double g=exp(-.5*(x-c)*(x-c)/s/s);
   return (ct*g);
}

void Yields::setParametersAverageTKE(void){
   
   einc_ref=2.53e-8;
  
   // Madland parameters for averageTKE=f(Einc)
   // NPA 772 (2006) 113
   
   // also, load the data about avTKE(A) and sTKE[A]
   
   fill_n(TKE_fit_param, 5, 0.);
   fill_n(avTKE, 300, 140.);
   fill_n(sTKE, 300, 7.);

   switch (id_cn0){
         
      case 94240: // TODO: neutron-induced only? Ionel.
				TKE_fit_param[0]=178.6;
				// orig        TKE_fit_param[1]=-0.3489;
				TKE_fit_param[1]=-0.361;
         n_poly=1;
         for (int a=pu240_minHFFmass; a<pu240_minHFFmass+pu240_nmass; ++a) {
            avTKE[a]=pu240_avTKE[a-pu240_minHFFmass];
            sTKE[a]=pu240_sTKE[a-pu240_minHFFmass];
            avTKE[Acn0-a]=pu240_avTKE[a-pu240_minHFFmass];
            sTKE[Acn0-a]=pu240_sTKE[a-pu240_minHFFmass];
         }
         break;
      case 92236:
				TKE_fit_param[0]=171.48;
				// original         TKE_fit_param[1]=-0.1544;
				TKE_fit_param[1]=-0.1873;
				n_poly=1;
         for (int a=u236_minHFFmass; a<u236_minHFFmass+u236_nmass; ++a) {
            avTKE[a]=u236_avTKE[a-u236_minHFFmass];
            sTKE[a]= u236_sTKE[a-u236_minHFFmass];
            avTKE[Acn0-a]=u236_avTKE[a-u236_minHFFmass];
            sTKE[Acn0-a]= u236_sTKE[a-u236_minHFFmass];
         }
         break;
      case 92239:
         TKE_fit_param[0]=171.7;
         TKE_fit_param[1]=-0.2396;
         TKE_fit_param[2]=0.003434;
         n_poly=2;
         for (int a=u236_minHFFmass; a<u236_minHFFmass+u236_nmass; ++a) {
            avTKE[a]=u239_avTKE[a-pu240_minHFFmass];
            sTKE[a]=u239_sTKE[a-pu240_minHFFmass];
            avTKE[Acn0-a]=u239_avTKE[a-pu240_minHFFmass];
            sTKE[Acn0-a]=u239_sTKE[a-pu240_minHFFmass];
         }
         break;
      case 92238:
         TKE_fit_param[0]=SystViola1966();  //spontaneous fisstion, no energy dependence
         n_poly=0;
         einc_ref=0.;
         break;
      case 98252:
         TKE_fit_param[0]=SystViola1966();  //spontaneous fisstion, no energy dependence
         n_poly=0;
         einc_ref=0.;
         break;
      default:
         cout << "This compound system not implemented yet";
   }
   return;
}

double Yields::SystViola1966(void){
   
   double x= Zcn0*Zcn0 / pow((double) Acn0, 0.333333333);
   return (0.1189*x+7.3);
   
}

double Yields::SystViola1966(int id){

   int Zcn=id/1000;
   int Acn=id%1000;
   
   double x= Zcn*Zcn / pow((double) Acn, 0.333333333);
   return (0.1189*x+7.3);
   
}

// Replaces Viola systematics for Pu isotopes, and reverts to Viola for all other isotopes
double Yields::NewSyst(int id){
	
	double tke0;
	int Zcn=id/1000;
	int Acn=id%1000;
	double slope=0.;
	
	switch( id_cn0 ){
			
		case 94240:
			
			tke0=178.6;
			slope=-0.4498;
			break;
			
		default:
			
			tke0=SystViola1966(id);
			break;
			
	}
	
	double x= Zcn*Zcn / pow((double) Acn, 0.333333333);
	double x0= Zcn0*Zcn0 / pow((double) Acn0, 0.333333333);
	return (tke0+slope*(x-x0));
	
}


double Yields::averageTKE(double e){
   
   double sum=TKE_fit_param[0];
   double e1=e;
   for (int i=0; i<n_poly; ++i) {
      sum += TKE_fit_param[i+1]*e1;
      e1*=e;
   }
   return (sum);
}


double Yields::sampleTKE(int A){
   
   int a;
   if(A>=a_cn/2+a_cn%2)
      a= A;
   else {
      a=a_cn-A ;
   }
   
   return(rescaleAverageTKE*avTKE[a] + sTKE[a]*sampleGaussian() );
   
}

double Yields::get_avTKE(int A){
   
   int a;
   if(A>=a_cn/2+a_cn%2)
      a= A;
   else {
      a=a_cn-A ;
   }
   
   return(rescaleAverageTKE*avTKE[a]);
   
}

double Yields::sampleGaussian()
{
   
   static bool compute = true;
   double rand1;
   static double rand2;
   
   if(!compute)
   {
      compute=true;
      return (rand2);
   }
   
   compute = false;
   double rsq;
   
   do{
      rand1 = 2.*genrand_real3()-1.;
      rand2 = 2.*genrand_real3()-1.;
      rsq=rand1*rand1+rand2*rand2;
   }while (rsq==0. || rsq >=1.);
   
   rsq=sqrt(-(double)2.*log(rsq)/rsq);
   rand2*=rsq;
   return (rand1*rsq);
}
