// $Id: specmc.cpp 5 2011-08-09 02:08:58Z kawano $
/******************************************************************************/
/*  specmc.cpp                                                                */
/*        Monte Carlo for neutron emission and gamma-ray cascading            */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstdio>

#include <fstream>

using namespace std;

#include "cgm.h"
#include "global.h"
#include "ripl2levels.h"
#include "mchf.h"
#include "mt19937ar.h"
#include "terminate.h"
#include "config.h"

#include "global-mcnp.h"

#ifdef CGMF
#include "FissionEvents.h"
void recordEmittedParticles (void);
void recordEmittedParticles2 (fissionEventType *);
extern FissionEvents *fissionEvents;
#endif

static void    mcHistoryLoop            (MCPoint, Pdata *, double **);
static MCPoint mcFindNext               (int, int, int, int);
static int     mcGammaCascade           (int, int, double **);
static MCPoint mcSampleStart            (double, Parity **);
static MCPoint mcFixedStart             (double, Parity **);
static void    mcRestorePopulation      (MCPoint, Nucleus *, Nucleus *);
static void    mcNormalizeSpectrum      (unsigned long, double **);
static void    mcWriteHistory           (void);

static inline void   mcStoreHistory     (bool, MCPoint, MCPoint, double **);
static inline double mcCmsToLabEnergy   (double, double);

static int       np = 0;           // number of emitted gamma-rays and neutrons
static MCHistory hs[MAX_CASCADE];  // emitted energy and particle identifier
static double    cnKE = 0.0;       // compound nucleus kinetic energy in Lab

#undef DEBUG_POPCHECK
#ifdef DEBUG_POPCHECK
static void    specPopCheck        (Nucleus *); 
#endif


/**********************************************************/
/*      Monte Carlo Main Calculation                      */
/**********************************************************/

//
// targE :: compound nucleus kinetic energy (set to 0.0 if not used?)
// initJ :: initial spin in CN
// initP :: initial parity in CN
// spinf :: spin cut-off parameter scaling factor (set to 0.0 to not change anything)

void specMCMain(double initJ, int initP, double targE, double spinf, unsigned long nsim, double **spc)
{
  Transmission  tin;
  Pdata         pdt[MAX_CHANNEL];
  Parity        **pin = NULL;
	
  try{
    tin.tran = new double [3*MAX_J];
    pin      = new Parity *[MAX_ENERGY_BIN];
    for(int j=0 ; j<MAX_ENERGY_BIN ; j++) pin[j] = new Parity[MAX_J];
  }
  catch(bad_alloc){
    cgmTerminateCode("memory allocation error in specMain");
  }
  	
  //---------------------------------------
  //      Parameter Setting
	
  /*** how many neutrons can be emitted */
  int nemit = 0;
  if(INCLUDE_NEUTRON_EMISSION){
    nemit = statTotalNeutrons(ncl[0].max_energy,&ncl[0].za);
  }
  
  /*** initialize system parameters */
  statSetupInitSystem(nemit,pdt);
  
  /*** setup discrete levels, energy bins and level densities */
	riplReadDiscreteLevels(ncl, reassign, nemit);
  for(int i=0 ; i<=nemit ; i++){
//    ncl[i].ndisc = riplReadDiscreteLevels(&ncl[i].za,ncl[i].lev,reassign);
    statFixDiscreteLevels(&ncl[i]);
    statSetupEnergyBin(&ncl[i]);
    statSetupLevelDensityParameter(&ncl[i],&ncl[i].ldp);
//    statClearPopulation(&ncl[i]);
  }
  
  /*** GDR parameters */
  GDR gdr[MAX_GDR];
  statSetupGdrParameter(&ncl[0],ncl[0].ldp.a,gdr,0.0);
  
  /*** special setting for evaporation calc. */
  if(EVAPORATION_SPECTRUM){
    initJ = 0.0;
    initP = 1;
    for(int i=0 ; i<=nemit ; i++){
      ncl[i].ndisc = 1;
      ncl[i].ncont = ncl[i].ntotal;
    }
  }
  
  /*** clear spectrum array */
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    for(int i=0 ; i<SPECTRA_OUTPUT ; i++) spc[i][k] = 0.0;
  }
  	
  /*** start Monte Carlo */
  if(nsim == 0) nsim = MAX_SIMULATION;
  
  /*** decay start at a given discrete level */
  if(ncl[0].ncont == 0){
		
		/*** find nearest discrete level for given Ex */
    int nstart = statStoreLevelExcite(&ncl[0]);
    
    for(unsigned long t=0 ; t<nsim ; t++){
      np = 0;
      mcGammaCascade(0,nstart,spc);
      if(ctl.print_history) mcWriteHistory();
#ifdef CGMF
			recordEmittedParticles();
#endif
    }
  }
  /*** start with the continuum bin */
  else{
    /*** set up initial population */
    if(ctl.init_pop == TRANSMISSION){
      CrossSection cx;
      double mu = ncl[1].mass * pdt[neutron].mass / (ncl[1].mass + pdt[neutron].mass);
      tin.ecms = ncl[0].excitation[0] -  ncl[0].cdt[1].binding_energy;
      if(tin.ecms<0.0) cgmTerminateCode("closed neutron channel");
      tin.lmax = omCalc(tin.ecms, &pdt[neutron], &ncl[1].za, mu, tin.tran, &cx);
    }
    statSetupInitialPopulationSpec(initJ,initP,0,spinf,&tin);
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
			for(int j=0 ; j<MAX_J ; j++) {
				pin[k][j] = ncl[0].pop[k][j];
				
			}
    }
    
    MCPoint minit;
    double  emax[MAX_COMPOUND];
    for(int i=0 ; i<=nemit ; i++) emax[i] = ncl[i].max_energy;
		
    for(unsigned long t=0 ; t<nsim ; t++){
      /*** choose initial state */
      for(int i=0 ; i<=nemit ; i++) ncl[i].max_energy = emax[i];
      
      if(ctl.init_pop == SINGLE) minit = mcFixedStart(emax[0],pin);
      else                       minit = mcSampleStart(emax[0],pin);
      			
      /*** set CN kinetic energy if provided */
      if(targE > 0.0) cnKE = targE;
      
      /*** perform simulation */
      np = 0;
      mcHistoryLoop(minit,pdt,spc);
      if(ctl.print_history) mcWriteHistory();
#ifdef CGMF
			recordEmittedParticles();
#endif
		}
  }

  /*** restore spectra array */
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    spc[0][k] = spc[2][k];
    spc[1][k] = spc[3][k];
  }
  mcNormalizeSpectrum(nsim,spc);
  
  for(int j=0 ; j<MAX_ENERGY_BIN ; j++) delete [] pin[j];
  delete [] pin;
  delete [] tin.tran;
}

/*********************************************/

void specMCMain2(int Z, int A, double U, double initJ, int initP, double targE, double spinf, double **spc, fissionEventType *fe)
{
	
	Transmission  tin;
	Pdata         pdt[MAX_CHANNEL];
	Parity        **pin = NULL;
	
	try{
		tin.tran = new double [3*MAX_J];
		pin      = new Parity *[MAX_ENERGY_BIN];
		for(int j=0 ; j<MAX_ENERGY_BIN ; j++) pin[j] = new Parity[MAX_J];
	}
	catch(bad_alloc){
		cgmTerminateCode("memory allocation error in specMain");
	}
	
	//---------------------------------------
	//      Parameter Setting

	ncl[0].za.setZA (Z, A);
	ncl[0].max_energy = U;
	
	/*** how many neutrons can be emitted */
	int nemit = 0;
	if(INCLUDE_NEUTRON_EMISSION){
		nemit = statTotalNeutrons(ncl[0].max_energy,&ncl[0].za);
	}
	
	/*** initialize system parameters */
	statSetupInitSystem(nemit,pdt);
	
	/*** setup discrete levels, energy bins and level densities */
	riplReadDiscreteLevels(ncl, reassign, nemit);
	for(int i=0 ; i<=nemit ; i++){
		//    ncl[i].ndisc = riplReadDiscreteLevels(&ncl[i].za,ncl[i].lev,reassign);
		statFixDiscreteLevels(&ncl[i]);
		statSetupEnergyBin(&ncl[i]);
		statSetupLevelDensityParameter(&ncl[i],&ncl[i].ldp);
		//    statClearPopulation(&ncl[i]);
	}
	
	/*** GDR parameters */
	GDR gdr[MAX_GDR];
	statSetupGdrParameter(&ncl[0],ncl[0].ldp.a,gdr,0.0);
	
	/*** special setting for evaporation calc. */
	if(EVAPORATION_SPECTRUM){
		initJ = 0.0;
		initP = 1;
		for(int i=0 ; i<=nemit ; i++){
			ncl[i].ndisc = 1;
			ncl[i].ncont = ncl[i].ntotal;
		}
	}
	
	/*** clear spectrum array */
	for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
		for(int i=0 ; i<SPECTRA_OUTPUT ; i++) spc[i][k] = 0.0;
	}
	
	/*** start Monte Carlo */
	
	/*** decay start at a given discrete level */
	if(ncl[0].ncont == 0){
		
		/*** find nearest discrete level for given Ex */
		int nstart = statStoreLevelExcite(&ncl[0]);
		
		np = 0;
		mcGammaCascade(0,nstart,spc);
		if(ctl.print_history) mcWriteHistory();
#ifdef CGMF
		recordEmittedParticles2 (fe);
#endif
	}
	/*** start with the continuum bin */
	else{
		/*** set up initial population */
		if(ctl.init_pop == TRANSMISSION){
			CrossSection cx;
			double mu = ncl[1].mass * pdt[neutron].mass / (ncl[1].mass + pdt[neutron].mass);
			tin.ecms = ncl[0].excitation[0] -  ncl[0].cdt[1].binding_energy;
			if(tin.ecms<0.0) cgmTerminateCode("closed neutron channel");
			tin.lmax = omCalc(tin.ecms, &pdt[neutron], &ncl[1].za, mu, tin.tran, &cx);
		}
		statSetupInitialPopulationSpec(initJ,initP,0,spinf,&tin);
		for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
			for(int j=0 ; j<MAX_J ; j++) {
				pin[k][j] = ncl[0].pop[k][j];
				
			}
		}
		
		MCPoint minit;
		double  emax[MAX_COMPOUND];
		for(int i=0 ; i<=nemit ; i++) emax[i] = ncl[i].max_energy;
		
		/*** choose initial state */
		for(int i=0 ; i<=nemit ; i++) ncl[i].max_energy = emax[i];
		
		if(ctl.init_pop == SINGLE) minit = mcFixedStart(emax[0],pin);
		else                       minit = mcSampleStart(emax[0],pin);
			
		/*** set CN kinetic energy if provided */
		if(targE > 0.0) cnKE = targE;
			
		/*** perform simulation */
		np = 0;
		mcHistoryLoop(minit,pdt,spc);
		if(ctl.print_history) mcWriteHistory();
#ifdef CGMF
		recordEmittedParticles2 (fe);
#endif
	}
	
	/*** restore spectra array */
	for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
		spc[0][k] = spc[2][k];
		spc[1][k] = spc[3][k];
	}
	mcNormalizeSpectrum(1,spc);
	
	for(int j=0 ; j<MAX_ENERGY_BIN ; j++) delete [] pin[j];
	delete [] pin;
	delete [] tin.tran;
}



/**********************************************************/
/*      Monte Carlo Neutron and Gamma Emission            */
/**********************************************************/
void mcHistoryLoop(MCPoint minit, Pdata *pdt, double **spc)
{
  //---------------------------------------
  //      Monte Carlo Cascade Calculation
  MCPoint mcurr,mnext;
  
  mcurr = minit;
	
  while(true){
    /*** prepare energy bins and level densities */
    int c0 = mcurr.getN();
    int c1 = ncl[mcurr.getN()].cdt[neutron].next;
		
    /*** restore initial population in the top bin,
     and reconstruct bins and level densities */
    mcRestorePopulation(mcurr,&ncl[c0],&ncl[c1]);
    
    /*** calculate Hauser-Feshbach decay */
    if(EVAPORATION_SPECTRUM) specEvaporationBin(c0,pdt,spc);
    else                     spectraBinary(c0,pdt,spc);
        
    /*** move current position to the next state */
    mnext = mcFindNext(c0,c1,mcurr.getJ(),mcurr.getP());

		if(PERTURB_EXCITATON_ENERGY && !mnext.ifDiscrete()){
      int c = mnext.getN();
      int k = mnext.getK();
      double emax = (k == 0) ? mnext.getE() : (ncl[c].excitation[k-1] +  ncl[c].excitation[k])*0.5;
      double emin = (ncl[c].excitation[k+1] +  ncl[c].excitation[k])*0.5;
      if(emin < ncl[c].lev[ncl[c].ndisc-1].energy) emin = ncl[c].lev[ncl[c].ndisc-1].energy;
      
      /*** perturbation to the final state energy within its energy bin */
      double ex = emin + (emax-emin)*genrand_real2();
      mnext.setE(ex);
    }
    
    
    /*** terminate gamma-decay, if too low continuum energy */
    // Ionel and TK
    //#ifdef HAVE_PRIVATE_ENERGY_GRID
    //    if( (CONTINUUM_LOWER_CUT > 0.0) && (c0 == c1) ){
    
    // PRIOR TO AUGUST 2012
    /*    if( (CONTINUUM_LOWER_CUT > 0.0) && (-1 == c1) ){
     double elcut = ncl[c0].lev[ncl[c0].ndisc-1].energy + CONTINUUM_LOWER_CUT;
     if( (mnext.getE() < elcut) && (!mnext.ifDiscrete()) ){
     mnext.set(c0,0,0,0,0.0);
     mnext.setDiscrete();
     }
     }*/
    
    
    // ***** FROM IONEL ON SEP. 5, 2012
//    if( (CONTINUUM_LOWER_CUT > 0.0) && (-1 == c1) ){
      if( (CONTINUUM_LOWER_CUT > 0.0) && (c0 == mnext.getN()) ){
      double elcut = ncl[c0].lev[ncl[c0].ndisc-1].energy + CONTINUUM_LOWER_CUT;
      if( (mnext.getE() < elcut) && (!mnext.ifDiscrete()) ){
        int index0=-1 ;
        int ndisc = ncl[c0].ndisc ;
        double xpr = -1.;
        for( int kk2 = 0 ; kk2 <  ndisc ; kk2++ ){
          if( ncl[c0].lev[kk2].parity * mcurr.getP() > 0 ) continue ; //E1 selection rules only
          if( fabs( ncl[c0].lev[kk2].spin - mcurr.getJ() ) < 2. && 1. <= ncl[c0].lev[kk2].spin + mcurr.getJ() ){
            double xpr1 = genrand_real1() ;
            if( xpr1 > xpr ){
              xpr = xpr1 ;
              index0 = kk2;
            }
          }
        }
        
        // Ionel's correction, Jan. 2013
        /*        if( (CONTINUUM_LOWER_CUT > 0.0) && ( c0 == mnext.getN() ) ){
         double elcut = ncl[c0].lev[ncl[c0].ndisc-1].energy + CONTINUUM_LOWER_CUT;
         if( (mnext.getE() < elcut) && (!mnext.ifDiscrete()) ){
         int index0=-1 ; double delJ = 100. ;
         int ndisc = ncl[c0].ndisc ;
         double xpr = -1.;
         for( int kk2 = 0 ; kk2 <  ndisc ; kk2++ ){
         if( ncl[c0].lev[kk2].parity * mcurr.getP() > 0 ) continue ; //E1 selection rules only
         if( fabs( ncl[c0].lev[kk2].spin - mcurr.getJ() ) < 2. && 1. <= ncl[c0].lev[kk2].spin + mcurr.getJ() ){
         double xpr1 = genrand_real1() ;
         if( xpr1 > xpr ){
         xpr = xpr1 ;
         index0 = kk2;
         }
         }
         }
         */
        
        /* if such a transition does not exist, force transition to the closest-j state */
        if(index0 == -1 ){
          double delJ=100.;
          for( int kk2 = 0 ; kk2 <  ndisc ; kk2++ ){
            if( fabs( ncl[c0].lev[kk2].spin - mcurr.getJ() ) < delJ ){
              delJ = fabs( ncl[c0].lev[kk2].spin - mcurr.getJ() );
              index0 = kk2 ;
            }
          }
          //	  cout << " finished, index0 = " << index0 << endl;
        }
        mnext.set(c0, (int)ncl[c0].lev[index0].parity,(int)ncl[c0].lev[index0].spin,index0,ncl[c0].lev[index0].energy);
        mnext.setDiscrete();
      }
    }
    
    //#endif
    
    mcStoreHistory(false,mcurr,mnext,spc);
    mcurr = mnext;
    
    /*** if discrete levels reached */
    if( mnext.ifDiscrete() ) break;
  }
  
  /*** discrete level case */
  if(mnext.getK() > 0) mcGammaCascade(mnext.getN(),mnext.getK(),spc);
  
  /*** insert final state */
  int c = mcurr.getN();
  mcurr.set(c,ncl[c].lev[0].parity,(int)ncl[c].lev[0].spin,0,ncl[c].lev[0].energy);
  mcurr.setDiscrete();
  hs[np++].set(gammaray, 0.0, mcurr);
}


/**********************************************************/
/*      Check Integrant Becomes Larger than Rand          */
/**********************************************************/
MCPoint mcFindNext(int c0, int c1, int jspin, int parity)
{
  MCPoint m;
  double  r = genrand_real3();
  double  s = 0.0;
  
  /*** skip k=0, because there are initial populations */
  for(int k=1 ; k<ncl[c0].ncont ; k++){
    for(int j=0 ; j<=ncl[c0].jmax ; j++){
      s += ncl[c0].pop[k][j].even;
      if(s>=r){
        m.set(c0, 1,j,k,ncl[c0].excitation[k]);
        return(m);
      }
      s += ncl[c0].pop[k][j].odd;
      if(s>=r){
        m.set(c0,-1,j,k,ncl[c0].excitation[k]);
        return(m);
      }
    }
  }
  
  for(int k=0 ; k<ncl[c0].ndisc ; k++){
    s += ncl[c0].lpop[k];
    if(s>=r){
      m.set(c0,ncl[c0].lev[k].parity,(int)ncl[c0].lev[k].spin,k,ncl[c0].lev[k].energy);
      m.setDiscrete();
      return(m);
    }
  }
  
  if(ncl[c0].cdt[neutron].status){
    for(int k=0 ; k<ncl[c1].ncont ; k++){
      for(int j=0 ; j<=ncl[c1].jmax ; j++) {
        s += ncl[c1].pop[k][j].even;
        if(s>=r){
          m.set(c1, 1,j,k,ncl[c1].excitation[k]);
          return(m);
        }
        s += ncl[c1].pop[k][j].odd;
        if(s>=r){
          m.set(c1,-1,j,k,ncl[c1].excitation[k]);
          return(m);
        }
      }
    }
    
    for(int k=0 ; k<ncl[c1].ndisc ; k++){
      s += ncl[c1].lpop[k];
      if(s>=r){
        m.set(c1,ncl[c1].lev[k].parity,(int)ncl[c1].lev[k].spin,k,ncl[c1].lev[k].energy);
        m.setDiscrete();
        return(m);
      }
    }
  }
  
  /*** if no transition, force decay to the ground state */
  //IS: Ionel :: modified to force the transition into a state allowed by multiplicity 1 transition
  int index0=-1 ;
  int ndisc = ncl[c0].ndisc ;
  double xpr = -1.;
  for( int kk2 = 0 ; kk2 <  ndisc ; kk2++ ){
    if( ncl[c0].lev[kk2].parity * parity > 0 ) continue ; //E1 selection rules only
    if( fabs( ncl[c0].lev[kk2].spin - jspin ) < 2. && 1. <= ncl[c0].lev[kk2].spin + jspin ){
      double xpr1 = genrand_real1() ;
      if( xpr1 > xpr ){
        xpr = xpr1 ;
        index0 = kk2;
      }
    }
  }
  /* if such a transition does not exist, force transition to the closest-j state */
  if(index0 == -1 ){
    double delJ=100.;
    for( int kk2 = 0 ; kk2 <  ndisc ; kk2++ ){
      if( fabs( ncl[c0].lev[kk2].spin - jspin ) < delJ ){
        delJ = fabs( ncl[c0].lev[kk2].spin - jspin );
        index0 = kk2 ;
      }
    }
  }
  
  m.set(c0, (int)ncl[c0].lev[index0].parity,(int)ncl[c0].lev[index0].spin,index0,ncl[c0].lev[index0].energy);
  m.setDiscrete();
  //  m.set(c0,0,0,0,0.0);
  //m.setDiscrete();
  
  return(m);
}


/**********************************************************/
/*      Monte Carlo Cascading From Each Level             */
/**********************************************************/
int mcGammaCascade(int c0, int i0, double **spc)
{
  MCPoint mcurr,mnext;
  
  int m = 0;
  while(1){
    mcurr.set(c0,ncl[c0].lev[i0].parity,(int)ncl[c0].lev[i0].spin,i0,ncl[c0].lev[i0].energy);
    mcurr.setDiscrete();
		
		/** Introduced in v.1.0.5, experimental time coincidence window compared to level half-life **/
    //		if (ncl[c0].lev[i0].halflife>0.0 && EXPERIMENTAL_TIME_WINDOW>0.0) {
    //			if (genrand_real3()>1.0-exp(-EXPERIMENTAL_TIME_WINDOW/ncl[c0].lev[i0].halflife*LOG_2)) return(i0);
    //		}

		if (ncl[c0].lev[i0].halflife>0.0 && timeGate>0.0) {
			if (genrand_real3()>1.0-exp(-timeGate/ncl[c0].lev[i0].halflife*LOG_2)) return(i0);
		}

		
    double r = genrand_real3();
    double s = 0.0;
    
    for(int j=0 ; j<ncl[c0].lev[i0].ngamma ; j++){
      int i1 = ncl[c0].lev[i0].fstate[j];
      s += ncl[c0].lev[i0].branch[j];
      if(s >= r){
        /*** if Internal Concversion happens, store no gamma infor */
        bool ic = false;
        if( INCLUDE_INTERNAL_CONVERSION && (genrand_real3() > ncl[c0].lev[i0].gratio[j]) ) ic = true;
        /*** point next state */
        mnext.set(c0,ncl[c0].lev[i1].parity,(int)ncl[c0].lev[i1].spin,i1,ncl[c0].lev[i1].energy);
        mnext.setDiscrete();
        /*** save history */
        mcStoreHistory(ic,mcurr,mnext,spc);
        m++;
        i0 = i1;
        mcurr = mnext;
        break;
      }
    }
    if(i0 == 0) break;
  }
  
  return(m);
}


/**********************************************************/
/*      Sample JP from Initial Distribution               */
/**********************************************************/
MCPoint mcSampleStart(double e0, Parity **p)
{
  int     n0 = 0, k0 = 0;
  MCPoint m(n0,0,0,k0,e0);
  	
  double  x = 0.0;
  double  r = genrand_real3();
  for(int j=0 ; j<MAX_J ; j++){
    x += p[k0][j].even;
    if(x >= r){ m.set(n0, 1,j,k0,e0); return(m); }
    x += p[k0][j].odd ;
    if(x >= r){ m.set(n0,-1,j,k0,e0); return(m); }
  }
  
  m.set(n0,1,MAX_J-1,k0,e0);
  return(m);
}


/**********************************************************/
/*      Search for Starting JP State                      */
/**********************************************************/
MCPoint mcFixedStart(double e0, Parity **p)
{
  int     n0 = 0, k0 = 0;
  MCPoint m(n0,0,0,k0,e0);
  	
  /*** starting point, if fixed */
  bool found = false;
  for(int j=0 ; j<MAX_J ; j++){
    if(     p[k0][j].even > 0.0){
      m.set(n0, 1,j,k0,e0); found = true; break;
    }
    else if(p[k0][j].odd  > 0.0){
      m.set(n0,-1,j,k0,e0); found = true; break;
    }
  }
  if(!found) cgmTerminateCode("initial population zero");
  
  return(m);
}


/**********************************************************/
/*      Store Initial Population in the Top Bin           */
/**********************************************************/
void mcRestorePopulation(MCPoint m, Nucleus *n0, Nucleus *n1)
{
  n0->max_energy = m.getE();
  statSetupEnergyBin(n0);
  statSetupLevelDensity(n0,&n0->ldp);
//  statClearPopulation(n0);
  
  if(n0->cdt[neutron].status){
    n1->max_energy = n0->max_energy - n0->cdt[neutron].binding_energy;
    statSetupEnergyBin(n1);
    statSetupLevelDensity(n1,&n1->ldp);
//    statClearPopulation(n1);
  }
  
  if(m.getP() > 0) n0->pop[0][m.getJ()].even = 1.0;
  else             n0->pop[0][m.getJ()].odd  = 1.0;
}


/**********************************************************/
/*      Store Decay Gamma-Ray Energies                    */
/**********************************************************/
void mcStoreHistory(bool ic, MCPoint m0, MCPoint m1, double **spc)
{
  double ecms = 0.0;
  int k = 0;
  
  if(m0.getN() == m1.getN()){
    ecms = m0.getE() - m1.getE();
    if(ic) hs[np++].set(unknown , ecms, m0);
    else   hs[np++].set(gammaray, ecms, m0);
  }else{
    /*** if neutron, subtract separation energy */
    ecms = m0.getE() - m1.getE() - ncl[m0.getN()].cdt[neutron].binding_energy;
    hs[np++].set(neutron, ecms, m0);
  }
  if(np == MAX_CASCADE) cgmTerminateCode("too many particle decays");
  
#ifdef HAVE_PRIVATE_ENERGY_GRID
  k = specFindCustomGrid(NUMBER_OF_PRIVATE_GRID,ecms,ncl[0].binwidth);
#else
  k = specFindEnergyBin(ecms,ncl[0].de);
#endif
  
  if(k > 0){
    if(m0.getN() == m1.getN()) spc[2][k] += 1.0;
    else                       spc[3][k] += 1.0;
  }
}





/**********************************************************/
/*      Write Decay History                               */
/**********************************************************/
void mcWriteHistory()
{
  double  etot = 0.0;
  char    name[9] = {'g','n',' ',' ',' ',' ',' ',' ','e'};
  
  /*** convert neutron energy into lab-energy */
  if( cnKE > 0.0 ){
    for(int i=0 ; i<np-1 ; i++){
      if(hs[i].getIndex() == 1){
        MCPoint mp = hs[i].getMCPoint();
        double  ac = ncl[ mp.getN() ].za.getA();
        double  ec = hs[i].getEnergy();
        double  el = mcCmsToLabEnergy(ec,ac);
        hs[i].set(neutron,el,mp);
      }
    }
  }
  
  /*** the last element of hs[] contains ground state spin and parity only */
  if(EVENT_OUTPUT_FORMAT == 1){
    int j = np-1;
    for(int i=0 ; i<np-1 ; i++) if(hs[i].getIndex() == unknown) j--;
    
    cout << setw(5) << j << endl;
    cout << setprecision(5) << setiosflags(ios::scientific);
    
    /*** print out particle energies */
    for(int i=0 ; i<np-1 ; i++){
      if(hs[i].getIndex() == unknown) continue;
      cout << setw(13) << hs[i].getEnergy();
      etot += hs[i].getEnergy();
    }
    cout << endl;
  }
  
  else if(EVENT_OUTPUT_FORMAT == 2){
    cout << setprecision(5) << setiosflags(ios::scientific);
    
    /*** print out particle energies */
    for(int i=0 ; i<np-1 ; i++){
      cout << setw(4)  << name[hs[i].getIndex()];
      cout << setw(13) << hs[i].getEnergy();
      etot += hs[i].getEnergy();
    }
    cout << endl;
    
    /*** print out parent state */
    cout << "     ";
    for(int i=0 ; i<np ; i++){
      MCPoint m = hs[i].getMCPoint();
      cout << setw(3)  << m.getJ();
      cout << setw(1)  << ((m.getP() == -1) ? '-' : '+');
      cout << setw(13) << m.getE();
    }
    cout << endl;
  }
  
  else{
    //    cout << setprecision(5) << setiosflags(ios::scientific);
    //    cout << setw(5) << np-1;
    
    /*** print out particle energies */
    for(int i=0 ; i<np-1 ; i++){
      //      cout << setw(4)  << name[hs[i].getIndex()];
      //      cout << setw(13) << hs[i].getEnergy();
      etot += hs[i].getEnergy();
    }
    //  cout << " : " << etot;
    //    cout << endl;
  }
	
}






/**********************************************************/
/*      Normalize MC Spectrum                             */
/**********************************************************/
void mcNormalizeSpectrum(unsigned long n, double **spc)
{
  for(int i=0 ; i<2 ; i++){
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) spc[i][k] = spc[i][k]/n;
  }
}


/**********************************************************/
/*      Determine Lab Energy of Emitted Neutron           */
/**********************************************************/
double mcCmsToLabEnergy(double ecms, double mass)
{
  /*** sample emission angle */
  double costh = cos(PI*genrand_real1());
  double elab  = ecms + cnKE/mass + 2.0*sqrt(ecms*cnKE/mass)*costh;
  
  /*** set new recoil energy */
  double cosph = sqrt(1.0 - ecms/elab * (1.0 - costh*costh));
  cnKE = ( elab + mass*cnKE - 2.0*sqrt(elab*cnKE*mass)*cosph )/(mass-1.0);
  //cnKE = ecms/(mass-1.0)+(mass-1.0)/mass*cnKE - 2.*sqrt(ecms*cnKE/mass)*costh;
  //cout << ecms <<"  "<< elab <<"  :" << acos(costh) / PI * 180<<endl;
  
  return(elab);
}


/**********************************************************/
/*      For Monitoring Only                               */
/**********************************************************/
#ifdef DEBUG_POPCHECK
#include <iomanip>
void specPopCheck(Nucleus *n)
{
  double s1 = 0.0, s2 = 0.0;
  
  cout << "#  ";
  cout << setw(3) << setfill('0') << n->za.getZ() << '-';
  cout << setw(3) << setfill('0') << n->za.getA() << endl;
  cout << setfill(' ');
  
  for(int k=0 ; k<n->ncont ; k++){
    cout.setf(ios::fixed, ios::floatfield);
    cout << setw(7) << setprecision(3) <<  n->excitation[k];
    cout.setf(ios::scientific, ios::floatfield);
    
    double s3 = 0.0;
    for(int j=0 ; j<19 ; j++) {
      s1 += n->pop[k][j].even+n->pop[k][j].odd;
      s3 += n->pop[k][j].even+n->pop[k][j].odd;
      cout << setprecision(2) << setw(9) << n->pop[k][j].even+n->pop[k][j].odd;
    }
    cout << setprecision(2) << setw(9) << s3 << endl;
  }
  
  s2 = 0.0;
  for(int k=0;k<n->ndisc;k++){
    s2 += n->lpop[k];
    printf("%5d %10.6f %10.3e\n",k,n->lev[k].energy,n->lpop[k]);
  }
  
  cout << "# SUMcont " << s1 << endl;
  cout << "# SUMdisc " << s2 << endl;
  cout << "# SUMall  " << s1+s2 << endl;
  cout << endl;
  cout << endl;
}
#endif


/********************************************************************************
 CGMF routine to record emitted neutrons and photons from fission event.
 *******************************************************************************/
#ifdef CGMF
void recordEmittedParticles (void) {

	int index;
	MCPoint mc_curr, mc_next;
	
  double neutronEnergies [MAX_NUMBER_PARTICLES];
  double gammaEnergies [MAX_NUMBER_PARTICLES];
	int gammaTypes [MAX_NUMBER_PARTICLES];
	double icEnergies [MAX_NUMBER_PARTICLES];
	
  int neutronMultiplicity=0;
  int gammaMultiplicity=0;
	int icMultiplicity=0;
	
	std::fill_n (neutronEnergies, MAX_NUMBER_PARTICLES, 0.0);
	std::fill_n (gammaEnergies, MAX_NUMBER_PARTICLES, 0.0);
		
  MCPoint mp0 = hs[0].getMCPoint();
  int A = ncl[mp0.getN()].za.getA();
	int Z = ncl[mp0.getN()].za.getZ();
	
	float oddSpin=0.0; if (A%2) { oddSpin=0.5; }
	
	float  Ji = mp0.getJ() + oddSpin; // retrieve initial spin index
	int    Pi = mp0.getP(); // retrieve initial parity index
	double Ui = mp0.getE(); // retrieve initial excitation energy

	for (int i=0; i<np-1; i++) { // loop over emitted particles
		index = hs[i].getIndex();
		mc_curr = hs[i].getMCPoint();
		mc_next = hs[i+1].getMCPoint();
		if (index==1) { // neutrons
			neutronEnergies[neutronMultiplicity++] = hs[i].getEnergy();
		} else if (index==0) { // gammas (no IC)
			gammaTypes[gammaMultiplicity] = 1; // continuum-to-continuum
			if (mc_curr.ifDiscrete()==1) {
				gammaTypes[gammaMultiplicity] = 3; // discrete-to-discrete
			} else if (mc_next.ifDiscrete()==1) {
				gammaTypes[gammaMultiplicity] = 2; // continuum-to-discrete
			}
			gammaEnergies[gammaMultiplicity++] = hs[i].getEnergy();
		} else { // Internal Conversion
			icEnergies[icMultiplicity++] = hs[i].getEnergy();
		}
	}

	
	// make fission counter global, then no need for saveFragmentEvent!!
	// fissionEvents->lightFragments[eventCounter].emissions[0].multiplicity = neutronMultiplicity;
	// ...
		
//	fissionEvents->saveFragmentEvent (A, Z, Ui, Ji, Pi, neutronMultiplicity, neutronEnergies,
//																		gammaMultiplicity, gammaEnergies, gammaTypes, icMultiplicity, icEnergies);

}
#endif

/********************************************************************************
 CGMF routine to record emitted neutrons and photons from fission event.
 *******************************************************************************/
#ifdef CGMF
void recordEmittedParticles2 (fissionEventType *fe) {
	
	int index;
	MCPoint mc_curr, mc_next;
	
	double neutronEnergies [MAX_NUMBER_PARTICLES];
	double gammaEnergies [MAX_NUMBER_PARTICLES];
	int gammaTypes [MAX_NUMBER_PARTICLES];
	double icEnergies [MAX_NUMBER_PARTICLES];
	
	int neutronMultiplicity=0;
	int gammaMultiplicity=0;
	int icMultiplicity=0;
	
	std::fill_n (neutronEnergies, MAX_NUMBER_PARTICLES, 0.0);
	std::fill_n (gammaEnergies, MAX_NUMBER_PARTICLES, 0.0);
	
	MCPoint mp0 = hs[0].getMCPoint();
	int A = ncl[mp0.getN()].za.getA();
	int Z = ncl[mp0.getN()].za.getZ();
	
	float oddSpin=0.0; if (A%2) { oddSpin=0.5; }
	
	float  Ji = mp0.getJ() + oddSpin; // retrieve initial spin index
	int    Pi = mp0.getP(); // retrieve initial parity index
	double Ui = mp0.getE(); // retrieve initial excitation energy
	
	for (int i=0; i<np-1; i++) { // loop over emitted particles
		index = hs[i].getIndex();
		mc_curr = hs[i].getMCPoint();
		mc_next = hs[i+1].getMCPoint();
		if (index==1) { // neutrons
			neutronEnergies[neutronMultiplicity++] = hs[i].getEnergy();
		} else if (index==0) { // gammas (no IC)
			gammaTypes[gammaMultiplicity] = 1; // continuum-to-continuum
			if (mc_curr.ifDiscrete()==1) {
				gammaTypes[gammaMultiplicity] = 3; // discrete-to-discrete
			} else if (mc_next.ifDiscrete()==1) {
				gammaTypes[gammaMultiplicity] = 2; // continuum-to-discrete
			}
			gammaEnergies[gammaMultiplicity++] = hs[i].getEnergy();
		} else { // Internal Conversion
			icEnergies[icMultiplicity++] = hs[i].getEnergy();
		}
	}
	
	fe->nu  = neutronMultiplicity;
	fe->nug = gammaMultiplicity;
	
	for (int i=0; i<fe->nu; i++) { fe->neutronEnergies[i] = neutronEnergies[i]; }
	for (int i=0; i<fe->nug; i++) { fe->photonEnergies[i] = gammaEnergies[i]; }
  
//	for (int i=0; i<fe->nug; i++) { cout << gammaEnergies[i] << " "; } cout << endl;
	
	// make fission counter global, then no need for saveFragmentEvent!!
	// fissionEvents->lightFragments[eventCounter].emissions[0].multiplicity = neutronMultiplicity;
	// ...
	
	fissionEvents->saveFragmentEvent (A, Z, Ui, Ji, Pi,
                                      neutronMultiplicity, neutronEnergies,
                                      gammaMultiplicity, gammaEnergies,
                                      gammaTypes, icMultiplicity, icEnergies);
	
}
#endif

