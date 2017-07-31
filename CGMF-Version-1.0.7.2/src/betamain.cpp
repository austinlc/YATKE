// $Id: betamain.cpp 5 2011-08-09 02:08:58Z kawano $
/******************************************************************************/
/*  betamain.cpp                                                              */
/*        main betadecay gamma and neutron spectra calculation                */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "cgm.h"
#include "ripl2levels.h"
#include "masstable.h"
#include "terminate.h"


static void betaParentSpin            (int, Level *, ZAnumber *, Nucleus *);

#undef DEBUG_POPCHECK
#ifdef DEBUG_POPCHECK
static void specPopCheck           (Nucleus *);
#endif


/**********************************************************/
/*      Main Beta Spectra Calculation                     */
/**********************************************************/
void betaMain(int isoflag, double **spc)
{
  Pdata      pdt[MAX_CHANNEL];
  double     *rcnt = NULL, *rdsc = NULL;
  
//---------------------------------------
//      Memory Allocation

  try{
    rcnt     = new double [MAX_ENERGY_BIN];
    rdsc     = new double [MAX_LEVELS];
  }
  catch(bad_alloc){
    cgmTerminateCode("memory allocation error in betamain.cpp");
  }


//---------------------------------------
//      Precursor Nucleus
 
  /*** precursor Z and A */
  ZAnumber precZA(ncl[0].za.getZ()-1, ncl[0].za.getA());

  Level prec;
  betaParentSpin(isoflag,&prec,&precZA,&ncl[0]);

  if(prec.spin<0.0) cgmTerminateCode("ground state spin not assigned");


//---------------------------------------
//      Parameter Setting
  
  /*** maximal excitation of daughter nucleus, Qbeta + Isomar energy */
  ncl[0].max_energy = mass_excess(precZA.getZ(),precZA.getA())
                    - mass_excess(ncl[0].za.getZ(),ncl[0].za.getA());

  /*** how many neutrons can be emitted */
  int nemit = statTotalNeutrons(ncl[0].max_energy,&ncl[0].za);

  /*** initialize system parameters */
  statSetupInitSystem(nemit,pdt);

  /*** setup discrete levels, energy bins and level densities */    
  for(int i=0 ; i<=nemit ; i++){
    ncl[i].ndisc = riplReadDiscreteLevels(&ncl[i].za,ncl[i].lev,reassign);
    statFixDiscreteLevels(&ncl[i]);
    statSetupEnergyBin(&ncl[i]);
    statSetupLevelDensityParameter(&ncl[i],&ncl[i].ldp);
    statClearPopulation(&ncl[i]);
  }

  /*** GDR parameters */
  GDR gdr[MAX_GDR];
  statSetupGdrParameter(&ncl[0],ncl[0].ldp.a,gdr,0.0);

  /*** beta strength profile in the daughter nucleus */
  betaStrengthProfile(&precZA,rcnt,rdsc,&ncl[0]);


//---------------------------------------
//      Beta Decay Start

  /*** store initial population in the top bin */
  statSetupInitialPopulationBeta(prec.spin,prec.parity,ncl[0].ncont);

  /*** initial population in the continuum and discrete levels */
  betaInitialPopulation(rcnt,rdsc,&ncl[0]);

  /*** Hauser-Feshbach calculation */
  if(ncl[0].ncont > 0) spectra(nemit+1,pdt,spc);
  else                 specGammaCascade(ncl[0].ndisc,spc[0],&ncl[0]);


//---------------------------------------
//      Electron and Neutrino Spectra

  betaElectronSpectrum(ncl[0].max_energy,rcnt,rdsc,&ncl[0],spc);


//---------------------------------------
//      Free Allocated Memory

  delete [] rcnt;
  delete [] rdsc;
}


/***********************************************************/
/*      Find Precursor Spin and Parity                      */
/***********************************************************/
void betaParentSpin(int isoflag, Level *prec, ZAnumber *precZA, Nucleus *n)
{
  int m = riplReadDiscreteLevels(precZA,n->lev,reassign);

  *prec = n->lev[0];

  if(isoflag != 0){
    bool found = false;
    int  isomer = 0;
    for(int i=1 ; i<m; i++){
      if(n->lev[i].halflife>0.0){
        isomer++;
        if(isoflag == isomer){
          *prec = n->lev[i];
          found = true;
          break;
        }
      }
    }
    if(!found) cgmTerminateCode("isomeric state not found");
  }
}


/**********************************************************/
/*      For Monitoring Only                               */
/**********************************************************/
#ifdef DEBUG_POPCHECK
void specPopCheck(Nucleus *n)
{
  double s1 = 0.0, s2 = 0.0;

  printf("#  %03d-%03d\n",n->za.getZ(),n->za.getA());
  for(int k=0;k<n->ncont;k++){
    printf("%7.3f",n->excitation[k]);

    double s3 = 0.0;
    for(int j=0 ; j<9 ; j++) {
      s1 += n->pop[k][j].even+n->pop[k][j].odd;
      s3 += n->pop[k][j].even+n->pop[k][j].odd;
      printf("%9.2e",n->pop[k][j].even+n->pop[k][j].odd);
    }
    printf(" %9.2e",s3);
    cout << endl;
  }

  s2 = 0.0;
  for(int k=0;k<n->ndisc;k++){
    s2 += n->lpop[k];
    printf("%5d %10.6f %10.3e\n",k,n->lev[k].energy,n->lpop[k]);
  }

  printf("# SUMcont %11.4e\n",s1);
  printf("# SUMdisc %11.4e\n",s2);
  printf("# SUMall  %11.4e\n",s1+s2);
  cout << endl;
  cout << endl;
}
#endif
