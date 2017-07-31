// $Id: betaprofile.cpp 5 2011-08-09 02:08:58Z kawano $
/******************************************************************************/
/*  betaprofile.cpp                                                           */
/*        calculate strength distribution after beta-decay                    */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "cgm.h"
#include "bstrength.h"
#include "config.h"
#include "terminate.h"

static bool   betaStrengthENSDF         (Beta *);
static void   betaGaussianBroadening    (Beta *, double *, Nucleus *);
static double betaDiscreteStrength      (Beta *, double *, Nucleus *);
static double betaDiscreteRemap         (double *, double *, Nucleus *);
static double betaDiscreteScale         (double, double *, double *, Nucleus *);
static double betaRenormalizeContinuum  (double, double *, Nucleus *);
static double betaRenormalizeDiscrete   (double, double *, Nucleus *);


/*** if ENSDF sum is more than this number, we use ENSDF only */
static const double THRESHOLD_COMPLETENESS = 0.90 ;

/*** continuum beta strength zero if below this value */
static const double THRESHOLD_STRENGTH = 1.0e-10;

/*** highest discrete level energy */
static double elevhigh = 0.0;

/*** scaling factor for low lying GT discrete transitions */
static double pfactor = 1.0;


/***********************************************************/
/*      Read in GT-Strength and ENSDF, and Merge Them      */
/***********************************************************/
void betaStrengthProfile(ZAnumber *precZA, double *rc, double *rd, Nucleus *n)
{
  Beta gts,ens;

//---------------------------------------
//      Memory Allocation

  try{
    gts.br = new Branch [MAX_BETA_BRANCH];
    ens.br = new Branch [MAX_BETA_BRANCH];
  }
  catch(bad_alloc){
    cgmTerminateCode("memory allocation error in betaprofile.cpp");
  }

  for(int k=0 ; k<MAX_LEVELS     ; k++) rd[k] = 0.0;
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) rc[k] = 0.0;

//---------------------------------------
//      Read Beta Strength Files

  betaDataRead(precZA, &gts);
  betaENSDFRead(&n->za, &ens);


  if( (gts.nstate == 0) && (ens.nstate == 0) )
    cgmTerminateCode("beta-decay strength = 0");

//---------------------------------------
//      Calculate Strength Distribtion

  elevhigh = n->lev[n->ndisc-1].energy;
  double td = 0.0, tc = 0.0;
  string datasource = "";

  /*** if ENSDF is thought to be complete, use all the data */
  if(betaStrengthENSDF(&ens)){
    betaGaussianBroadening(&ens,rc,n);
    td = betaDiscreteStrength(&ens,rd,n);
    tc = betaRenormalizeContinuum(td,rc,n);
    datasource = "ENSDF";
  }
  else{
    betaGaussianBroadening(&gts,rc,n);

    /*** if ENSDF not found, use Gamow-Teller only */
    if(ens.nstate == 0){
      td = betaDiscreteRemap(rc,rd,n);
      datasource = "GT";
    }
    /*** mix ENSDF and Gamow-Teller */
    else{
      td = betaDiscreteStrength(&ens,rd,n);
      datasource = "Mixed";
    }
    if(pfactor>1.0) td = betaDiscreteScale(pfactor,rc,rd,n);
    tc = betaRenormalizeContinuum(td,rc,n);
  }


  /*** renormalize again, if no continuum */
  if( (tc == 0.0) && ((tc+td)<1.0) ){
    td = betaRenormalizeDiscrete(td,rd,n);
  }

//cout << "# "<< datasource << " : discrete fraction : " << td << endl;

//---------------------------------------
//      Free Allocated Memory

  delete [] gts.br;
  delete [] ens.br;
/*
  cout << td << "  " << tc << "  " << td + tc << endl;
  cout << elevhigh << endl;
  cout << n->excitation[0] << endl;

  for(int i=0 ; i<n->ndisc ; i++){
    cout << "D " <<  i <<"  " <<n->lev[i].energy;
    cout << "  " << n->lev[i].spin <<"  " << n->lev[i].parity;
    cout << "  " << rd[i] << endl;
  }
  for(int k=0 ; k<n->ncont ; k++){
    cout << "C " <<  k <<"  " <<n->excitation[k];
    cout << "  " <<rc[k] << endl;
  }
*/
}


/***********************************************************/
/*      Check If ENSDF Is Complete                         */
/***********************************************************/
bool betaStrengthENSDF(Beta *ensdf)
{
  bool c = false;
  /*** get a fraction of decay to discrete levels from ENSDF */
  double x = 0.0;
  for(int i=0 ; i<ensdf->nstate ; i++) x += ensdf->br[i].getR();

  if(x > THRESHOLD_COMPLETENESS){
    x = 1.0/x;
    for(int i=0 ; i<ensdf->nstate ; i++) ensdf->br[i].scaleR(x);
    c = true;
  }

  return(c);
}


/***********************************************************/
/*      Gaussian Broadening Beta Strength Data             */
/***********************************************************/
void betaGaussianBroadening(Beta *b, double *r, Nucleus *n)
{
  /*** Gaussian broadening for all beta transitions,
       generate normalized continuum final distribution */

  double g = PROFILE_BROADENING_WIDTH;
  double c = 1.0 / (sqrt(PI2)*g);

  double s = 0.0;
  /*** start k=1 to skip zero energy transition */
  for(int k=1 ; k<n->ntotal ; k++){
    r[k] = 0.0;
    for(int j=0 ; j<b->nstate ; j++){
      double d = b->br[j].getE() - n->excitation[k];
      r[k] += c * b->br[j].getR() * exp( -d*d/(g*g) );
    }

    if(r[k] < THRESHOLD_STRENGTH) r[k] = 0.0;

    s += r[k];
  }
  if(s==0.0) cgmTerminateCode("beta-decay strength zero");

  s = 1.0/s;
  for(int k=0 ; k<n->ntotal ; k++) r[k] *= s;
}


/***********************************************************/
/*      Discrete Beta Strength for Same Level Scheme       */
/***********************************************************/
double betaDiscreteStrength(Beta *b, double *r, Nucleus *n)
{
  /*** first, calculate sum of discrete transition strengths
       that are above Elevhigh (cont/disc boundary) */
  double s = 0.0;
  for(int k=0 ; k<b->nstate ;  k++){
    if(b->br[k].getE() > elevhigh) break;
    s += b->br[k].getR();
  }

  /*** find the same energy level */
  for(int k=0 ; k<b->nstate ;  k++){
    if(b->br[k].getE() > elevhigh) break;

    /*** RIPL and ENSDF level energies should be consistent */
    bool found = false;
    for(int i=0 ; i<n->ndisc ; i++){
      if( fabs(b->br[k].getE() - n->lev[i].energy)<0.001 ){
        r[i] = b->br[k].getR();
        found = true;
        break;
      }
    }
    /*** if not, force feed to the nearest level */
    if(!found){
      int i0 = 0;
      double d0 = fabs(b->br[k].getE() - n->lev[i0].energy);
      for(int i1=1 ; i1<n->ndisc ; i1++){
        double d1 = fabs(b->br[k].getE() - n->lev[i1].energy);
        if( d1<d0 ){
          d0 = d1;
          i0 = i1;
        }
      }
      r[i0] = b->br[k].getR();
    }
  }

  /*** sum strength in the discrete levels */
  double t = 0.0;
  for(int i=0 ; i<n->ndisc ; i++) t += r[i];

  return(t);
}


/***********************************************************/
/*      Renormalize Discrete Strength                      */
/***********************************************************/
double betaRenormalizeDiscrete(double t, double *r, Nucleus *n)
{
  double s = 0.0;
  if(t>0.0){
    t = 1.0/t;
    for(int i=0 ; i<n->ndisc ; i++){
      r[i] *= t;
      s += r[i];
    }
  }

  return(s);
}


/***********************************************************/
/*      Renormalize Continuum Strength                     */
/***********************************************************/
double betaRenormalizeContinuum(double t, double *r, Nucleus *n)
{
  /*** renormalize continuum part to be 1 - sum(discrete) */
  double s = 0.0;
  for(int k=0 ; k<n->ncont ; k++) s += r[k];

  if(s>0.0){
    s = (1.0-t)/s;
    for(int k=0 ; k<n->ncont ; k++) r[k] *= s;
  }

  /*** clear overlapped region */
  for(int k=n->ncont ; k<n->ntotal ; k++) r[k] = 0.0;

  s = 0.0;
  for(int k=0 ; k<n->ntotal ; k++) s += r[k];

  return(s);
}


/***********************************************************/
/*      Remap Beta Strength Data onto Discrete Levels      */
/***********************************************************/
double betaDiscreteRemap(double *rc, double *rd, Nucleus *n)
{
  /*** distribute the strength in the discrete region
       onto the nearest discerete levels,
       regardless the spin and parity.
       when no final state exists in the energy bin,
       this flux is distributed to all other level population. */

  double nolevpop = 0.0, levpop = 0.0;
  /*** scan discrete region */
  for(int k=n->ncont ; k<n->ntotal ;  k++){
    if(rc[k]==0.0) continue;

    double e1 = n->excitation[k+1];
    double e2 = n->excitation[k  ];
    int    m  = 0; // number of levels in [E1,E2]
    for(int i=1 ; i<n->ndisc ; i++){
      if( (e1 <= n->lev[i].energy) && (n->lev[i].energy < e2) ){
        m++;
        }
    }

    /*** when levels exist between E1 and E2 */
    if(m>0){
      for(int i=1 ; i<n->ndisc ; i++){
        if( (e1 <= n->lev[i].energy) && (n->lev[i].energy < e2) ){
          rd[i] += rc[k] / (double)m;
        }
      }
    }
    /*** no level case */
    else if(m==0) nolevpop += rc[k];

    levpop += rc[k];
  }

  /*** when no-level bins are found, scale-up others */
  if(levpop>nolevpop){
    double scaling = levpop/(levpop-nolevpop);
    for(int i=0 ; i<n->ndisc ; i++) rd[i] *= scaling;
  }

  double t = 0.0;
  for(int i=0 ; i<n->ndisc ; i++) t += rd[i];
  
  return(t);
}


/***********************************************************/
/*      Renormalize Beta Strength in Discrete Levels       */
/***********************************************************/
double betaDiscreteScale(double f, double *rc, double *rd, Nucleus *n)
{
  double s = 0.0;

  for(int i=0 ; i<n->ndisc ; i++) s += f*rd[i];

  /*** if scaled rd[] exceed 1.0, set zero in continuum,
       and remap the 1.0 strenght on the discrete states */
  if(s>1.0){
    for(int k=0 ; k<n->ncont ; k++) rc[k] = 0.0;
    for(int i=0 ; i<n->ndisc ; i++) rd[i] = f * rd[i] / s;
  }else{
    for(int k=0 ; k<n->ncont ; k++) rc[k] *= 1-s;
    for(int i=0 ; i<n->ndisc ; i++) rd[i] *= f;
  }

  s = 0.0;
  for(int i=0 ; i<n->ndisc ; i++) s += rd[i];

  return(s);
}


/***********************************************************/
/*      Store Beta Strength Data into Pop Arrays           */
/***********************************************************/
void betaInitialPopulation(double *rc, double *rd, Nucleus *n)
{
  /*** multiply k0-th bin population by the strength profile */
  for(int k=0 ; k<n->ncont ; k++){
    for(int j=0;j<=n->jmax;j++){
      n->pop[k][j].even *= rc[k];
      n->pop[k][j].odd  *= rc[k];
    }
  }

  /*** discrete population */
  for(int i=0 ; i<n->ndisc ; i++) n->lpop[i] = rd[i];
}

