// $Id: cgm.cpp 5 2011-08-09 02:08:58Z kawano $
/******************************************************************************/
/**                                                                          **/
/**   C G M  :  Cascading Gamma-ray and Multiplicity                         **/
/**                                            Version 3.5 Lysithea (2012)   **/
/**                                                              T. Kawano   **/
/**   History                                                                **/
/**   1.0  2008 May  (Io)        Hauser-Feshbach for gamma spectra           **/
/**   2.0  2008 Jul. (Europa)    Gamma cascade by Monte Carlo                **/
/**   -    2008 Dec. (Ganymede)  CoH3beta branch full Hauser-Feshbach code   **/
/**   3.0  2009 Sep. (Amalthea)  Backport version from CoH3-Callisto         **/
/**   3.1  2010 Feb. (Himalia)   Beta decay enhanced                         **/
/**   3.2  2010 Apr. (Elara)     Monte Carlo recovery version                **/
/**   3.3  2011 Jan. (Pasiphae)  Multiple neutron emission version           **/
/**   3.4  2011 Jun. (Sinope)    Monte Carlo with neutron emission version   **/
/**   3.5  2012 Feb. (Lysithea)  Neutron in Lab Frame version                **/
/**                                                                          **/
/**                                         Los Alamos National Laboratory   **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "cgm.h"
#include "global.h"
#include "terminate.h"
#include "config.h"
#include "mt19937ar.h"

#define CGM_TOPLEVEL

static int      cgmCheckRange       (int, int, double);
static void     cgmOutputOptions    (unsigned int);
static void     cgmHelp             (void);
static void     cgmAllocateMemory   (void);
static void     cgmDeleteAllocated  (void);

Calculation   ctl;
Nucleus      *ncl;                 // 0: parent, 1 - MAX_COMPOUND-1: daughters etc
double       *spc[SPECTRA_OUTPUT]; // 0: gamma, 1: neutron, 2: electron, 3: neutrino
ostringstream oserr;

#ifdef HAVE_PRIVATE_ENERGY_GRID
#include  GRID_STRUCTURE_FILE 
#endif


/**********************************************************/
/*      CGM                                               */
/**********************************************************/
int main(int argc, char *argv[])
{
//---------------------------------------
//      Command Line Options

  double        exciE = 0.0; // nuclear excitation energy
  double        initJ = -1.; // initial state spin, J
  double        targE = 0.0; // target kinetic energy, when moving
  double        spinf = 0.0; // multiplication factor for initial spin distribution (sigma)
  int           initP = 0;   // initial state parity, -1 or 1
  int           targZ = 0;   // target Z number
  int           targA = 0;   // target A number
  int           isomr = 0;   // isomeric state for beta decay
  unsigned long nsim  = 0;   // number of Monte Carlo sampling

  int           p;
  unsigned int  o=0xffff,c=0;
  while((p=getopt(argc,argv,"z:a:e:k:i:j:p:f:s:o:mbth"))!=-1){
    switch(p){
    case 'z': targZ      = atoi(optarg);   break;
    case 'a': targA      = atoi(optarg);   break;
    case 'e': exciE      = atof(optarg);   break;
    case 'k': targE      = atof(optarg);   break;
    case 'i': isomr      = atoi(optarg);   break;
    case 'j': initJ      = atof(optarg);   break;
    case 'p': initP      = atoi(optarg);   break;
    case 'f': spinf      = atof(optarg);   break;
    case 's': nsim       = atoi(optarg);   break;
    case 'o': o          = atoi(optarg);   break;
    case 'm': c          = c | CALC_MC;    break;
    case 'b': c          = c | CALC_BETA;  break;
    case 't': c          = c | CALC_TRAN;  break;
    case 'h': cgmHelp();                   break;
    case ':': cerr << "ERROR     :need a value for option" << p << endl;
              cgmHelp();                   break;
    case '?': cerr << "ERROR     :unknown optin" << p << endl;
              cgmHelp();                   break;
    }
  }


//---------------------------------------
//      Check Z, A, and E Range, spin and parity

  if(!(c & CALC_BETA)){
    if( (cgmCheckRange(targZ,targA,exciE))<0 ){
      cerr << oserr.str() << endl;  return -1;
    }

    if(initJ>=0.0 && initP==0){
      oserr << "ERROR     :parity of the state not given";
      cerr << oserr.str() << endl;  return -1;
    }
    if(initP!=0 && initJ<0.0){
      oserr << "ERROR     :spin of the state not given";
      cerr << oserr.str() << endl;  return -1;
    }
  }

   unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
   if(RANDOM_SEED_BY_TIME){
   for(int i=0 ; i<4 ; i++) init[i] += (unsigned long)time(NULL);
   }
   init_by_array(init, length);


//---------------------------------------
//      Output Items and Calc Flow  Initialize

  cgmOutputOptions(o);

  if(c != 0){
    ctl.calc_montecarlo  = c & CALC_MC;
    ctl.calc_betadecay   = c & CALC_BETA;
    ctl.calc_entrance    = c & CALC_TRAN;
  }

  if(ctl.calc_betadecay) ctl.init_pop = BETA;
  else{
    if(initJ >= 0.0) ctl.init_pop =  SINGLE;
    else{
      ctl.init_pop = (ctl.calc_entrance) ? TRANSMISSION : LEVDEN;
    }
  }  

//---------------------------------------
//      Main Part

  /***  Allocate Memory */
  cgmAllocateMemory();

  /***  Main Loop */

  if(ctl.calc_betadecay){
    ncl[0].za.setZA(targZ+1,targA);
    betaMain(isomr,spc);
  }
  else{
    ncl[0].za.setZA(targZ,targA);
    ncl[0].max_energy = exciE;

    if(ctl.calc_montecarlo) specMCMain(initJ,initP,targE,spinf,nsim,spc);
    else                    specMain(initJ,initP,spinf,spc);
  }

  /***  Print Results */
#ifdef HAVE_PRIVATE_ENERGY_GRID
  if(ctl.print_spectrum ) cgmPrintCustomGrid(NUMBER_OF_PRIVATE_GRID,custom_energy_grid,ncl[0].binwidth,spc);
#else
  if(ctl.print_spectrum ) cgmPrintSpectra(ctl.calc_betadecay,ncl[0].de,spc);
#endif

  /*** Spectrum in the LAB frame */
  if(targE > 0.0) cgmLabSpectrum(targE/targA,ncl[0].binwidth,spc);

  /*** Free Allocate Memory */
  cgmDeleteAllocated();

  return 0;
}


/**********************************************************/
/*      Setting Output Options                            */
/**********************************************************/
void  cgmOutputOptions(unsigned int p)
{
  ctl.print_spectrum     = DEFAULT_CAL & OUT_SPEC;
  ctl.print_init_pop     = DEFAULT_CAL & OUT_INIPOP;
  ctl.print_history      = DEFAULT_CAL & OUT_HIST;
  ctl.print_indivspec    = DEFAULT_CAL & OUT_INDIVSPEC;

  if(p != 0xffff){
    ctl.print_spectrum   = p & OUT_SPEC;
    ctl.print_init_pop   = p & OUT_INIPOP;
    ctl.print_history    = p & OUT_HIST;
    ctl.print_broadened  = p & OUT_BROADENED;
    ctl.print_indivspec  = p & OUT_INDIVSPEC;
  }

  if(ctl.print_broadened) ctl.print_spectrum = true;
}


/**********************************************************/
/*      Check Z/A/E Range                                 */
/**********************************************************/
int cgmCheckRange(int z, int a, double e)
{
  if(z < 6 || z>118){
    oserr << "ERROR     :Z number " << z << " out of range";
    return -1;
  }
  if(a < 12 || a>300){
    oserr << "ERROR     :A number " << a << " out of range";
    return -1;
  }
  if(e <= 0.0){
    oserr << "ERROR     :Excitation energy " << e << " negative or zero";
    return -1;
  }

  double em = 0.0;
#ifdef   HAVE_PRIVATE_ENERGY_GRID
  em = custom_energy_grid[NUMBER_OF_PRIVATE_GRID-1];
#else
  em = (MAX_ENERGY_BIN-1) * ENERGY_BIN;
#endif

  if(e >= em){
    oserr << "ERROR     :Excitation energy " << e << " too high";
    return -1;
  }

  return(0);
}


/**********************************************************/
/*      Help                                              */
/**********************************************************/
void cgmHelp()
{
  cout << "cgm -z Znumber -a Anumber -e Excitation Energy ... " << endl;
  cout << "option synopsis" << endl;
  cout << "    -j J : spin of the excited state" << endl;
  cout << "    -p P : parity of the excited state" << endl;
  cout << "           if J and P not given, gaussian dist. is assumed" << endl;
  cout << "    -f F : initial spin distribution scaled by factor F" << endl;
  cout << "    -t   : initial population calculated with neutron Tj (level density, otherwise)" << endl;
  cout << "    -k K : initial kinetic energy of compound nucleus" << endl;
  cout << "           emitted neutron energy will be converted in the lab frame" << endl;
  cout << "    -m   : do Monte Carlo (toggle)" << endl;
  cout << "    -s N : number of simulation" << endl;
  cout << "    -o N : print control; sum of the following numbers" << endl;
  cout << "           1 print energy spectrum" << endl;
  cout << "           2 initial population" << endl;
  cout << "           4 Monte Carlo history (if -m given)" << endl;
  cout << "           8 Gaussian broadening for the output spectrum" << endl;
  cout << "          16 print individual spectra from each neutron emission stage" << endl;
  exit(-1);
}


/**********************************************************/
/*     Memory Allocation                                  */
/**********************************************************/
void cgmAllocateMemory()
{
  try{
    /*** compound nucleus */
    ncl = new Nucleus [MAX_COMPOUND];

    /*** calculated results */
    for(int i=0 ; i<4 ; i++){
      spc[i] = new double[MAX_ENERGY_BIN];
    }
  }
  catch(bad_alloc){
    cerr << "ERROR     :memory allocation error";
    exit(-1);
  }
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void cgmDeleteAllocated()
{
  delete [] ncl;
  for(int i=0 ; i<SPECTRA_OUTPUT ; i++) delete [] spc[i]; 
}


/**********************************************************/
/*     Emergency Stop                                     */
/**********************************************************/
int cgmTerminateCode(string msg)
{
  /*** Release global storage */
  cgmDeleteAllocated();

  /*** Exit code */
  cerr << "ERROR     :" << msg << endl;
  exit(-1);
}

int cgmTerminateCode(string msg, int n)
{
  /*** Release global storage */
  cgmDeleteAllocated();

  /*** Exit code */
  cerr << "ERROR     :" << msg << " : " << n << endl;
  exit(-1);
}

int cgmTerminateCode(string msg, double x)
{
  /*** Release global storage */
  cgmDeleteAllocated();

  /*** Exit code */
  cerr << "ERROR     :" << msg << " : " << x << endl;
  exit(-1);
}
