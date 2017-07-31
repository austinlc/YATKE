/******************************************************************************/
/**                                                                          **/
/**                              C G M F                                     **/
/**                                                                          **/
/**--------------------------------------------------------------------------**/
/**                                                                          **/
/**                          [ LA-CC-13-063 ]                                **/
/**                                                                          **/
/**--------------------------------------------------------------------------**/
/**                                                                          **/
/** Based on CGM, T.Kawano  [ LA-CC-11-018 ]                                 **/
/**          FFD, P.Talou   [ LA-CC-10-003 ]                                 **/
/**                                                                          **/
/** Monte Carlo Hauser-Feshbach (MCHF) Code for Fission Fragments            **/
/**                                                                          **/
/**--------------------------------------------------------------------------**/
/**                                                                          **/
/** USAGE (WITHOUT MPI):                                                     **/
/**                                                                          **/
/**   ./CGMF -i $ZAIDt -e $Einc -n $nevents -y $ny                           **/
/**                                                                          **/
/** USAGE (WITH MPI): [ MPI is turned on in config.h, variable MPIRUN ]      **/
/**                                                                          **/
/**   mpirun -n $np ./CGMF -i $ZAIDt -e $Einc -f $nevents -y $ny             **/
/**                                                                          **/
/**--------------------------------------------------------------------------**/
/**                                                                          **/
/**   $ZAIDt: 1000*Z+A of target nucleus, or parent nucleus (sf)             **/
/**    $Einc: incident neutron energy in MeV (0.0 for spontaneous fission)   **/
/**      $ny: number of yields desired (produces FF initial conditions only) **/
/** $nevents: number of fission events (i.e., number Monte Carlo samplings)  **/
/**      $np: number of processors (if MPIRUN defined in config-ff.h         **/
/**                                                                          **/
/**--------------------------------------------------------------------------**/
/**                                                                          **/
/** Latest Version: 1.0.7, Sep. 1, 2015                                      **/
/**                                                                          **/
/**--------------------------------------------------------------------------**/
/** History:                                                                 **/
/**                                                                          **/
/**          1.0.7, 09/01/2015: Speed optimizations (for now)... IN PROGRESS **/
/**          1.0.6, 08/30/2015: Extension to 20 MeV, including multi-chance  **/
/**                             fission and pre-equilibrium contributions    **/
/**          1.0.5, 11/17/2014: Time coincidence window for gamma spectra    **/
/**          1.0.4, 04/22/2014: Added option to read Y(A,Z,TKE) file         **/
/**          1.0.3, ???                                                      **/
/**          1.0.2, 11/05/2013: Removed main input file option               **/
/**          1.0.1, 09/09/2013: Added reading history file option            **/
/**            1.0, 05/31/2013: Official release, LA-CC-13-063               **/
/**                                                                          **/
/**--------------------------------------------------------------------------**/
/**                                                                          **/
/**                                         T.Kawano, P.Talou and I.Stetcu   **/
/**                                                                          **/
/**                                                        kawano@lanl.gov   **/
/**                                                         talou@lanl.gov   **/
/**                                                        stetcu@lanl.gov   **/
/**                                                                          **/
/**                       T-2, Nuclear Physics Group, Theoretical Division   **/
/**                                                                          **/
/**                    kawano@lanl.gov && talou@lanl.gov && stetcu@lanl.gov  **/
/**                                                                          **/
/**                      Copyrights @ 2013, Los Alamos National Laboratory   **/
/**                                                                          **/
/******************************************************************************/

/* Copyright abstract *********************************************************
 
 The CGMF code is based on two previously codes developed at LANL: CGM
 (LA-CC-11-018) and FFD (LA-CC-10-003). It performs Monte Carlo simulations of
 the decay of excited fission fragments by emission of prompt neutrons and
 gamma rays. The Hauser-Feshbach statistical theory of compound nuclear
 reactions is used to compute the emission probabilities at each step of the
 cascade. Monte Carlo histories are recorded and analyzed.
 
 Average prompt fission neutron multiplicity (PFNM), PFNM distribution P(nu),
 average PFNM as a function of fragment characteristics <PFNM>(A,Z) can all be
 extracted from CGMF calculations. Similar quantities can also be obtained for
 prompt gamma rays. In addition, n-n, n-g, and g-g correlations can be studied
 both in energy and angle.
 
 -- References --
 
 "Monte Carlo Simulation for Particle and Gamma-Ray Emissions in Statistical
 Hauser-Feshbach Model," T.Kawano, P.Talou, M.B.Chadwick, and T.Watanabe,
 J. Nucl. Sci. Tech. 47, No.5, 462 (2010).
 
 "Advanced Monte Carlo Modeling of Prompt Fission Neutrons for Thermal and
 Fast Neutron-Induced Fission Reaction on Pu-239," P.Talou, B.Becker,
 T.Kawano, M.B.Chadwick, and Y.Danon. Phys. Rev. C83, 064612 (2011).
 
 "Monte Carlo Hauser-Feshbach Predictions of Prompt Fission Gamma Rays-
 Application to $n_{th}+^{235}$U, $n_{th}+^{239}$Pu and $^{252}$Cf (sf),"
 B.Becker, P.Talou, T.Kawano, Y.Danon, and I.Stetcu, Phys. Rev. C 87,
 014617 (2013).
 
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <netdb.h>
#include <unistd.h>

#include <netinet/in.h>    // only needed for section 8
#include <arpa/inet.h>     // only needed for section 8
#include <sys/socket.h>    // only needed for definition of AF_INET for part 5

using namespace std;

#include "mchf.h"
#include "cgm.h"
#include "global.h"
#include "terminate.h"
#include "config.h"
#include "mt19937ar.h"

#define CGM_TOPLEVEL

#ifdef HAVE_PRIVATE_ENERGY_GRID
#include  GRID_STRUCTURE_FILE
#endif

// << 1.0.6 >>
//#include "Yields.h"                      // << Any need for this here?
//#include "FissionFragments1.h"

#include "FissionFragments.h"
#include "FissionEvents.h"

static int      cgmCheckRange       (int, int, double);
static void     cgmOutputOptions    (unsigned int);
static void     cgmfHelp            (void);

static void readHistoryFile (string, int, int, double, double);
static void generateYields  (int, string, int, double, double, double);

Calculation   ctl;
Nucleus      *ncl;                 // 0: parent, 1 - MAX_COMPOUND-1: daughters etc
double       *spc[SPECTRA_OUTPUT]; // 0: gamma, 1: neutron, 2: electron, 3: neutrino
ostringstream oserr;

#define CHOICE 5
#define PARTICLE 0

bool too_many_gammas = false ;

#include "global-mcnp.h"

FissionEvents *fissionEvents;
//double **pfnInfo; // << 1.0.6 >>

string fileExt = "";   // extension to output filenames
double timeGate = 0.0; // time coincidence window (in sec)


/**********************************************************/
/*                         MAIN                           */
/**********************************************************/
int main (int argc, char *argv[])
{
  
  FissionFragments* ff;
  
  //-- Command Line Options --
  
  double Einc      = 0.0; // incident neutron energy in MeV (2.53e-8 for thermal; 0.0 for spontaneous fission)
  int    ZAIDt     = 0;   // ZAID (1000Z+A) of target nucleus (or parent nucleus in case of spontaneous fission)
  int    ZAIDf     = 0;   // 1000*Z+A of fission fragment (single fragment decay option '-s' only)
  double Eexc      = 0.0; // Excitation energy of fission fragment (single fragment decay option '-s' only)
  int    nevents   = 0;   // number of fission events per processor
  double alphaI    = 0.0; // factor to modify average initial angular momentum
  double RTin      = 0.0; // factor to modify energy sharing between fragments
      
  string historyFile = ""; // filename for histories to be read in
  string yieldsFile  = ""; // filename for input fission fragment yields given as (A, Z, TKE, Yield)

  int p;
  int ip=0; // processor index
  
  

#ifdef MPIRUN
  int np; // number of processors used
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
#endif
  if(ip == 0){
    cout << endl << string(60,'*') << endl;
    cout << "***                   C G M F - " << CGMF_VERSION << "                    ***\n";
    cout << string(60,'*') << endl << endl;
  }
  // initialize here the random number generator
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
  if (RANDOM_SEED_BY_TIME) {
    for (int i=0; i<4; i++) init[i] += (unsigned long) (time(NULL)+ip*3938*(1+i));
  }
  init_by_array(init, length);

  
  // read command line arguments
  while ((p=getopt(argc,argv,"i:e:y:h:n:u:s:a:r:t:x:b:"))!=-1) {
    switch(p){
      case 'i': ZAIDt       = atoi(optarg);   break; // ZAID target = 1000*Zt+At
      case 'e': Einc        = atof(optarg);   break; // equals fission fragment excitation energy with '-s' option
      case 'y': yieldsFile  = optarg;         break; // fission fragment yields file
      case 'n': nevents     = atoi(optarg);   break; // number of Monte Carlo events; if nf is negative, performs Yield calculations only
      case 'u': ZAIDf       = atoi(optarg);   break; // single fission fragment option
      case 'a': alphaI      = atof(optarg);   break; // factor to modify average initial angular momentum
      case 'r': historyFile = optarg;         break; // name of file containing fission histories to be read
      case 't': timeGate    = atof(optarg);   break; // time gate (sec) for counting prompt gamma rays
      case 'x': fileExt     = optarg;         break; // output filename extension
      case 'b': RTin        = atof(optarg);   break; // define RT in if selected
      case 'h': cgmfHelp();                   break;
      case ':': cerr << "ERROR     :need a value for option" << p << endl;
        cgmfHelp();                   break;
      case '?': cerr << "ERROR     :unknown option" << p << endl;
        cgmfHelp();                   break;
    }
  }

  if (timeGate==0.0) timeGate = EXPERIMENTAL_TIME_WINDOW;
  
  //-- Initializations [ hardwired CGM options for CGMF runs ] -----------------
  cgmOutputOptions(4); // o=4
  ctl.calc_montecarlo = CALC_MC;
  ctl.init_pop        = SINGLE;
  cgmAllocateMemory();
  
  
  // Special use for generating fission fragment yields for sensitivity studies
  /*  ff = new FissionFragments (98252, 0.0, &alphaI);
   string yatkefile = "YATKE1.dat";
   ff->generateYieldsForSensitivityStudies (yatkefile);
   cout << "The end!\n";
   exit(0);
   */
  
  // Reading existing Monte Carlo history file ---------------------------------
  if (historyFile!="") readHistoryFile (historyFile, nevents, ZAIDt, Einc, alphaI);
  
  // If ZAIDf==0 then it initializes the fission fragment yields Y(A,Z,TKE) for
  // the reaction ZAIDt (n,f). Otherwise, it initializes 'ff' for a single
  // fragment with a single initial excitation energy; in this case, Einc is NOT
  // an incident neutron energy, but the initial excitation energy itself.
  
  
  pfnInfo = new double *[abs(nevents)]; // << 1.0.6 >>
  
  // produces initial fission fragment yields Y(A,Z,KE,U,J,pi) only ---------------
  if (nevents<0) generateYields (nevents, "yields", ZAIDt, Einc, alphaI, RTin);
  
  if (ZAIDf==0) {
    if (yieldsFile != "") {
      ff = new FissionFragments (ZAIDt, Einc, &alphaI, yieldsFile);
    } else if (RTin != 0.0){
      ff = new FissionFragments (ZAIDt, Einc, &alphaI, RTin);
    } else {
      ff = new FissionFragments (ZAIDt, Einc, &alphaI);
    }
  } else {
    Eexc = Einc;
    ff = new FissionFragments (ZAIDf, Eexc, "single");
  }
  
  fissionFragmentType *lightFragments, *heavyFragments;
  
  lightFragments = new fissionFragmentType [nevents];
  heavyFragments = new fissionFragmentType [nevents];
  
  fissionFragmentType lf, hf;
  
  // Sample the fission fragment yields (Monte Carlo) to produce initial
  // conditions for CGM cascade calculations. The number of events is fixed by nf.
  
  fissionEvents = new FissionEvents (nevents);
  
  if (ZAIDf==0) {
    ff->generateInitialFissionFragmentHistories (lightFragments, heavyFragments, nevents);
    fissionEvents->setZAIDcn (ZAIDt, Einc); // neutron-induced fission or spontaneous fission
  } else {
    ff->generateSingleFissionFragments (lightFragments, heavyFragments, ZAIDf, Eexc, nevents);
    fissionEvents->setZAIDcn(98252, 0.0);
  }
  
  
  
  // ---------------------------------------------------------------------------
  // BEGIN LOOP OVER FISSION EVENTS
  // ---------------------------------------------------------------------------
  int fivePercent = 5*nevents/100;
  
  
  fissionEventType eventLF, eventHF;
  
  for (int ievent=0; ievent<nevents; ievent++) {
    
    if (nevents>=100){
      if (ip ==0){
      if ((ievent+1)%fivePercent==0) cout << ip << " : " << double(ievent+1)/double(nevents)*100 << "%\n";
    }
    }
    
    lf = lightFragments[ievent];
    hf = heavyFragments[ievent];
    
    fissionEvents->addFragments (lf, hf);
    
    //-- light fragment calc. --------------------------------------------------
    
    ncl[0].za.setZA (lf.Z, lf.A);
    ncl[0].max_energy = lf.U;
    
    //    specMCMain (lf.spin, lf.parity, 0.0, 0.0, 1, spc);
    specMCMain2 (lf.Z, lf.A, lf.U, lf.spin, lf.parity, 0.0, 0.0, spc, &eventLF);
    
    //-- heavy fragment calc. --------------------------------------------------
    
    ncl[0].za.setZA( hf.Z, hf.A);
    ncl[0].max_energy = hf.U;
    specMCMain2 (hf.Z, hf.A, hf.U, hf.spin, hf.parity, 0.0, 0.0, spc, &eventHF);
    //    specMCMain (hf.spin, hf.parity, 0.0, 0.0, 1, spc);
    
  }
  // ---------------------------------------------------------------------------
  // END LOOP OVER FISSION EVENTS
  // ---------------------------------------------------------------------------

  fissionEvents->writeHistories ("histories.CGMF"+fileExt, ZAIDt, Einc, nevents, alphaI);
  
  if (Einc>0.0) {
    fissionEvents->analyzeResults (ZAIDt+1); // neutron-induced
  } else {
    fissionEvents->analyzeResults (ZAIDt);   // spontaneous fission
  }
  
#ifdef MPIRUN
  if (ip==0) {
    fissionEvents->computeFinalResults();
    fissionEvents->printSummaryResults ("", ZAIDt, Einc, nevents, alphaI);
    fissionEvents->printSummaryResults ("summary.CGMF"+fileExt, ZAIDt, Einc, nevents, alphaI);
  }
  MPI_Finalize();
#else
  fissionEvents->computeFinalResults();
  fissionEvents->printSummaryResults("", ZAIDt, Einc, nevents, alphaI);
  fissionEvents->printSummaryResults ("summary.CGMF"+fileExt, ZAIDt, Einc, nevents, alphaI);
#endif
  
  fissionEvents->saveResultsToGnuplot ();
  
  // Use only to save results for GEANT simulations of DANCE (9/6/2012)
#if defined(GEANT) && defined(MPIRUN)
  fissionEvents.saveResultsGEANT ();
#endif
  
  /*** Free Allocate Memory */
  cgmDeleteAllocated();
  
  delete [] lightFragments;
  delete [] heavyFragments;
  
  if (ip==0) {
    cout << "\n*-*-* NORMAL STOP OF CGMF\n";
    cout << "*-*-* Results to be retrieved in *.CGMF files\n\n";
  }
  
  return 0;
}



/***********************************************************/
/* If nevents is given as a negative number, this routine  */
/* will be used to produce initial fission fragment yields */
/* Y(A,Z,KE,U,J,pi) and write them out.                    */
/***********************************************************/
void generateYields (int nevents, string outputFilename, int ZAIDt, double Einc, double alphaI, double RTin) {
  
  FissionFragments *ff;
  if (RTin == 0.0){
    ff = new FissionFragments (ZAIDt, Einc, &alphaI);
  } else {
    ff = new FissionFragments (ZAIDt, Einc, &alphaI, RTin);
  }
  cout << "[CGMF] Generating Fission Fragment Yields... ";
  ff->generateInitialFissionFragmentHistories ( "yields", -nevents);
  ff->checkDistributions ("yields", "yields.out");
  cout << "Done\n";
  exit(0);
  
}


/**********************************************************/
/* Read and analyze a Monte Carlo history file already    */
/* generated by CGMF.                                     */
/**********************************************************/
void readHistoryFile (string historyFile, int nevents, int ZAIDt, double Einc, double alphaI) {
  
  if (nevents==0) nevents=1000000; // set default value to 1 million events
  
  fissionEvents = new FissionEvents (nevents);
  
  fissionEvents->readHistories (historyFile, &nevents, &ZAIDt, &Einc, &alphaI);
  fissionEvents->setZAIDcn (ZAIDt, Einc);
  
  if (Einc>0.0) {
    fissionEvents->analyzeResults (ZAIDt+1); // neutron-induced fission
  } else {
    fissionEvents->analyzeResults (ZAIDt);   // spontaneous fission
  }
  
  fissionEvents->computeFinalResults();
  
  fissionEvents->printSummaryResults("", ZAIDt, Einc, nevents, alphaI);
  fissionEvents->printSummaryResults("summary.CGMF", ZAIDt, Einc, nevents, alphaI);
  
  fissionEvents->saveResultsToGnuplot ();
  
  exit(0);
  
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
  
  return(0);
}


/**********************************************************/
/*      Help                                              */
/**********************************************************/
void cgmfHelp()
{
  cout << "CGMF -i ZAnumber -e Incident Energy ... " << endl;
  cout << "option synopsis" << endl;
  cout << "    -n N : number of Monte Carlo events" << endl;
  cout << "    -y N : number of FF yields (no CGM decay is done)" << endl;
  exit(-1);
}

