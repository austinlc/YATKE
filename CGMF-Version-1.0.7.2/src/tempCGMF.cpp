/*
 *  tempCGMF.cpp
 *
 * Special version of CGMF driver used to study temperatures and excitation
 * energy sorting in fission fragments near scission.
 *
 * Author: P.Talou, talou@lanl.gov
 * May 2012
 *
 * USAGE:
 * 
 * time mpirun -np $NSLOTS ./mpiCGMF -s 1 -m -o 4 -n 128
 * 
 *  np: number of processors
 *   n: number of fission events per processor
 *
 */


#include "resultsMC.h"
#include "FissionFragments.h"

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

static int      cgmCheckRange       (int, int, double);
static void     cgmOutputOptions    (unsigned int);
static void     cgmHelp             (void);
static void     cgmAllocateMemory   (void);
static void     cgmDeleteAllocated  (void);

extern "C" {
  void return_bspin2c_( double * , int * , int * ) ;
}

Calculation   ctl;
Nucleus      *ncl;                 // 0: parent, 1 - MAX_COMPOUND-1: daughters etc
double       *spc[SPECTRA_OUTPUT]; // 0: gamma, 1: neutron, 2: electron, 3: neutrino
ostringstream oserr;

#ifdef HAVE_PRIVATE_ENERGY_GRID
#include  GRID_STRUCTURE_FILE 
#endif

#define CHOICE 5
#define PARTICLE 0

resultsMC   resMC[ 2 ] ;  // global variable to store the MC output 
bool too_many_gammas = false ;

int irun_times = 0 ;


/**********************************************************/
/*      CGM                                               */
/**********************************************************/
int main(int argc, char *argv[])
{
  unsigned long nsim  = 0;   // number of Monte Carlo sampling
	
  int           samples=0;
  unsigned int  o=0xffff,c=0;
  int           max_nucl, max_nucl_pp ;
  
  double spinf = 0.0 ;
  
  double bspin[300] ; int amax2 = 300 , n_s[300] ;
  double bspin2[300] , tmp[ 300 ] ; int n_s2[300] ;
  
  resultsMC resMC_store[ 2 ] ;
  
  long t_start , t_start_loop , t_stop_loop , t_io , time0 , time1 , time_ff = 0L ;
	
  int ip , np ;
	
  ifstream data_file , data_file_h ;
	
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  
  double x_clock = 1.0 / CLOCKS_PER_SEC ;
  
  if (ip==0) {
    t_start = clock();
    cout << "reading data: " << t_start * x_clock << " s" << endl ;
  }
  
  MPI_Barrier( MPI_COMM_WORLD ) ;
	
  // initialize here the random number generator
	
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
  if(RAMDOM_SEED_BY_TIME){
    for(int i=0 ; i<4 ; i++) init[i] += (unsigned long) (time(NULL) + ip * 3938 * ( 1 + i ) ) ;
  }
  irun_times = 1 ;
  init_by_array(init, length);
	
  for (int i = 0 ; i < 2 ; i++) {
    resMC[ i ].ipart = i ;
    resMC[ i ].setMPI_info( MPI_COMM_WORLD , ip , np ) ;
    resMC_store[ i ].ipart = i ;
    resMC_store[ i ].setMPI_info( MPI_COMM_WORLD , ip , np ) ;
  }
  
  // bypass user input
  
  nsim = 1; // -s 1
  o = 4; // -o 4
  c = 256; // -m
  max_nucl = 128; // -n 128
  
  //---------------------------------------
  //      Output Items and Calc Flow  Initialize
  
  cgmOutputOptions(o);
  
  if(c != 0){
    ctl.calc_montecarlo  = c & CALC_MC;
    ctl.calc_betadecay   = c & CALC_BETA;
    ctl.calc_entrance    = c & CALC_TRAN;
  }
  
  ctl.init_pop = SINGLE;

  /***  Allocate Memory */
  cgmAllocateMemory();
  
  max_nucl_pp = max_nucl ;
  
  cout << " ip = " << ip << " max= " << max_nucl_pp << endl ;
  
  int ZAIDt = 98252;
  double incidentEnergy = 0.0; // 2.53e-8; // 0.0 for spontaneous fission! 2.53e-8 for thermal energy
  
  FissionFragments ff (ZAIDt, incidentEnergy);

/*  int Zmax [NUMA];
  std::fill_n (Zmax, NUMA, 0);
  
  double b2 [NUMA];
  std::fill_n (b2, NUMA, 0.0);
  
  int jmax;
  for (int i=ff.Amin; i<ff.Amax; i++) {
  	double Ymax=0.0;
    double sum=0.0;
    for (int j=ff.Zmin; j<=ff.Zmax; j++) {
      if (ff.YZA2[j][i]!=0.0) {
        b2[i] = b2[i] + ff.beta2[j][i] * ff.YZA2[j][i];
        sum += ff.YZA2[j][i];
      }
    	if (ff.YZA2[j][i]>Ymax) { Ymax=ff.YZA2[j][i]; jmax=j; }
    }
    //    cout << i << " " << jmax << " " << b2[i] << "\n";
    if (sum!=0.0) b2[i] /= sum;
    Zmax[i]=jmax;
  }
  
  double db2 [NUMA];
  std::fill_n (db2, NUMA, 0.0);
  
  for (int i=ff.Amin+2; i<=ff.Amax-2; i++) {
  	db2[i] = (-b2[i+2]+8*b2[i+1]-8*b2[i-1]+b2[i-2])/12.;
    cout << i << " " << Zmax[i] << " " << b2[i] << " " << db2[i] << "\n";
  }
  
  exit(0);
*/
  
  ff.studyEnergySorting();
  cout << "STOP MAIN studyEnergySorting\n";
  exit(0);  
  
  fissionFragmentType *lightFragments, *heavyFragments;
  
  lightFragments = new fissionFragmentType [max_nucl_pp];
  heavyFragments = new fissionFragmentType [max_nucl_pp];
  
  ff.generateInitialFissionFragmentHistories (lightFragments, heavyFragments, max_nucl_pp);
    
  fissionFragmentType lf, hf;
    
  if (ip==0) {
   ff.generateInitialFissionFragmentHistories ( "yields2", 10000); // << Y(Z) OK!
   ff.checkDistributions ("yields2", "yields2.out");
   exit(0);
   }
   
  
  
  for (int i=0; i<2; i++) {
    resMC[i].init(ZAIDt, incidentEnergy);
    resMC[ i ].geant_init( max_nucl_pp );
    resMC_store[ i ] = resMC[ i ];
  }

  MPI_Barrier (MPI_COMM_WORLD);
  
  t_start_loop = clock() ;
  
  if (ip==0) {
    cout << "starting loop: " << t_start_loop * x_clock << " s" << endl  ;
    cout << "time to read and i/o: " << ( t_start_loop -t_start ) * x_clock << " s" << endl ;
  }
  
  MPI_Barrier( MPI_COMM_WORLD ) ;
	
  cout << "ip = " << ip << " max_nucl_pp = " << max_nucl_pp << endl ;
    
  for (samples=0; samples<max_nucl_pp; samples++) {
    
    lf = lightFragments [samples];
    hf = heavyFragments [samples];    

    // light fragment calc.
    for (int i=0; i<2; i++) resMC[i].setFragment(lf.A, lf.Z, 0, lf.KE);
    
    ncl[0].za.setZA (lf.Z, lf.A);
    ncl[0].max_energy = lf.U;
    
    specMCMain (lf.doubleSpin/2., lf.parity, 0.0, spinf, nsim, spc);
    
    for (int i=0; i<2; i++) resMC[i].update_int();
		
    // heavy fragment calc.
    for( int i = 0 ; i < 2 ; i++)	resMC[ i ].setFragment(hf.A, hf.Z, 1, hf.KE);
    
    ncl[0].za.setZA( hf.Z, hf.A);
    ncl[0].max_energy = hf.U;
    
    specMCMain (hf.doubleSpin/2., hf.parity, 0.0, spinf, nsim, spc);
    
    
    // For now, too_many_gammas ALWAYS SET TO FALSE
    if( too_many_gammas ){
      for( int i = 0 ; i < 2 ; i++ ) resMC[ i ] = resMC_store[ i ] ;
      too_many_gammas = false ;
    }
    else{
      for( int i = 0 ; i < 2 ; i++ ){		
				resMC[ i ].update() ;
				resMC_store[ i ] = resMC[ i ] ;
      }
    }
    
    // reset fission event
    resMC[1].addFissionEvent (lf.A, lf.Z, lf.U, hf.U); // neutrons
    resMC[0].addFissionEvent (lf.A, lf.Z, lf.U, hf.U); // gammas
    
  } // end loop over FF samples
  
  MPI_Barrier( MPI_COMM_WORLD ) ;
  t_stop_loop = clock() ;
  time1 = t_stop_loop - t_start_loop ;
  
  if( ip == 0 ){
    cout << "time spent in loop: " << ( t_stop_loop - t_start_loop ) * x_clock << " s" << endl ;
  }
  
  MPI_Reduce( &time1 , &time0 , 1 , MPI_LONG , MPI_MIN , 0 , MPI_COMM_WORLD ) ;
  MPI_Reduce( &time1 , &time_ff , 1 , MPI_LONG , MPI_MAX , 0 , MPI_COMM_WORLD ) ;
  
  if (ip==0) {
    cout << "min time spent in loop: " << time0 * x_clock << " s" << endl ;
    cout << "max time spent in loop: " << time_ff * x_clock << " s" << endl ;
  }
  
  
  
  /***  Print Results  */
#ifdef HAVE_PRIVATE_ENERGY_GRID
  if(ctl.print_spectrum ) cgmPrintCustomGrid(NUMBER_OF_PRIVATE_GRID,custom_energy_grid,ncl[0].binwidth,spc);
#else
  if(ctl.print_spectrum ) cgmPrintSpectra(ctl.calc_betadecay,ncl[0].de,spc);
#endif
		
  if(ctl.print_spectrum ) cgmPrintSpectra(ctl.calc_betadecay,ncl[0].de,spc);
  
  /*** Free Allocate Memory */
  cgmDeleteAllocated();
	
  if( ip == 0 ) 
    remove( "output.txt" ) ;
  
  resMC[ 0 ].printInfo( "test_g" ) ;
  resMC[ 1 ].printInfo( "test_n" ) ;
  
  //33  resMC[0].printAllResults ("gammas.CGMF.out");
  resMC[1].printAllResults ("neutrons.CGMF.out");
	resMC[0].printAllResults ("gammas.CGMF.out");
  
  for( int i = 0 ; i < amax2 ; i++ ) tmp[ i ] = bspin2[ i ] * bspin2[ i ] ;
  
  MPI_Reduce( bspin2 , bspin , amax2 , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD ) ;
  MPI_Reduce( tmp , bspin2 , amax2 , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD ) ;
  MPI_Reduce( n_s , n_s2 , amax2 , MPI_INT , MPI_SUM , 0 , MPI_COMM_WORLD ) ;
  
  if (ip==0) {
    t_io = clock() ;
    cout << "Final i/o time: " << ( t_io - t_stop_loop ) * x_clock << " s" << endl ;
  }
  
  delete [] lightFragments;
  delete [] heavyFragments; 
  
  MPI_Finalize() ;
  
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
	cout << "    -t   : initial population calculated with neutron Tj (level density, otherwise)" << endl;
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
