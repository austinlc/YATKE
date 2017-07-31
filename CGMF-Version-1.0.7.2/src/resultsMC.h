/*
 *  resultsMC.h
 *  cgm
 *
 *  Created by Ionel Stetcu on 11/9/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#define AMAX_TOTAL 300
#define ZMAX_TOTAL 150
#define MAX_ENERG 1300

#define DDE 0.01

#define MAX_MULT 50 

#define HAVE_PRIVATE_ENERGY_GRID

#define NUM_SPEC_A 8

#ifdef MPIRUN
  #include <mpi.h>
#endif

#include <string>
#include <ostream>
#include <fstream>
#include <string>
#include <vector>

#ifndef __FFD_H__
	#include "ffd.h"
#endif

using namespace std;

class resultsMC {

private:
  
  int ZAIDc, Zc, Ac;
  int Asym;
  
  int ichoice ;
  int A_fragment , Z_fragment ;
  int * n , * n2 ;
  int * mult , * mult_lh[ 2 ] , * mult_mult[ MAX_MULT ] , * mult_cm , * mult_cm_lh[2] ;
  int pr[ MAX_MULT ] , pr_lh[ 2 ][ MAX_MULT ] ;
  int * y ;
  int n2_sum ;

  double **twoNeutronEnergyCorrelations;
  int ** corr_e1_e2[ 3 ] ; // this index is for the multimplicity of heavy: 0,1,2
  double * energy_av , * energy2 ;	
  double energ[ 100 ] ;
  int imult , imult_lh[ 2 ] , n_events ;
  
  int totalMultiplicity, lfMultiplicity, hfMultiplicity;

  void addN( void ) ;	
  void updateMult( double , double ) ;
  void print_phi_A( char * , int * ) ;
	
  double e_av , e2 ;
  int n_runs ;
  std::ofstream out_geant ;
  int n_part[ 2 ] , n_part2[ 2 ] ;
  double e_part[ 2 ] , e_part2[ 2 ] ;
  int lh , n2lh ;
	
  void boostlab( double e_n ) ;
  void cmToLab (double *Ecm, double *Elab, int *fragmentMass, double *fragmentKineticEnergy);
  
	
  double kinEn ;
  int a_frag ;
  double e_lab ;
	
  int ip , np ;
    
  int n_geant_sim ;	
  int * num_geant_per_sim ;
  double ** en_geant ;

  int * phi[ NUM_SPEC_A ] , * mult_za ;
  int map_a[ AMAX_TOTAL ] ;
  int parentE[ MAX_ENERG ] ;
  void print_spectra( int * , std::ofstream & ) ;
  void print_multiplicity( std::ofstream & );
  
  // MPI-specific
#ifdef MPIRUN
  MPI_Comm comm ;	
  void mpiReduceAllResults (void);
  void reduceArray (int *array, int sizeArray);
  void reduceArray (double *array, int sizeArray);
#endif
  
  double *cmSpectrum;
  double *labSpectrum;
  double **cmExclusiveSpectra;

  double **cmSpectrumVsMass;
  
  // Multiplicity Probability Distributions P(nu) + light fragment + heavy fragment P(nu)
  double *Pnu, *lfPnu, *hfPnu;

  double *nubarA;
  double *EcmA;
  double *yieldsA;
  int *countsA;
  
  double *YA, *YZ;
  double *initialSpinA;
  double *initialExcitationEnergyA;

  double *YUl, *YUh, *YUtot; // distributions of initial excitation energies in light and heavy fragments, and total
  
  // To use vectors instead... should be used everywhere though...
  vector<double> Pnu2;
  
  template<int N> int findEnergyIndex (double x0, double (&xarray)[N]);
  
public:

  void addFragmentEvent (int A, int Z, double Ui, int Ji, int Pi, int multiplicity, double *energies);
  void addFissionEvent (int Al, int Zl, double Ul, double Uh);
 
  int eventMultiplicity;
  int numberFissionEvents;
  
  void main_resultsMC( double = 0. ) ;	

#ifdef MPIRUN
  void setMPI_info( MPI_Comm , int , int ) ;
#endif
  
  resultsMC( void ) ;
  ~resultsMC( void ) ;	
  void init (int ZAIDt, double incidentEnergy);  
  void setFragment( int , int , int , double ) ;
  void printInfo( string ) ;
  void update( void ) ;
  void update_int( void ) ;
  void get_geant_fn( char * ) ;
  void geant_init( int ) ;
  double errorMC( double , double , int ) ;
  double errorMC( int , int , int ) ;
  void parentState( double , int , int ) ;
  int getFragmentA( void ) { return A_fragment ; } ;
  int getFragmentZ( void ) { return Z_fragment ; } ;
  int getProcId( void ) { return ip ; } ;

  void printAllResults (string);
  void printSummary (string);
    
  const resultsMC &operator=( const resultsMC & ) ;
  void printCorrelationResults (string filename) ;	
  int ipart ;
	
} ;

