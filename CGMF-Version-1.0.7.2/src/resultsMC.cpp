/*
 *  resultsMC.cpp
 *  cgm
 *
 *  Created by Ionel Stetcu on 11/9/11.
 *  Copyright 2011 LANL. All rights reserved.
 *
 */

#include "resultsMC.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include "mt19937ar.h"
#include <string>
#include <sstream>
//#include "ffd.h"

using namespace std ;

const int NUMBER_OF_PRIVATE_GRID = 565; // Ionel's using 535 points instead
#define   GRID_STRUCTURE_FILE    "privategrid1.h"
#include GRID_STRUCTURE_FILE
//const int NUMBER_OF_PRIVATE_GRID = 272;
//#define   GRID_STRUCTURE_FILE    "privategrid2.h"

#define Z_TEST 50
#define A_TEST 100

static int a_spectra[ NUM_SPEC_A ] = { 88, 94, 100, 106, 136, 142, 148, 154 } ;


/*******************************************************************************
 *
 * Constructor
 *
 ******************************************************************************/
resultsMC::resultsMC()

{
  
  mult = new int[ NUMBER_OF_PRIVATE_GRID ] ;
  mult_za = new int[ NUMBER_OF_PRIVATE_GRID ] ;
  mult_cm = new int[ NUMBER_OF_PRIVATE_GRID ] ;
  
  cmSpectrum  = new double [NUMBER_OF_PRIVATE_GRID];
	labSpectrum = new double [NUMBER_OF_PRIVATE_GRID];

  cmExclusiveSpectra = new double * [MAX_MULT];
  for (int i=0; i<MAX_MULT; i++) {
    cmExclusiveSpectra[i] = new double [NUMBER_OF_PRIVATE_GRID];
    for (int j=0; j<NUMBER_OF_PRIVATE_GRID; j++) {
      cmExclusiveSpectra[i][j] = 0.0;
    }
  }

  cmSpectrumVsMass = new double * [AMAX_TOTAL];
	for (int i=0; i<AMAX_TOTAL; i++) {
    cmSpectrumVsMass [i] = new double [NUMBER_OF_PRIVATE_GRID];
  }
  
  
  for( int i = 0 ; i < MAX_MULT ; i++ ){
    if( i < 2 ){
      mult_lh[ i ] = new int[ NUMBER_OF_PRIVATE_GRID ] ;
      mult_cm_lh[ i ] =new int[ NUMBER_OF_PRIVATE_GRID ] ;
    }
    mult_mult[ i ] = new int[ NUMBER_OF_PRIVATE_GRID ] ;

  }
	
  for ( int i = 0 ; i < NUMBER_OF_PRIVATE_GRID ; i++ ) {
				
    mult[ i ] = 0 ;

    mult_cm[ i ] = 0 ;

    cmSpectrum[i]  = 0.0;
    labSpectrum[i] = 0.0;

    mult_za[ i ] = 0 ;
    for( int j = 0 ; j < MAX_MULT ; j++ ){
      if( j < 2 ){
	mult_lh[ j ][ i ] = 0 ;
	mult_cm_lh[ j ][ i ] = 0 ;
      }
      mult_mult[ j ][ i ] = 0 ;
    }
  }

  for (int i=0; i<NUM_SPEC_A; i++) {
    phi[i] = new int [NUMBER_OF_PRIVATE_GRID];
    for (int j=0; j<NUMBER_OF_PRIVATE_GRID; j++) {
      phi[i][j]=0;
    }
  }

  Pnu   = new double [MAX_MULT];
  lfPnu = new double [MAX_MULT];
  hfPnu = new double [MAX_MULT];
  			
  for (int i=0; i<MAX_MULT; i++){
 		Pnu[i]=0.0;
    lfPnu[i]=0.0;
    hfPnu[i]=0.0;
    pr[i]=0;
    for(int j=0; j< 2 ; j++)pr_lh[j][i]=0 ;
  }
			
  imult=0;
  n_runs=0;
			
  y = new int [AMAX_TOTAL];
  n = new int [AMAX_TOTAL];
	
  n2 = new int [AMAX_TOTAL];
	
  nubarA  = new double [AMAX_TOTAL];
  EcmA    = new double [AMAX_TOTAL];
  yieldsA = new double [AMAX_TOTAL];  
  countsA = new int [AMAX_TOTAL];

  YA = new double [AMAX_TOTAL];
  YZ = new double [ZMAX_TOTAL];

  initialSpinA = new double [AMAX_TOTAL];
	initialExcitationEnergyA = new double [AMAX_TOTAL];
  
  eventMultiplicity = 0;
	numberFissionEvents = 0;
  
  energy_av = new double [AMAX_TOTAL];
	
  energy2 = new double [AMAX_TOTAL];  
  
  for (int i=0; i<AMAX_TOTAL; i++) {
    y[i]=0;
    n[i]=0;
    n2[i]=0;
    
    nubarA[i]  = 0.0;
    EcmA[i]    = 0.0;
    yieldsA[i] = 0.0;
    countsA[i] = 0;
    
    YA[i] = 0.0;
    initialSpinA[i] = 0.0;
    initialExcitationEnergyA[i] = 0.0;
        
    energy_av[i]=0.;
    energy2[i]=0.;
    map_a[i]=-1;
    
    std::fill_n (cmSpectrumVsMass[i], NUMBER_OF_PRIVATE_GRID, 0.0);
        
  }

  YUl   = new double [NUMBER_OUTPUT_ENERGY_GRID];
  YUh   = new double [NUMBER_OUTPUT_ENERGY_GRID];
  YUtot = new double [NUMBER_OUTPUT_ENERGY_GRID];
  
  std::fill_n(YZ, ZMAX_TOTAL, 0.0);
  
  twoNeutronEnergyCorrelations = new double * [NUMBER_OF_PRIVATE_GRID];
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) {
  	twoNeutronEnergyCorrelations[i] = new double [NUMBER_OF_PRIVATE_GRID];
    for (int j=0; j<NUMBER_OF_PRIVATE_GRID; j++) {
      twoNeutronEnergyCorrelations[i][j] = 0.0;
    }
  }
    
  for(int i=0; i < 3; i++){
    corr_e1_e2[i]= new (int * [NUMBER_OF_PRIVATE_GRID]);
    for(int j=0 ; j < NUMBER_OF_PRIVATE_GRID ; j++ ) corr_e1_e2[i][j]= new int[NUMBER_OF_PRIVATE_GRID];
  }
  
  for (int i=0; i<NUM_SPEC_A; i++) map_a [a_spectra[i]]=i;
	
  n_geant_sim=0;
  e_av=0.;
  e2=0.;
  n2_sum=0;
	
  for (int i=0; i<2; i++) {
    e_part[i]=0.;
    n_part[i]=0;
    e_part2[i]=0.;
    n_part2[i]=0;
    imult_lh[i]=0;
  }
	
  n2lh=0;

	lfMultiplicity = 0;
	hfMultiplicity = 0;
  
  
}

/*******************************************************************************
 *
 * Destructor
 *
 ******************************************************************************/
resultsMC::~resultsMC( void )
{

  delete [] Pnu;
  delete [] lfPnu;
  delete [] hfPnu;

  delete [] cmSpectrum;
  delete [] labSpectrum;
  
  delete [] n;
  delete [] y;
  delete [] mult;
  delete [] mult_cm;
  for( int i = 0 ; i < MAX_MULT ; i++ ){
    delete [] mult_mult[i];
    if(i<2) delete [] mult_lh[i]; 
  }
  
  delete [] nubarA;
  delete [] countsA;
  delete [] yieldsA;
  delete [] EcmA;
  
  delete [] energy_av;
  delete [] energy2;
  delete [] n2;
  delete [] mult_za;
  for(int i=0;i<3;i++){
    for(int j=0; j < NUMBER_OF_PRIVATE_GRID ; j++ ) delete [] corr_e1_e2[i][j];
    delete [] corr_e1_e2[i];
  }

  return ;

  for (int i=0; i<n_geant_sim; i++) {
    delete [] en_geant[i];
  }

  for (int i=0; i<NUM_SPEC_A; i++) delete [] phi[i];

  delete [] num_geant_per_sim;
  delete [] en_geant;
	
}


/*******************************************************************************
 * init()
 *******************************************************************************
 * Set up (Z,A) of fissioning compound nucleus, and symmetric mass split (Asym).
 * To be called as soon as ZAIDt is defined.
 ******************************************************************************/
void resultsMC::init (int ZAIDt, double incidentEnergy)
{
	if (incidentEnergy==0.0) { // Spontaneous fission
    ZAIDc = ZAIDt;
  } else { // neutron-induced fission
    ZAIDc = ZAIDt+1;
  }
  Ac = ZAIDc%1000;
  Zc = int(ZAIDc/1000);
  Asym = Ac/2;
}

/*******************************************************************************
 *
 ******************************************************************************/
inline void resultsMC::addN (void)
{
  n[A_fragment]++;
  n2_sum++;	
}

/*******************************************************************************
 *
 ******************************************************************************/
void resultsMC::setFragment (int A, int Z, int light_heavy, double ke) 
{
	
  A_fragment=A;
  Z_fragment=Z;
  lh=light_heavy; // light (0) or heavy (1)
  kinEn=ke;
  a_frag=A;
	
}


/*******************************************************************************
 * print out final results
 ******************************************************************************/
void resultsMC::printInfo (string filename)
{
	
  string fn;
  stringstream s;
	
  int *itemp, *itemp2, *ytot;
  double *temp1, *temp2;
		
  ofstream out_file;
	
  double sum=0.;
  int n_part1=0;

  double *tmp1, *tmp2;

  itemp = new int[MAX_ENERG];

#ifdef MPIRUN
  MPI_Reduce (parentE, itemp, MAX_ENERG, MPI_INT, MPI_SUM, 0, comm);
#endif
  
  if (ip==0) { // processor rank 0
    
    s << filename << "_eParent.dat"; // to identify Ui when low energy gammas are emitted
    fn=s.str();
    out_file.open (&filename[0], ios::out);		
    sum=0.;
    for (int i=0; i<MAX_ENERG; i++) sum += (double) itemp[i];
    for (int i=0; i<MAX_ENERG; i++) out_file << i * DDE << " " << itemp[ i ] / sum << endl;
    out_file.close();
    s.str("");

  }

  delete [] itemp;

  // Spectrum in lab
  
  tmp1 = new double [NUMBER_OF_PRIVATE_GRID];
  tmp2 = new double [NUMBER_OF_PRIVATE_GRID];
  temp1 = new double [NUMBER_OF_PRIVATE_GRID];
  temp2 = new double [NUMBER_OF_PRIVATE_GRID];

  sum=0;
	
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) sum += (double) mult[i];

  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) {
      temp1[i] = mult[i] / sum / ( custom_energy_grid[i+1] - custom_energy_grid[i] );
      temp2[i] = temp1[i]*temp1[i];
  }

#ifdef MPIRUN
  MPI_Reduce (temp1, tmp1, NUMBER_OF_PRIVATE_GRID, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce (temp2, tmp2, NUMBER_OF_PRIVATE_GRID, MPI_DOUBLE, MPI_SUM, 0, comm);
#endif
  
  if (ip==0) {

    // spectrum
    s << filename << "_spectr.dat";
    fn=s.str();
    cout << endl << filename << endl;
    out_file.open (fn.c_str(), ios::out);
    out_file << setprecision(4) << setiosflags(ios::scientific);
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) {
      out_file << custom_energy_grid[i] << " " << tmp1[i] / np << "  " << errorMC (tmp2[i], tmp1 [i], np) << endl;			
    }
    out_file.close();
  }

  delete [] temp1;
  delete [] temp2;
  delete [] tmp1;
  delete [] tmp2;

  itemp = new int [NUMBER_OF_PRIVATE_GRID];

  print_phi_A (&filename[0], itemp);
	
  delete [] itemp;
	
  temp1 = new double [MAX_MULT];
  temp2 = new double [MAX_MULT];

  sum=0.;
  for (int i=0; i<MAX_MULT; i++) sum += (double) pr[i];

  for (int i=0; i<MAX_MULT; i++) {
    temp1[i] = pr[i]/sum;
    temp2[i] = temp1[i]*temp1[i];
  }

  tmp1 = new double [MAX_MULT];
  tmp2 = new double [MAX_MULT];	

#ifdef MPIRUN
  MPI_Reduce (temp1, tmp1, MAX_MULT, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce (temp2, tmp2, MAX_MULT, MPI_DOUBLE, MPI_SUM, 0, comm);
#endif
  
  if (ip==0) {

    // P(nu)
    s.str("");
    s << filename << "_multPr.dat";
    fn=s.str();
    out_file.open (fn.c_str(), ios::out);  	
        
    for (int i=0; i<MAX_MULT; i++) { // probablility to have a certain multiplicity
      out_file << i << " " << tmp1[i] / np << "  " << errorMC (tmp2[i], tmp1[i], np) << endl;		
    }
	
    out_file.close() ;
		
  }
	
  delete [] tmp1;
  delete [] tmp2;
  delete [] temp1;
  delete [] temp2;
	
  itemp = new int [AMAX_TOTAL];
  itemp2 = new int [AMAX_TOTAL];
  ytot = new int [AMAX_TOTAL];
  temp1 = new double [AMAX_TOTAL];
  temp2 = new double [AMAX_TOTAL];
	
#ifdef MPIRUN
  MPI_Reduce (n , itemp, AMAX_TOTAL, MPI_INT, MPI_SUM, 0, comm);
  MPI_Reduce (n2 , itemp2, AMAX_TOTAL, MPI_INT, MPI_SUM, 0, comm);
  MPI_Reduce (energy_av, temp1, AMAX_TOTAL, MPI_DOUBLE , MPI_SUM, 0, comm);
  MPI_Reduce (energy2, temp2, AMAX_TOTAL, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce (y, ytot, AMAX_TOTAL, MPI_INT, MPI_SUM, 0, comm);
#endif
  
  if (ip==0) {

		// <nu>(A)    
    s.str( "" );
    s << filename << "_avp.dat";
    fn=s.str();
    out_file.open (fn.c_str(), ios::out);
    
    n2_sum=0;
    out_file << setprecision(4) << setiosflags(ios::scientific);
    
    for (int k=0; k<AMAX_TOTAL; k++) { // average number of particles as a function of A
      n_part1 += itemp[k];
      n2_sum += itemp2[k];
      if (ytot[k]>0) out_file << k << " " << itemp[k] / ( (double) ytot[k] ) << " " << errorMC (itemp2[k], itemp[k], ytot[k]) << endl;
    }
	
		// <Ecm>(A)
    out_file.close();
    s.str("");
    s << filename << "_aven_A.dat";
    fn=s.str();
    out_file.open (fn.c_str(), ios::out);
		
    for (int k=0; k<AMAX_TOTAL; k++) { // average energy as a function of A
      if (itemp[k]>0) out_file << k << " " << temp1[k]/itemp[k] << " " << errorMC (temp2[k], temp1[k], itemp[k]) << endl;
    }
	
    out_file.close();
		
  }
	
  delete [] itemp;
  delete [] itemp2;
  delete [] temp1;
  delete [] temp2;
  delete [] ytot;
	
  double e_av_tot, e2_tot;
  int n_events_tot, n_runs_tot, n2_sum_tot;
	
#ifdef MPIRUN
  MPI_Reduce (&e_av, &e_av_tot, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce (&e2, &e2_tot, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	
  MPI_Reduce (&n_events, &n_events_tot, 1, MPI_INT, MPI_SUM, 0, comm);
  MPI_Reduce (&n_runs, &n_runs_tot, 1, MPI_INT, MPI_SUM, 0, comm);
  MPI_Reduce (&n2_sum, &n2_sum_tot, 1, MPI_INT, MPI_SUM, 0, comm);
#endif
  
  if (ip==0) {

    out_file.open ("output.txt", ios::out | ios::app);

    if (ipart==0) out_file << "****************" << endl << " gamma " << endl << "****************" << endl;

    if (ipart==1) out_file << endl << "****************" << endl << " neutrons " << endl << "****************" << endl;

    cout << "average energy: " << e_av_tot / ( (double ) n_events_tot ) << " +/- " << errorMC( e2_tot , e_av_tot , n_events_tot ) <<  endl ;	
    cout << "average number of events: " << n_part1 / ( (double ) n_runs_tot ) << " +/- " << errorMC( n2_sum_tot , n_part1 , n_events_tot ) << endl ;
    cout << "number of runs " << n_runs_tot << endl ;
    cout << "number of particles: " << n_part1 <<endl ;

    out_file << "average energy: " << e_av_tot / ( (double ) n_events_tot ) << " +/- " << errorMC( e2_tot , e_av_tot , n_events_tot ) <<  endl ;
    out_file << "average number of events: " << n_part1 / ( (double ) n_runs_tot ) << " +/- " << errorMC( n2_sum_tot , n_part1 , n_events_tot ) << endl ;
    out_file << "number of runs " << n_runs_tot << endl ;	
    out_file << "number of particles: " << n_part1 <<endl ;

    out_file.close() ;
	
    /*		cout << "Light fragment energy average: " << e_part[ 0 ] / n_runs << " +/- " << errorMC( e_part2[ 0 ] , e_part[ 0 ] , n_runs ) << endl ;
		cout << "Heavy fragment energy average: " << e_part[ 1 ] / n_runs << " +/- " << errorMC( e_part2[ 1 ] , e_part[ 1 ] , n_runs ) << endl ;
		cout << "Light fragment average # of particles: " << ( ( double ) n_part[ 0 ] ) / n_runs << " +/- " << errorMC( n_part2[ 0 ] , n_part[ 0 ] , n_runs  ) << endl ;
		cout << "Heavy fragment average # of particles: " << ( ( double ) n_part[ 1 ] ) / n_runs << "+/- " << errorMC( n_part2[ 1 ] , n_part[ 1 ] , n_runs  ) << endl ;
    */
		
  }
	
  // geant data

  /*
	
  s.str("") ;

  s << filename << "_geant.dat" ;

  fn = s.str() ;

  if( ip == 0 )
		
    remove( fn.c_str() ) ;
	
  for( int n = 0 ; n < np ; n++ ){
		
    if( n == ip ){
			
      out_geant.open( fn.c_str() , ios::out | ios::app ) ;
      out_geant << setprecision(6) << setiosflags(ios::fixed);
	    
      for ( int i = 0 ; i < n_geant_sim ; i++ ) {
	out_geant << num_geant_per_sim[ i ] << endl ;
	for( int j = 0 ; j < num_geant_per_sim[ i ] ; j++ )
	  out_geant << en_geant[ i ][ j ] << "   " ;
	out_geant << endl ;
      }
			
      out_geant.close() ;
			
    }
		
    MPI_Barrier( comm ) ;
		
  }

  */
	
}

/*******************************************************************************
 * Updating
 *
 * en: energy in center-of-mass
 *
 ******************************************************************************/
void resultsMC::main_resultsMC (double en)
{

  // counting number of events per A		
  n[A_fragment]++;
  
  n2_sum++;	// for error calc

  boostlab (en) ; // boost cm energy to lab energy
//  if (ipart==1) cmToLab (&en, &e_lab, &a_frag, &kinEn); // only used for neutrons
  
  updateMult (e_lab, en); // computes spectrum in lab
  energ[imult]=e_lab; // saved for GEANT file
  energy_av[A_fragment] += en; // <Ecm>(A)
  
  energy2 [A_fragment] += en*en; // for error calc.
  imult++; // increment number of emitted particles
  e_av += e_lab; // to calculate <Elab>
  e2 += e_lab*e_lab; // for error calc.
  n_events++; // total number of particles

  // not used for now
  e_part[lh] += en;
  e_part2[lh] += en*en;
  n_part[lh]++;
  n2lh++;
  imult_lh[lh]++;

  totalMultiplicity++;
  if (lh==0) {
    lfMultiplicity++;
  } else {
    hfMultiplicity++;
  }
  
}

/*******************************************************************************
 * Computes spectrum in lab and cm frames
 ******************************************************************************/
void resultsMC::updateMult (double en, double e_cm)
{
	
  for (int k=1; k<NUMBER_OF_PRIVATE_GRID; k++) {
    if (custom_energy_grid[k]>en) {
      mult[k-1]++;
      labSpectrum[k-1]++;
      break;
    }
  }

  for ( int k = 1 ; k < NUMBER_OF_PRIVATE_GRID ; k++ ) {    
    if ( custom_energy_grid[ k ] > en ) {
      mult_lh[ lh ][ k - 1 ]++ ;
      break ;			
    }	
  }

  for ( int k = 1 ; k < NUMBER_OF_PRIVATE_GRID ; k++ ) {    
    if ( custom_energy_grid[ k ] > e_cm ) {
      mult_cm_lh[ lh ][ k - 1 ]++ ;
      break ;			
    }	
  }

  for ( int k = 1 ; k < NUMBER_OF_PRIVATE_GRID ; k++ ) {
    if ( custom_energy_grid[ k ] > e_cm ) {
      mult_cm[ k - 1 ]++ ;
      break ;		
    }	
  }

  int m=map_a[A_fragment];

  if (m>=0) {
    for (int k=1; k<NUMBER_OF_PRIVATE_GRID; k++) {
      if (custom_energy_grid[k]>en) {
      	phi[m][k-1]++;
	break;
      }
    }
  }

  if (Z_TEST==Z_fragment && A_TEST==A_fragment) {
    for (int k=0; k<NUMBER_OF_PRIVATE_GRID; k++) {
      if (e_cm>custom_energy_grid[k]) {
        mult_za[k]++;
        break;
      }
    }
  }

  return;

}

/*******************************************************************************
 * addFragmentEvent
 *------------------------------------------------------------------------------
 * Stores the characteristics of the decay of a given fragment in a fission
 * event.
 * Input: (A,Z) of fragment, followed by initial excitation energy, spin and
 * parity, then followed by the number of particles emitted (multiplicity) and
 * their energies.
 ******************************************************************************/
// PT new, May 1st, 2012
void resultsMC::addFragmentEvent (int A, int Z, double Ui, int Ji, int Pi, int multiplicity, double *energies)
{
  int k;

  
  yieldsA[A]++;


  for (int i=0; i<multiplicity; i++) {
  	EcmA[A] += energies[i];
    k = findEnergyIndex(energies[i], custom_energy_grid);
//    cmSpectrum[k]++;
//    cmSpectrumVsMass[A][k]++;
  }

  if (multiplicity>=2) {
  	double Ecm1 = energies[0];
    double Ecm2 = energies[1];
    int k1 = findEnergyIndex (Ecm1, custom_energy_grid);
		int k2 = findEnergyIndex (Ecm2, custom_energy_grid);
//		twoNeutronEnergyCorrelations[k1][k2]++;    
  }

  nubarA[A] += multiplicity;
	if (A<Asym) {
    lfPnu[multiplicity]++;
  } else {
    hfPnu[multiplicity]++;
  }

  initialExcitationEnergyA[A] += Ui;
  initialSpinA[A] += Ji;
  
}

/*******************************************************************************
 * Finds the index corresponding to the value 'x0' in the array 'xarray'.
 ******************************************************************************/
template <int N> int resultsMC::findEnergyIndex (double x0, double (&xarray)[N])
{
	int k0;
  int numberElements = N;
  for (int k=1; k<numberElements; k++) {
  	if (xarray[k]>x0) {
      k0=k-1;
      break;
    }
  }
  return k0;
}

/*******************************************************************************
 * addFissionEvent
 *------------------------------------------------------------------------------
 * Stores the characteristics of a fission event, which depends on both
 * fragments, so it needs to be called AFTER the decay of both fragments.
 *
 * Input: light fragment mass (Al), charge (Zl), total kinetic energy 
 * (TKE [MeV]).
 ******************************************************************************/
void resultsMC::addFissionEvent (int Al, int Zl, double Ul, double Uh)
{ 
  
  Pnu[eventMultiplicity]++; // increment multiplicity probability
  eventMultiplicity=0;      // reset particle multiplicity for next fission event
  numberFissionEvents++;    // increment total number of fission events

  YA[Al]++;
  YA[Ac-Al]++;
  YZ[Zl]++;
  YZ[Zc-Zl]++;

  int kl = findEnergyIndex (Ul, outputEnergyGrid);    YUl[kl]++;
  int kh = findEnergyIndex (Uh, outputEnergyGrid);    YUh[kh]++;
  int k  = findEnergyIndex (Ul+Uh, outputEnergyGrid); YUtot[k]++;
  
//  cout << kl << " " << kh << " " << k << "\n";

}


/*******************************************************************************
 // Intermediate updating results MC for 1st fragment only 
 // (always light fragment in cgmFF.cpp algorithm!)
 ******************************************************************************/
void resultsMC::update_int ()
{
  y[A_fragment]++; // Y(A)
  countsA[A_fragment]++;
  n2[A_fragment]+=n2_sum*n2_sum; // N^2 (A) for error calc
  n2_sum=0;
  n_part2[lh]+=n2lh*n2lh; // not used
  n2lh=0;  
}

/*******************************************************************************
 * Updating resMC[].
 ******************************************************************************/
void resultsMC::update() // for both neutrons and gammas
{

  pr[imult]++;      // P(nu)
  y[A_fragment]++;  // mass yields	
  countsA[A_fragment]++;
  n_runs++;         // total number of counts
  n2[A_fragment] += n2_sum*n2_sum; // N^2 (A) for error calc
  n2_sum=0;

  for( int i=0 ; i < imult ; i++ ){ //spectrum function of multiplicity 
    for (int k=1; k<NUMBER_OF_PRIVATE_GRID; k++) {
      if (custom_energy_grid[k]>energ[i]) {
				mult_mult[imult][k-1]++;
        cmExclusiveSpectra[imult][k-1]++;
				break;
      }
    }
  }

  int idx[2] ;
  if(imult==2){ //correlations between the neutron energies

    for(int i=0; i<imult_lh[0];i++){
      for (int k=1; k<NUMBER_OF_PRIVATE_GRID; k++) {
	if (custom_energy_grid[k]>energ[i]) {
	  idx[i]=k;
	  break;
	}
      }
    }
    for(int i=0; i<imult_lh[1];i++){
      for (int k=1; k<NUMBER_OF_PRIVATE_GRID; k++) {
	if (custom_energy_grid[k]>energ[i+imult_lh[0]]) {
	  idx[i+imult_lh[0]]=k;
	  break;
	}
      }
    }
    corr_e1_e2[imult_lh[1]][idx[0]][idx[1]]++;      
  }

  imult=0; // multiplicity count
  n_part2[lh] += n2lh*n2lh; // not used
  n2lh=0;
  
  for(int i=0 ; i<2 ; i++) {
    pr_lh[i][imult_lh[i]]++ ;
    imult_lh[i]=0;
  }

  return;

  // not used anymore
  if (imult>0) {
    num_geant_per_sim[n_geant_sim] = imult;
    en_geant[n_geant_sim] = new double [imult];
    for (int i=0; i<imult; i++) {
      en_geant[n_geant_sim][i] = energ[i];
    }
    n_geant_sim++;
  }
	
}

/*******************************************************************************
 * ???
 ******************************************************************************/
void resultsMC::get_geant_fn (char* filename)
{
  out_geant.open (filename, ios::out);	
  out_geant << setprecision(6) << setiosflags(ios::fixed);
  return;
}

/*******************************************************************************
 * Computes sqrt [ 1/N * ( <x^2> - <x>^2 ) ] for error calculation.
 ******************************************************************************/
double resultsMC::errorMC (double x2, double x, int nsamp)
{
  double dx;
  dx = (x2/nsamp-pow(x/nsamp,2.))/nsamp;
  return (sqrt(dx));
}

/*******************************************************************************
 * Same as above for integers.
 ******************************************************************************/
double resultsMC::errorMC (int x2, int x, int nsamp)
{
  double dx, xsamp=(double) nsamp, xx2=(double) x2, xx=(double) x;
  dx = ((xx2)/xsamp-(xx*xx)/xsamp/xsamp)/xsamp;
  return (sqrt(dx));
}

/*******************************************************************************
 * Transforms a center-of-mass neutron energy into laboratory energy. Also,
 * computes the kinematic recoil, and remove one neutron from the nucleus with
 * mass A.
 * 
 * Follows Eq.(1) in J.Terrell, Phys. Rev. 113, 527 (1958).
 ******************************************************************************/
void resultsMC::boostlab( double e_n )
{
	
  if( ipart == 0 ) { // gammas
    e_lab = e_n;
    return;
  }
	
  double cos_th, a1;
  a1 = (double) a_frag;
  
  cos_th = 1. - 2. * genrand_real1();	
    
  e_lab = e_n + kinEn / a1 + 2. * sqrt( e_n * kinEn / a1 ) * cos_th ;
  kinEn = e_n / ( a1 - 1. )+ ( a1 - 1. ) * kinEn / a1 - 2. * sqrt( e_n * kinEn / a1 ) * cos_th ;
  a_frag-- ;
  
  return ;
}

/*******************************************************************************
 * Similar to boostlab() but specifies I/O arguments.
 * SHOULD only be used for neutrons.
 *
 * Follows Eq.(1) in J.Terrell, Phys. Rev. 113, 527 (1959).
 ******************************************************************************/
void resultsMC::cmToLab (double *Ecm, double *Elab, int *fragmentMass, double *fragmentKineticEnergy)
{
  double cosTheta;
  cosTheta = 1. - 2. * genrand_real1();

  double Af = (double) *fragmentMass;
    
  *Elab = *Ecm + *fragmentKineticEnergy / Af + 2.0 * sqrt (*Ecm * *fragmentKineticEnergy / Af ) * cosTheta ;
  *fragmentKineticEnergy = *Ecm / (Af-1.) + (Af-1.) * *fragmentKineticEnergy / Af - 2.0 * sqrt( *Ecm * *fragmentKineticEnergy / Af ) * cosTheta;
 
  *fragmentMass--;

  return ;
}

/*******************************************************************************
 * MPI instruction to set current processor and total number of processors used
 * in the present calculations.
 ******************************************************************************/
#ifdef MPIRUN
void resultsMC::setMPI_info( MPI_Comm commun , int ip1 , int np1) {	
  ip = ip1 ;
  np = np1 ;
  comm = commun ;
}
#endif

/*******************************************************************************
 * ???
 ******************************************************************************/
void resultsMC::geant_init( int events_perproc ) {
  num_geant_per_sim = new int [ events_perproc ] ;
  en_geant = new ( double * [ events_perproc ] ) ;
  return ;
}

/*******************************************************************************
 * Saves spectrum for a given fragment?
 ******************************************************************************/
inline void resultsMC::print_phi_A( char * filename , int * multipl ) {

  ofstream out_file ; 

  int itemp[ NUMBER_OF_PRIVATE_GRID ] ;

  for (int i=0; i<NUM_SPEC_A; i++) {
    for (int j=0 ; j<NUMBER_OF_PRIVATE_GRID; j++)
      itemp[j] = phi[i][j];

#ifdef MPIRUN
    MPI_Reduce( itemp , multipl , NUMBER_OF_PRIVATE_GRID , MPI_INT , MPI_SUM, 0 , comm ) ;
#endif
    
    if (ip == 0 ) {

      stringstream s ;

      s << filename << "_A" << a_spectra[ i ] << "_spec.dat" ;

      string fn = s.str();
      out_file.open( fn.c_str() , ios::out ) ;

      out_file << setprecision(4) << setiosflags(ios::scientific);

      double sum=0.;
      for (int j=0; j<NUMBER_OF_PRIVATE_GRID-1;j++)
				sum += (double) multipl[j];

      for (int j=0; j<NUMBER_OF_PRIVATE_GRID-1; j++)
        out_file << custom_energy_grid[ j ] << " " << multipl[ j ] / sum / ( custom_energy_grid[ j + 1 ] - custom_energy_grid[ j ] ) << endl ;

      out_file.close() ;

    }
  }

  return;
}

/*******************************************************************************
 * Stores the initial excitation energy?
 * >> jspin and nn NOT USED?
 ******************************************************************************/
void resultsMC::parentState (double e, int jspin, int nn ){

  for( int i = 0 ; i < MAX_ENERG ; i++ )
    if( ( i + 1 ) * DDE > e ){
      parentE[ i ]++ ;
      //      cout << "energy=" << e << " " << i << endl ;
      break ;
    }
  
}

/*******************************************************************************
 * Definition of "=" operator between resultsMC types.
 * >> Needed?
 ******************************************************************************/
const resultsMC &resultsMC::operator=( const resultsMC &right ){

  if( &right != this ){

    A_fragment = right.A_fragment ;
    Z_fragment = right.Z_fragment ;
    n2_sum = right.n2_sum ;
    e2 = right.e2 ;
    n_runs = right.n_runs ;
    ichoice = right.ichoice ;
    e_av = right.e_av ;
    n_events = right.n_events ;
    imult = right.imult ;
    n2lh = right.n2lh ;
    lh = right.lh ;

    for (int i=0; i<2; i++) {

      n_part[i]  = right.n_part[i];
      n_part2[i] = right.n_part2[i];
      e_part[i]  = right.e_part[i];
      e_part2[i] = right.e_part2[i];

    }

    for (int i=0; i<MAX_ENERG; i++) parentE[i] = right.parentE[i];

    for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) mult[i] = right.mult[i];
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) labSpectrum[i] = right.labSpectrum[i];

    for (int i=0; i<MAX_MULT; i++) pr[i] = right.pr[i];
    for (int i=0; i<MAX_MULT; i++) Pnu[i] = right.Pnu[i];

    for (int i=0; i<AMAX_TOTAL; i++) {
      n[i]         = right.n[i];
      n2[i]        = right.n2[i];
      energy_av[i] = right.energy_av[i];
      energy2[i]   = right.energy2[i];
      y[i]         = right.y[i];
      
      nubarA[i]    = right.nubarA[i];
      EcmA[i]      = right.EcmA[i];
      countsA[i]   = right.countsA[i];
      yieldsA[i]   = right.yieldsA[i];
      
    }

    /*

    if( n_geant_sim > right.n_geant_sim ){

      for( int n_gnt = right.n_geant_sim ; n_gnt < n_geant_sim )

	delete [] en_geant[ n_geant_sim ] ;

      n_geant_sim-- ;

    }else { 

      for( int n_gnt = n_geant_sim ; n_gnt < right.n_geant_sim ; n_gnt++ ){

	num_geant_per_sim[ n_gnt ] = right.num_geant_per_sim

    */

  }

  return *this ;

}


/*******************************************************************************
 * mpiReduceAllResults
 *------------------------------------------------------------------------------
 * "Reduce" all results from different processors and place the final results
 * in processor of rank 0.
 ******************************************************************************/
#ifdef MPIRUN
void resultsMC::mpiReduceAllResults (void)
{
  double sum;
  
  if (ip==0) cout << "Reducing all results...\n";
  
  // Transform histograms into spectra
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) {
  	labSpectrum[i] /= (custom_energy_grid[i+1]-custom_energy_grid[i]);
  	cmSpectrum[i]  /= (custom_energy_grid[i+1]-custom_energy_grid[i]);
    for (int j=1; j<MAX_MULT; j++) cmExclusiveSpectra[j][i] /= (custom_energy_grid[i+1]-custom_energy_grid[i]);
    for (int j=1; j<AMAX_TOTAL; j++) cmSpectrumVsMass[j][i] /= (custom_energy_grid[i+1]-custom_energy_grid[i]);
  }
  
  // reduce spectra
  reduceArray (cmSpectrum, NUMBER_OF_PRIVATE_GRID);
  reduceArray (labSpectrum, NUMBER_OF_PRIVATE_GRID);
  for (int j=1; j<10; j++) reduceArray (cmExclusiveSpectra[j], NUMBER_OF_PRIVATE_GRID);
	for (int j=1; j<AMAX_TOTAL; j++) reduceArray (cmSpectrumVsMass[j], NUMBER_OF_PRIVATE_GRID);
  
  // renormalize spectra
  sum=0.0;
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) sum += cmSpectrum[i]*(custom_energy_grid[i+1]-custom_energy_grid[i]);
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) cmSpectrum[i] /= sum;
  
  sum=0.0;
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) sum += labSpectrum[i]*(custom_energy_grid[i+1]-custom_energy_grid[i]);
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) labSpectrum[i] /= sum;
  for (int j=1; j<MAX_MULT; j++) {
    sum=0.0;
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) sum += cmExclusiveSpectra[j][i]*(custom_energy_grid[i+1]-custom_energy_grid[i]);
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) cmExclusiveSpectra[j][i] /= sum;
	}
  for (int j=1; j<AMAX_TOTAL; j++) {
    sum=0.0;
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) sum += cmSpectrumVsMass[j][i]*(custom_energy_grid[i+1]-custom_energy_grid[i]);
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) cmSpectrumVsMass[j][i] /= sum;
	}
  
  // reduce initial excitation energy distributions Y(U_l) and Y(U_h)
  reduceArray (YUl,   NUMBER_OUTPUT_ENERGY_GRID);
  reduceArray (YUh,   NUMBER_OUTPUT_ENERGY_GRID);
  reduceArray (YUtot, NUMBER_OUTPUT_ENERGY_GRID);
  
  // reduce multiplicity distributions
  for (int i=0; i<MAX_MULT; i++) {
    Pnu[i]   /= numberFissionEvents;
    lfPnu[i] /= numberFissionEvents;
    hfPnu[i] /= numberFissionEvents;
  }
  reduceArray (Pnu, MAX_MULT);
	reduceArray (lfPnu, MAX_MULT);
	reduceArray (hfPnu, MAX_MULT);

  double *dummy;
  dummy = new double [AMAX_TOTAL];
  std::fill_n (dummy, AMAX_TOTAL, 0.0);

  // <Ecm>(A)
  for (int i=0; i<AMAX_TOTAL; i++) {
    dummy[i] = EcmA[i];
    if (yieldsA[i]!=0) dummy[i] /= nubarA[i];
  }
  MPI_Reduce (dummy, EcmA, AMAX_TOTAL, MPI_DOUBLE , MPI_SUM, 0, comm);
  for (int i=0; i<AMAX_TOTAL; i++) EcmA[i] /= np;

  // <nu>(A)
  std::fill_n (dummy, AMAX_TOTAL, 0.0);
  for (int i=0; i<AMAX_TOTAL; i++) {
    dummy[i] = nubarA[i];
    if (yieldsA[i]!=0) dummy[i] /= yieldsA[i];
  }
  MPI_Reduce (dummy, nubarA, AMAX_TOTAL, MPI_DOUBLE, MPI_SUM, 0, comm);
  for (int i=0; i<AMAX_TOTAL; i++) nubarA[i] /= np;

  // <Ui>(A)
  std::fill_n (dummy, AMAX_TOTAL, 0.0);
  for (int i=0; i<AMAX_TOTAL; i++) {
    if (yieldsA[i]!=0) dummy[i] = initialExcitationEnergyA[i] / yieldsA[i];
  }
  MPI_Reduce (dummy, initialExcitationEnergyA, AMAX_TOTAL, MPI_DOUBLE , MPI_SUM, 0, comm);
	for (int i=0; i<AMAX_TOTAL; i++) initialExcitationEnergyA[i] /= np;
  
  // <Ji>(A)
  std::fill_n (dummy, AMAX_TOTAL, 0.0);
  for (int i=0; i<AMAX_TOTAL; i++) {
    if (yieldsA[i]!=0) dummy[i] = initialSpinA[i] / yieldsA[i];
  }
  MPI_Reduce (dummy, initialSpinA, AMAX_TOTAL, MPI_DOUBLE , MPI_SUM, 0, comm);
	for (int i=0; i<AMAX_TOTAL; i++) initialSpinA[i] /= np;
    
  // Y(A)
  sum = 0.0;
  std::fill_n (dummy, AMAX_TOTAL, 0.0);
  for (int i=0; i<AMAX_TOTAL; i++) sum += YA[i];
  for (int i=0; i<AMAX_TOTAL; i++) dummy[i] = YA[i] / sum;
  MPI_Reduce (dummy, YA, AMAX_TOTAL, MPI_DOUBLE , MPI_SUM, 0, comm);
  for (int i=0; i<AMAX_TOTAL; i++) YA[i] = YA[i] / np * 100.0;

  // Y(Z)
  delete [] dummy;
  dummy = new double [ZMAX_TOTAL];
  sum = 0.0;
  for (int i=0; i<ZMAX_TOTAL; i++) sum += YZ[i];
  for (int i=0; i<ZMAX_TOTAL; i++) dummy[i] = YZ[i] / sum;
  MPI_Reduce (dummy, YZ, ZMAX_TOTAL, MPI_DOUBLE , MPI_SUM, 0, comm);
  for (int i=0; i<ZMAX_TOTAL; i++) YZ[i] = YZ[i] / np * 100.0;
  
  // Mass Yields Y(A), from counts to percent -- CAUTION: this should be done after <nu>(A), <Ecm>(A), ...
  sum = 0.0;
  for (int i=0; i<AMAX_TOTAL; i++) sum += yieldsA[i];
  for (int i=0; i<AMAX_TOTAL; i++) yieldsA[i] /= sum;
  

}


/*******************************************************************************
 * reduceArray
 *------------------------------------------------------------------------------
 * MPI reduce an array to processor of rank 0.
 ******************************************************************************/
void resultsMC::reduceArray (double *array, int sizeArray)
{
	int i;
  double sum=0.0;
  double *dummy;
  dummy = new double [sizeArray];
  for (i=0; i<sizeArray; i++) sum += array[i];
  for (i=0; i<sizeArray; i++) dummy[i] = array[i]/sum;
  MPI_Reduce (dummy, array, sizeArray, MPI_DOUBLE, MPI_SUM, 0, comm);
  for (i=0; i<sizeArray; i++) array[i] /= np;
  return;
}

/*******************************************************************************
 * Same as above for an array of integers.
 ******************************************************************************/
void resultsMC::reduceArray (int *array, int sizeArray)
{
	double *newArray;
  newArray = new double [sizeArray];
  for (int i=0; i<sizeArray; i++) newArray[i] = (double) array[i];
  reduceArray(newArray, sizeArray);
  return;
}

#endif // endif MPIRUN


/*******************************************************************************
 * printAllResults
 *------------------------------------------------------------------------------
 * Saves all results in one output file (filename) to be used by GNUPLOT for
 * instance.
 ******************************************************************************/
void resultsMC::printAllResults (string filename)
{
  
  ofstream out;
  
  //*** TMP PT
  // Transform histograms into spectra
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) {
  	labSpectrum[i] /= (custom_energy_grid[i+1]-custom_energy_grid[i]);
  	cmSpectrum[i]  /= (custom_energy_grid[i+1]-custom_energy_grid[i]);
    for (int j=1; j<MAX_MULT; j++) cmExclusiveSpectra[j][i] /= (custom_energy_grid[i+1]-custom_energy_grid[i]);
    for (int j=1; j<AMAX_TOTAL; j++) cmSpectrumVsMass[j][i] /= (custom_energy_grid[i+1]-custom_energy_grid[i]);
  }
    
  // renormalize spectra
  double sum=0.0;
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) sum += cmSpectrum[i]*(custom_energy_grid[i+1]-custom_energy_grid[i]);
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) cmSpectrum[i] /= sum;
  
  sum=0.0;
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) sum += labSpectrum[i]*(custom_energy_grid[i+1]-custom_energy_grid[i]);
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) labSpectrum[i] /= sum;
  
  //*** END TMP PT
  
  
  // Testing new approach to reduce all results
#ifdef MPIRUN
  mpiReduceAllResults();
#endif  
  if (ip==0) { // only on main processor
    
    ofstream out;
    out.open(&filename[0]);
    out << setprecision(4) << setiosflags(ios::scientific);
		out << "#\n# CGMF Results\n#\n";
    
    int gnuplotIndex=-1;

    //-- Laboratory PFNS
    out << "# [gnuplot #" << ++gnuplotIndex << "] lab. spectrum\n\n";
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) 
      out << custom_energy_grid[i] << " " << labSpectrum[i] << "\n";
      
    //-- cm PFNS and multiplicity-dependent exclusive spectra
    out << "\n# [gnuplot #" << ++gnuplotIndex << "] c.m. spectrum\n";
    out << "# Energy    Total       nu=1      nu=2      nu=3       nu=4      nu=5\n";
    out << "# (MeV)     (1/MeV)\n\n";
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) 
      out << custom_energy_grid[i] << " " << cmSpectrum[i] << " " <<
      cmExclusiveSpectra[1][i] << " " << cmExclusiveSpectra[2][i] << " " << 
      cmExclusiveSpectra[3][i] << " " << cmExclusiveSpectra[4][i] << " " <<
      cmExclusiveSpectra[5][i] << "\n";

    //-- Multiplicity distributions P(nu), and same for light and heavy fragments resp.
		out << "\n# [gnuplot #" << ++gnuplotIndex << "] P(nu) P_LF(nu) P_HF(nu)\n\n";
    for (int i=0; i<20; i++) 
      out << i << " " << Pnu[i] << " " << lfPnu[i] << " " << hfPnu[i] << "\n";
    
    //-- Fragment mass-dependent results: mass Yields Y(A), <nu>(A), <Ecm>(A), etc.
    out << "\n# [gnuplot #" << ++gnuplotIndex << "] Y(A)   <nu>(A)   <Ecm>(A)  <Ui>(A)   <Ji>(A)\n\n";
  	for (int i=0; i<AMAX_TOTAL; i++) 
      out << i << " " << yieldsA[i] << " " << nubarA[i] << " " 
      	<< EcmA[i] << " " << initialExcitationEnergyA[i] << " " 
      	<< initialSpinA[i] << "\n";

    //-- Z-dependent results: Y(Z), ...
    out << "\n#  [gnuplot #" << ++gnuplotIndex << "] Z   Y(Z) \n\n";
    for (int i=0; i<ZMAX_TOTAL; i++)
      out << i << " " << YZ[i] << "\n";
    
    //-- Initial Excitation Energy Distributions Y(Ul) and Y(Uh)
    out << "\n# [gnuplot #" << ++gnuplotIndex << "] Y(Ul), Y(Uh), Y(Utot)\n\n";
    for (int i=0; i<NUMBER_OUTPUT_ENERGY_GRID; i++) {
    	out << outputEnergyGrid[i] << " " << YUl[i] << " " << YUh[i] << " " << YUtot[i] << "\n";
    }
    
    out.close();

  }

  //-- (e1,e2) two-neutron energy correlations (in center-of-mass)
  double tmp [NUMBER_OF_PRIVATE_GRID];
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) {
#ifdef MPIRUN
    MPI_Reduce (twoNeutronEnergyCorrelations[i], tmp, NUMBER_OF_PRIVATE_GRID, MPI_DOUBLE, MPI_SUM, 0, comm);
#endif
    if (ip==0) {
      for (int j=0; j<NUMBER_OF_PRIVATE_GRID; j++) {
        twoNeutronEnergyCorrelations[i][j] = tmp[j];
      }
    }
  }

  double **correlationMatrix;
  correlationMatrix = new double * [NUMBER_OF_PRIVATE_GRID];
  if (ip==0) {
  	for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) {
      correlationMatrix[i] = new double [NUMBER_OF_PRIVATE_GRID];
      for (int j=0; j<NUMBER_OF_PRIVATE_GRID; j++) {
        if (twoNeutronEnergyCorrelations[i][i] != 0.0 && twoNeutronEnergyCorrelations[j][j] != 0.0) {
	        correlationMatrix[i][j] = twoNeutronEnergyCorrelations[i][j] / 
  		      sqrt ( twoNeutronEnergyCorrelations[i][i] * twoNeutronEnergyCorrelations[j][j] );
        } else {
          correlationMatrix[i][j] = 0.0;
        }
      }
    }
  }
  

  // n-n Energy Correlations
  if (ip==0) {
    string correlationFilename;
    correlationFilename = filename + ".corr";
    out.open(&correlationFilename[0]);
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) {
      for (int j=0; j<NUMBER_OF_PRIVATE_GRID; j++) {
        out << .5 * ( custom_energy_grid[i] + custom_energy_grid[i+1] ) << " " 
        		<< .5 * ( custom_energy_grid[j] + custom_energy_grid[j+1] ) << " " 
        		<< correlationMatrix[i][j] << " " << twoNeutronEnergyCorrelations[i][j] << "\n";
      }
      out << endl ;
    }
    out << endl ;
    out.close();
  }


}

/*******************************************************************************
 * Print a general summary of CGMF results.
 ******************************************************************************/
void resultsMC::printSummary (string outputFilename) {

  ofstream out;
  
  if (ip==0) {
    
    out.open (&outputFilename[0], ios::out | ios::app);
    
    if (ipart==0) out << "****************" << endl << " gamma " << endl << "****************" << endl;
    
    if (ipart==1) out << endl << "****************" << endl << " neutrons " << endl << "****************" << endl;
    
/*    cout << "average energy: " << e_av_tot / ( (double ) n_events_tot ) << " +/- " << errorMC( e2_tot , e_av_tot , n_events_tot ) <<  endl ;	
    cout << "average number of events: " << n_part1 / ( (double ) n_runs_tot ) << " +/- " << errorMC( n2_sum_tot , n_part1 , n_events_tot ) << endl ;
    cout << "number of runs " << n_runs_tot << endl ;
    cout << "number of particles: " << n_part1 <<endl ;
    
    out << "average energy: " << e_av_tot / ( (double ) n_events_tot ) << " +/- " << errorMC( e2_tot , e_av_tot , n_events_tot ) <<  endl ;
    out << "average number of events: " << n_part1 / ( (double ) n_runs_tot ) << " +/- " << errorMC( n2_sum_tot , n_part1 , n_events_tot ) << endl ;
    out << "number of runs " << n_runs_tot << endl ;	
    out << "number of particles: " << n_part1 <<endl ;
*/    
    out.close();
    		
  }  
    
}



/*******************************************************************************
 * MPI reduce and print spectrum array.
 ******************************************************************************/
void resultsMC::print_spectra( int * inputArray , ofstream & out ){

  double * tmp1 , * tmp2 , * temp1 , * temp2 ;

  tmp1  = new double [NUMBER_OF_PRIVATE_GRID];
  tmp2  = new double [NUMBER_OF_PRIVATE_GRID];
  temp1 = new double [NUMBER_OF_PRIVATE_GRID];
  temp2 = new double [NUMBER_OF_PRIVATE_GRID];
    
  double sum=0.0;
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) sum += (double) inputArray[i];
  
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) {
    temp1[i] = inputArray[i] / sum / ( custom_energy_grid[i+1] - custom_energy_grid[i] );
    temp2[i] = temp1[i]*temp1[i];
  }
  
#ifdef MPIRUN
  MPI_Reduce (temp1, tmp1, NUMBER_OF_PRIVATE_GRID, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce (temp2, tmp2, NUMBER_OF_PRIVATE_GRID, MPI_DOUBLE, MPI_SUM, 0, comm);
#endif
  
  if (ip==0) { // processor rank 0
    out << setprecision(4) << setiosflags(ios::scientific);
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID-1; i++) {
      out << custom_energy_grid[i] << " " << tmp1[i] / np << "  " << errorMC (tmp2[i], tmp1 [i], np) << endl;			
    }

  }

  out << endl;

  delete [] tmp1 ; delete [] tmp2 ; delete [] temp1 ; delete [] temp2 ;

}

/*******************************************************************************
 * MPI reduce and print multiplicity.
 ******************************************************************************/
void resultsMC::print_multiplicity( ofstream &out ){

  double * temp1, *temp2, *tmp1, *tmp2, *tmp1L , *tmp2L , *tmp1H, *tmp2H;

  temp1 = new double [MAX_MULT];
  temp2 = new double [MAX_MULT];
  
  double sum=0.;
  for (int i=0; i<MAX_MULT; i++) sum += (double) pr[i];
  
  for (int i=0; i<MAX_MULT; i++) {
    temp1[i] = pr[i]/sum;
    temp2[i] = temp1[i]*temp1[i];
  }
  
  tmp1 = new double [MAX_MULT];
  tmp2 = new double [MAX_MULT];

#ifdef MPIRUN
  MPI_Reduce (temp1, tmp1, MAX_MULT, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce (temp2, tmp2, MAX_MULT, MPI_DOUBLE, MPI_SUM, 0, comm);
#endif
  
  sum = 0.;
  for (int i=0; i<MAX_MULT; i++) sum += (double) pr_lh[0][i];
  
  for (int i=0; i<MAX_MULT; i++) {
    temp1[i] = pr_lh[0][i]/sum;
    temp2[i] = temp1[i]*temp1[i];
  }
  
  tmp1L = new double [MAX_MULT];
  tmp2L = new double [MAX_MULT];	

#ifdef MPIRUN
  MPI_Reduce (temp1, tmp1L, MAX_MULT, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce (temp2, tmp2L, MAX_MULT, MPI_DOUBLE, MPI_SUM, 0, comm);
#endif
  
  sum=0.;
  for (int i=0; i<MAX_MULT; i++) sum += (double) pr_lh[1][i];
  
  for (int i=0; i<MAX_MULT; i++) {
    temp1[i] = pr_lh[1][i]/sum;
    temp2[i] = temp1[i]*temp1[i];
  }
  
  tmp1H = new double [MAX_MULT];
  tmp2H = new double [MAX_MULT];	

#ifdef MPIRUN
  MPI_Reduce (temp1, tmp1H, MAX_MULT, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce (temp2, tmp2H, MAX_MULT, MPI_DOUBLE, MPI_SUM, 0, comm);
#endif
  
  if (ip==0) { // processor rank 0

    for (int i=0; i<MAX_MULT; i++) {
      out << i << " " << tmp1[i] / np << "  " << errorMC (tmp2[i], tmp1[i], np) 
	  << " " << tmp1L[i] / np << "  " << errorMC (tmp2L[i], tmp1L[i], np) 
	  << " " << tmp1H[i] / np << "  " << errorMC (tmp2H[i], tmp1H[i], np) 
	  << endl;		
    }
    out << endl ;
  }
	
  delete [] tmp1;
  delete [] tmp2;
  delete [] temp1;
  delete [] temp2;
  delete [] tmp1L ; delete [] tmp1H ; delete [] tmp2L ; delete [] tmp2H ;

}

/*******************************************************************************
 * Print out correlation matrix rho(E1,E2) in output file.
 ******************************************************************************/
void resultsMC::printCorrelationResults (string filename)
{
  int gnuplotIndex=-1;

  int itmp[NUMBER_OF_PRIVATE_GRID];
  
  ofstream out;

  if(ip==0){
    out.open (&filename[0]);
    out << "#\n" << "# CGMF Results\n" << "#\n# [gnuplot# " << ++gnuplotIndex << "] Two particle correlations: energy\n" ;
  }
    out << setprecision(4) << setiosflags(ios::scientific);
  
  for(int i=0; i<3 ; i++ ){
    for(int j=0 ; j < NUMBER_OF_PRIVATE_GRID ; j++){
#ifdef MPIRUN
      MPI_Reduce (corr_e1_e2[i][j], itmp, NUMBER_OF_PRIVATE_GRID, MPI_INT, MPI_SUM, 0, comm);
#endif
      if( ip==0){
	for(int k=0; k< NUMBER_OF_PRIVATE_GRID; k++){
	  corr_e1_e2[i][j][k]=itmp[k];
	}
      }
    }
  }

  if(ip==0){
    for(int i=0; i< NUMBER_OF_PRIVATE_GRID; i++){
      for(int j=0; j< NUMBER_OF_PRIVATE_GRID; j++){
	out << .5 * ( custom_energy_grid[i] + custom_energy_grid[i+1] ) << " " << .5 * ( custom_energy_grid[j] + custom_energy_grid[j+1] ) << " " ;
	for(int k=0; k<3 ; k++ ) out << corr_e1_e2[k][i][j] << " "  ;
	out << endl ;
      }
      out << endl ;
    }
    out << endl ;
    out.close();
  }

}
