/*
 *  CGMF
 *  Version: 1.0.6
 *
 *  [ FissionFragments.h ]
 * 
 *  Generates primary fission fragment yields to be sampled in CGMF. It either
 *  reads data from input file, or produce them from systematics.
 *
 */

#ifndef __FISSIONFRAGMENTS_H__
#define __FISSIONFRAGMENTS_H__
#endif

#include <iostream>

// To include with XML I/O
//#include <libxml/parser.h>
//#include <libxml/tree.h>

#include "config-ff.h"
#include "Yields.h"
#include "BrosaYields.h"

using namespace std;

#ifndef __STRUCTUR_H__
#define __STRUCTUR_H__
#include "structur.h"
#endif

struct fissionFragmentType {
  int Z, A, N;
  double mass; // ground-state mass (MeV)
  double U;
  double spin;
  int parity;
  double KE;
  double TKE, energyRelease; // also depend on complementary fragment
  double p[4]; // 4-vector energy-momentum relativistic notation; p[0]=energy
  double pc[3]; // (classical) momentum components (x,y,z) in lab. frame
};


class FissionFragments
{
  
  /* PUBLIC */
  public :
	
	// FOR TESTING TEMPERATURES PURPOSE ONLY ---------------------------------------------------- /////////////////////
	void computeTemperatures (int, int);
	void readTemperatures (void);
	void readLevelDensityParameters (void);
	
	double temp [2000][32];
	double ldp [2000][32];
	int ZAIDlist [2000];
	
  FissionFragments (int ZAID, double incidentEnergy, double *alphaSpin);
  FissionFragments (int ZAID, double Einc, double *alphaSpin, double RTin);
  FissionFragments (int ZAID, double incidentEnergy, double *alphaSpin, string yieldsFile);
  FissionFragments (int ZAIDf, double excitationEnergy, const string);
  FissionFragments (void);
	
  ~FissionFragments (void);
  
  void writeYieldsInFile (string outputFilename);
  void generateInitialFissionFragmentHistories (string outputFilename, int numberEvents);
  void generateInitialFissionFragmentHistories (fissionFragmentType*, fissionFragmentType*, const int numberEvents);
  void generateSingleFissionFragments (fissionFragmentType*, fissionFragmentType*, int ZAIDf, double Eexc, int numberEvents);
  void checkDistributions (string inputFilename, string outputFilename);
  
  void generateYieldsForSensitivityStudies (string YATKEfilename); // for work with Ramona (LLNL) and Jorgen (LBNL)
  
  void studyEnergySorting (void);

  void sortExcitationEnergy (double Utot, int Al, int Ah, double *Ul, double *Uh, double RT); // for Carjan

  int Amin, Amax;       // min and max of fragment mass
  int Zmin, Zmax;       // min and max of fragment charge
  int TKEmin, TKEmax;   // min and max of total kinetic energy (MeV)
  int Zpmin, Zpmax;     // min and max of most probable charges
  
  int ZAIDc;            // ZAID of compound (fissioning) nucleus
  int ZAIDt;            // ZAID of target nucleus
  int Ac, Zc;           // mass and charge of compound (fissioning) nucleus
  int At, Zt;           // mass and charge of target nucleus
  int Ap, Zp;           // mass and charge of projectile
  int Asym, Zsym;       // mass and charge of symmetric fragment
	
	int ZAIDc0, Ac0, Zc0; // Original (before any pre-fission neutron emission) compound fissioning nucleus

  double incidentEnergy;  // incident neutron energy [MeV] or equivalent incident neutron energy if a neutron is emitted
	double incidentEnergy0; // incident neutron energy [MeV] BEFORE any pre-fission neutron emission
	
  double YA[NUMA];
  double YATKE[NUMA][NUMTKE];
  double YTKE[NUMTKE];
  double YZ[NUMZ];
  double YZA[NUMdZ][NUMA];
  double YZA2[NUMZ][NUMA];

  // Average Q-values as a function of TKE and A
  double QfTKE[NUMTKE];
  double QfA[NUMA];
  
  double maxYieldZ[NUMA];
  double maxYieldATKE;
	
	double maxYieldA;	// << 1.0.6 >>
	
  double beta2 [NUMZ][NUMA];

  double YZATKE[NUMZ][NUMA][NUMTKE];
  
  /* PRIVATE */
  private :
  
  void init (void);
  void setOptions (void);
	
  void readMainInputFile (string inputFilename);
//  void createDefaultInputFile (string inputFilename, int ZAID, double incidentEnergy, double alphaSpin);
  
  // To include with XML I/O
  //  void readXMLInputFile (string inputFilename);
  //	void createDefaultXMLInputFile (string inputFilename, int ZAID, double incidentEnergy);
  
  void buildYields (void);
  void readYieldsATKE (string);
	void readYieldsAZTKE (string);
  
  void sampleFissionFragments (fissionFragmentType *lf, fissionFragmentType *hf);

	// << 1.0.6 >>
	void sampleFissionFragments (fissionFragmentType *lf, fissionFragmentType *hf, int);
	void sampleFissionFragments (fissionFragmentType *lightFragment, fissionFragmentType *heavyFragment, Yields *ffy) ;

  // brosa sampling
  void sampleFissionFragmentsBROSA (fissionFragmentType *lightFragment, fissionFragmentType *heavyFragment, BrosaYields *BY) ;

	
  // << 1.0.7 >> added fragment momentum calculation
  void computeFragmentInitialConditions (double TKE, int Zl, int Al, int Zh, int Ah, double *Ul, double *Uh,
                                         double *spinl, double *spinh, int *parl, int *parh, double *energyRelease,
                                         double *momentum);
    
  void computeFragmentExcitationEnergy (double TXE, int Zl, int Al, int Zh, 
                                        int Ah, double *Ul, double *Uh, string sortingOption);
  
  void getInitialSpin (int Z, int A, double *U, double *spin);
  
  void computeWahlParametersForChargeDistribution (void);
  void buildWahlYA (void);
  void buildYATKEfromSystematics (void);
  void buildSystematicsTKEvsMass (void);
	
	// << 1.0.6 >>
  void buildYATKEfromSystematics (Yields *);
  void buildYATKE (Yields *);
	
  void buildYA (void);
  void buildYZA (void);
  void buildYZ (void);
  void buildYTKE (void);
  
  void resetYields (void);
  
  // energy sorting methods & variables
  double Pfactor (int Z, int N);
  void computePfactors (void);
  void computeLevelDensityTables (void);
  void computeLevelDensityParameterTables (void);
  
  void readLevelDensityParameterTables (string filename);
	void readLevelDensityTables (string filename);
    
  double computeEnergyFromMaxEntropy (int Zl, int Al, double TIXE);
  
  void readRTAh (int Zt, int At);
  void readRTAh (int Zt, int At, string filename);
  void readMollerDeformationEnergies (string filename);  
  
  double Pfactors[NUMA][NUMZ];
  double RTAh[NUMA];
  
  // level density tables f(A,Z,E)
  double ***levelDensities;
  double ***levelDensityParameters;
  double ldEnergyGrid [NUME];
	double ldEnergyGrid2 [NUME2];
  
  //-- IONEL
  //	double beta2 [NUMZ][NUMA]; // put in PUBLIC 
  void readDeformations (void);
  
  double ***temperatures;
  double spinParameter (double, int, int);
	void setSpinStrength (double); // Ionel, new version of high-Einc yields
  void computeTemperatureTables (void);

  
  Nucleus *fragments; // light and heavy fragments following CGM structure
  
  void readNeutronBindingEnergies (void);
  void readMasses (void);
  void computeNeutronSeparationEnergies (void);
  
  // Wahl's parameters
  double Z0[NUMA];
  double sZ0[NUMA];
  double FZZ0[NUMA];
  double FNZ0[NUMA];
  
  double sigmaZ;
  int dZ;
  
  double meanTKE[NUMA];
  double sigTKE[NUMA];
  
  string energySortingOption;
  double RT;
	double ratioU;
  double alphaI;

  double EdefA [NUMA]; // deformation energies (MeV)
  
  int lightFragmentSpinParameter;
  int heavyFragmentSpinParameter;
  
  string yieldsOption;
  string yieldsFilename;

  double B1n[NUMZ][NUMA];
  double B2n[NUMZ][NUMA];
  
  double S1n[NUMZ][NUMA];
  double S2n[NUMZ][NUMA];
  
  double ExpMassMeV[NUMZ][NUMA]; // excess masses in MeV
  double ExpMassAMU[NUMZ][NUMA]; // excess masses in atomic mass units
  
  double massExcessCompound; // excess mass in MeV for the compound fissioning nucleus
  double SnCompound; // neutron separation energy for the compound fissioning nucleus

	// << 1.0.6 >>
	
	double SnCompound0; // neutron separation energy for the original (before any pre-fission neutron emission) compound fissioning nucleus
	
	void constructYields (Yields *ffy);        // constructs yields for multi-chance fission
	void cgmBinaryNeutron (int nemit,int c0);  // this is a wrapper similar to specMain
	
	// pre-fission neutron emission << 1.0.6 >>
	
	double *sp1;   // total spectrum for the first neutron emitted
	double *sp;    // evaporation spectrum from the neutrons
	Pdata pdt[MAX_CHANNEL];
	Transmission *tc;
	int nemit;
	double * barrier, * emissprob, emissprob_max;

	void readMultichanceFissionData (void); // void read_multifission_data(void);
	int getPrefissionNeutrons (double *); // int get_prefission_neutrons(double *);

	double *energy_grid1;
	double *sp_pfn;
	double *energy_pfn ;// the energy of the generated prefission neutrons
	double prCheck[10];
	int **sp_pfn_mult;
	
	int index_save;
	
  template<int N> int findEnergyIndex (double x0, double (&xarray)[N]);
  
}; // end class Fission Fragments
