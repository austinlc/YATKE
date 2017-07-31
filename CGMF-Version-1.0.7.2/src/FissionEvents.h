/*
 *  CGMF
 *  [ FissionEvents.h ]
 *
 *  Version: 1.0.3
 *  April 26, 2013
 *
 *  P.Talou, talou@lanl.gov
 *
 *  Reads and analyzes output from CGMF.
 *
 */

#ifndef CGMF_FissionEvents_h
#define CGMF_FissionEvents_h
#endif

#ifdef MPIRUN
#include <mpi.h>
#endif

#ifndef __FISSIONFRAGMENTS_H__
#include "FissionFragments.h"
#endif

#ifndef __FFD_H__
#include "ffd.h"
#endif

#include <iostream>
#include <fstream>

using namespace std;

struct emittedParticleType
{
  
  char particleType; // 'n' (neutron) or 'g' (gamma) or 'i' (internal conversion)
  
  // to keep track of total number of emitted particles, also for light and heavy fragments
  int counts, lfCounts, hfCounts;
  
  // multiplicity and multiplicity distributions
  double nubar, lfNubar, hfNubar;
  double *Pnu, *lfPnu, *hfPnu;   // total, light fragment and heavy fragment P(nu)
  
  // average energies in center-of-mass and lab frames
  double Ecm, lfEcm, hfEcm;
  double Elab, lfElab, hfElab;
  
  // Angle between neutron and light fragment
  double *nLF_angles;
  double *nn_angles;
  
  // results as a function of fragment mass
  double *EcmA, *ElabA, *nuA;
  double *EcmAcut, *nuAcut;
  
  // results as a function of TKE
  double *EcmTKE, *ElabTKE, *nuTKE;
  
  // spectra in the center-of-mass and lab frames
  double *cmSpectrum;
  double *labSpectrum;
  
  double *ccLabSpectrum; // continuum-to-continuum
  double *cdLabSpectrum; // continuum-to-discrete
  double *ddLabSpectrum; // discrete-to-discrete
  
  double **ccLabSpectrumA; // continuum-to-continuum, per fission fragment mass
  double **cdLabSpectrumA; // continuum-to-discrete, per fission fragment mass
  double **ddLabSpectrumA; // discrete-to-discrete, per fission fragment mass
  
  // average spectra for light and heavy fragments
  double *cmLightFragmentSpectrum, *cmHeavyFragmentSpectrum;
  double *labLightFragmentSpectrum, *labHeavyFragmentSpectrum;
  
  // c.m. spectrum as a function of mass of primary fragment
  double **cmSpectrumA;
  
  // c.m. histograms as a function of mass of primary fragment
  double **cmHistogramA;
  
  // average c.m. energy as a function of mass and TKE
  double **EcmATKE;
  int **countsATKE;
  double **nubarATKE;
  
  // exclusive spectra for given multiplicities
  double **cmExclusiveSpectra; // function of multiplicity and Eout
  
  // exclusive c.m. spectra as a function of (A,nu)
  double ***cmSpectrumNuA;
  
  // exclusive c.m. histograms as a function of (A,nu)
  double ***cmHistogramNuA;
  
  // two-particle energy correlations
  double **energyCorrelations;
  
  // two-particle angular correlations
  double **angularCorrelations;
  
};

struct emissionType {
  int multiplicity;
  double cmEnergies   [MAX_NUMBER_PARTICLES];
  double labEnergies  [MAX_NUMBER_PARTICLES];
  double cmAngles     [MAX_NUMBER_PARTICLES];
  double labAngles    [MAX_NUMBER_PARTICLES];
  int transitionTypes [MAX_NUMBER_PARTICLES];
  double Vlab         [MAX_NUMBER_PARTICLES][3]; // (x,y,z) components in lab. frame
  double Vcm          [MAX_NUMBER_PARTICLES][3]; // (x,y,z) components of velocity in c.m. of fragment
};


struct fragmentEventType {
  int    A, Z;
  double mass; // ground-state mass (MeV)
  double KE;  // kinetic energy (MeV)
  double Ui;  // initial excitation energy (MeV)
  int    Pi;  // initial parity
  double Ji;  // initial spin
  double preMomentum[3];  // (x,y,z) components of initial fragment momentum
	double postMomentum[3]; // (x,y,z) components of final fragment momentum
  emissionType emissions[3]; // neutrons [0], gammas [1] and internal conversion [2]
};


class FissionEvents
{
  
  //------------------------------------------------------------------------------
public:
  //------------------------------------------------------------------------------
  
  FissionEvents (int);
  ~FissionEvents (void);
  
  int lfCounter, hfCounter;
  
  void setZAIDcn (int ZAIDt, double Einc);
  
  void addFragments (fissionFragmentType lightFragment, fissionFragmentType heavyFragment);
  
  void addFragments (int Zl, int Al, double Uil, double Jil, int Pil, double KEl,
                     int Zh, int Ah, double Uih, double Jih, int Pih, double KEh);
  
  void saveFragmentEvent (int Af, int Zf, double initialEnergy, double initialSpin, int initialParity,
                          int neutronMultiplicity, double *neutronEnergies,
                          int gammaMultiplicity, double *gammaEnergies, int *gammaTypes,
                          int icMultiplicity, double *icEnergies);
  
  void computeFinalResults (void);
  void computeFinalResults (emittedParticleType *);
  
  void writeHistories (string outputFilename, int, double, int, double);
  void readHistories (string inputFilename, int *numberEvents, int *ZAIDt, double *Einc, double *alphaI);
  void analyzeResults (int ZAID);
  
  void printSummaryResults (string outputFilename, int, double, int, double);
  void printCompactSummaryResults (void);
  void saveResultsToGnuplot (void);
  
  void printNubarA (string outputFilename);
  
#ifdef MPIRUN
  void setMPI ();
#endif
  
#if defined(MPIRUN) && defined(GEANT)
  void saveResultsGEANT (void);
#endif
  
  fragmentEventType *lightFragments, *heavyFragments;
  double **preFissionNeutrons;
  
  // used only if HISTORIES defined
  ofstream historyFile;
  
  //------------------------------------------------------------------------------
private:
  //------------------------------------------------------------------------------
  
  void init                 (int);
  void deleteParticles      (emittedParticleType *);
  void initEmittedParticles (emittedParticleType *);
  void saveParticles        (emittedParticleType *particles, int lfnu, double *lfcmEnergies, int hfnu, double *hfcmEnergies);
  
  void saveParticleResultsToDataFile (emittedParticleType, string);
  void saveMonteCarloYieldsToFile    (string);
  void saveGeneralResults            (string);
  
  template<int N> int findEnergyIndex (double x0, double (&xarray)[N]);
  
  int ip, np; // index and number of processors
  
  int eventCounter;
  
  emittedParticleType neutrons, gammas, ic;
  fragmentEventType fragmentEvents[2]; // events for light [0] and heavy [1] fragments
  
  double incidentNeutronEnergy; // [MeV]
  
  int Ac, Zc;     // compound (fissioning) nucleus
  int Asym, Zsym; // symmetric fragment
  
  int Zl, Al, Zh, Ah; // light and heavy fragments
  double Ml, Mh;      // ground-state masses of the fragments
  double KEl, KEh;    // kinetic energies
  double Ul, Uh;      // excitation energies
  double Jl, Jh;      // spins
  int pl, ph;         // parities
  
  double momentumLF[3], momentumHF[3]; // (x,y,z) components of light and heavy fragment momenta
  
  double Ef; // fission fragment kinetic energy (follows the decay chain and recoils)
  
  double aveMass, aveCharge, aveTKE, aveTXE;
  double aveMassLF, aveMassHF;
  double aveChargeLF, aveChargeHF;
  double aveKEl, aveKEh;
  
  double YZ [NUMZ];
  double YA [NUMA];
  double YApost [NUMA]; // post-neutron emission mass yields
  
  double YUl   [NUMBER_EXCITATION_ENERGY_GRID];
  double YUh   [NUMBER_EXCITATION_ENERGY_GRID];
  double YUtot [NUMBER_EXCITATION_ENERGY_GRID];
  
  double YKEl [NUMBER_EXCITATION_ENERGY_GRID];
  double YKEh [NUMBER_EXCITATION_ENERGY_GRID];
  double YTKE [NUMBER_EXCITATION_ENERGY_GRID];
  
  double UiA [NUMA];
  double initialSpinA [NUMA];
  
  double initialSpinTKE [NUMBER_EXCITATION_ENERGY_GRID];
  
  double **YUiA; // Initial excitation energy distribution as a function of fragment mass
  
  // angular distributions of evaporated neutrons
  double YThetalab [NUMANGLES];
  double YThetacm  [NUMANGLES];
  
  // as a function of (int) 2*spin
  double lfPspin[100];
  double hfPspin[100];
  
  int neutronLFMultiplicity, gammaLFMultiplicity;
  int numberFissionEvents;
  
  double totalGammaEnergyA [NUMA];
  double totalGammaEnergyTKE [NUMBER_EXCITATION_ENERGY_GRID];
  
  double totalGammaEnergyNu [20];  // <Eg>(nu)
  double gammaMultiplicityNu [20]; // <Ng>(nu)
  int countsNu[20]; // number of counts per nu (to compute <Eg> and <Ng> vs. neutron multiplicity)
  
  double spectrumEnergyGrid[NUMBER_SPECTRUM_ENERGY_GRID];
  double spectrumEnergySteps[NUMBER_SPECTRUM_ENERGY_GRID];
  
  // for GEANT-4 simulation of DANCE (LANSCE)
  double **gammaGEANT;
  int indexGEANT;
  
  int numberPFN; // number of pre-fission neutrons
  double aveEnergyPFN; // average energy of pre-fission neutrons
  
};

