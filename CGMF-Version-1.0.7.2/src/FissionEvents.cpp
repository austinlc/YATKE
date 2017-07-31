//
//  fissionEvents.cpp
//  CGMF
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>

#include "mt19937ar.h"
#include "FissionFragments.h"
#include "FissionEvents.h"
#include "physics.h"
#include "masstable.h"

#ifndef __FFD_H__
#include "ffd.h"
#endif

#include "global-mcnp.h"
#include "global.h"

const int NUMBER_OF_PRIVATE_GRID = 565;
#define GRID_STRUCTURE_FILE "privategrid1.h"
#include GRID_STRUCTURE_FILE

FissionEvents::FissionEvents (int maxNumberEvents) { init(maxNumberEvents); }

FissionEvents::~FissionEvents (void) {
  
  //deleteParticles (&neutrons);
  //deleteParticles (&gammas);
  
}


/*******************************************************************************
 Initialization of an instance of FissionEvents.
 
 - Sets all arrays to zero
 - Defines energy grids on which results are to be saved
 - Initializes arrays to save neutron and gamma data
 
 ******************************************************************************/

void FissionEvents::init (int maxNumberEvents) {
  
#ifdef MPIRUN
  setMPI ();
#else
  ip=0; np=1;
#endif
  
  numberFissionEvents = 0;
  
  eventCounter = -1;
  lfCounter = -1;
  hfCounter = -1;
  
  try {
    lightFragments     = new fragmentEventType [maxNumberEvents];
    heavyFragments     = new fragmentEventType [maxNumberEvents];
    preFissionNeutrons = new double * [maxNumberEvents]; // << 1.0.6 >>
    for (int i=0; i<maxNumberEvents; i++) {
      preFissionNeutrons[i] = new double [5]; // max. of 5 pre-fission neutrons at this point
    }
  }
  catch(bad_alloc){
    cerr << "ERROR     :memory allocation error";
    exit(-1);
  }
  
  // define spectrum energy grid [original "Madland-Nix model" grid; 90 steps per decade]
  // modified to get finer grid at highest outgoing energies
  // nen = 641
  // numberStepsPerDecade = 90
  // Emax=30.0_rk
  // -- less number of points per decade
  // nen=231
  // numberStepsPerDecade=30
  
  /*  int numberStepsPerDecade = 90;
   
   int k;
   k=-1;
   double Emin = 1.0e-5;
   for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) {
   double E = pow(10.0,i)*Emin;
   for (int j=1; j<=numberStepsPerDecade; j++) {
   k++;
   spectrumEnergySteps[k] = (pow(10.0,i+1)-pow(10.0,i))*Emin/float(numberStepsPerDecade);
   spectrumEnergyGrid[k]  = E + float(j-1) * spectrumEnergySteps[k] ;
   if (spectrumEnergyGrid[k]>10.0) { // finer grid (dE=200 keV)
   spectrumEnergySteps[k]=0.2;
   spectrumEnergyGrid[k] = spectrumEnergyGrid[k-1]+0.2;
   }
   if (k>NUMBER_SPECTRUM_ENERGY_GRID-2) break;
   }
   if (k>NUMBER_SPECTRUM_ENERGY_GRID-2) break;
   }
   */
  
  spectrumEnergyGrid[0]=1.0e-5;
  for (int i=1; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) {
    spectrumEnergyGrid[i] = spectrumEnergyGrid[i-1]*1.025;
    spectrumEnergySteps[i] = spectrumEnergyGrid[i]-spectrumEnergyGrid[i-1];
    //cout << i << " " << spectrumEnergyGrid[i] << " " << spectrumEnergySteps[i] << "\n";
  }
  spectrumEnergySteps[0] = spectrumEnergySteps[1];
  //exit(0);
  
  for (int i=0; i<NUMZ; i++) { YZ[i]=0.0; }
  for (int i=0; i<NUMA; i++) {
    YA[i] = 0.0;
    UiA[i] = 0.0;
    initialSpinA[i] = 0.0;
  }
  
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    YTKE[i]=YKEl[i]=YKEh[i]=0.0;
    YUtot[i]=YUl[i]=YUh[i]=0.0;
    initialSpinTKE[i]=0.0;
  }
  
  for (int i=0; i<NUMANGLES; i++) {
    YThetalab[i]=0.0;
    YThetacm[i]=0.0;
  }
  
  YUiA = new double * [NUMA];
  for (int i=0; i<NUMA; i++) {
    YUiA[i] = new double [NUMBER_EXCITATION_ENERGY_GRID];
    std::fill_n (YUiA[i], NUMBER_EXCITATION_ENERGY_GRID, 0.0);
  }
  
  // Proba(2 * initial spin)
  for (int i=0; i<100; i++) {
    hfPspin[i] = 0.0;
    lfPspin[i] = 0.0;
  }
  
  // Initialize emitted particle structures (neutrons and gammas)
  initEmittedParticles (&gammas);   gammas.particleType='g';
  initEmittedParticles (&neutrons); neutrons.particleType='n';
  initEmittedParticles (&ic); ic.particleType='i';
  
  std::fill_n (totalGammaEnergyTKE, NUMBER_EXCITATION_ENERGY_GRID, 0.0);
  std::fill_n (totalGammaEnergyA, NUMA, 0.0);
  std::fill_n (totalGammaEnergyNu, 20, 0.0);
  std::fill_n (gammaMultiplicityNu, 20, 0.0);
  std::fill_n (countsNu, 20, 0.0);
  
  
  
#if defined(GEANT) && defined(MPIRUN)
  // Initialize array of results for GEANT simulations (J.Ullmann)
  indexGEANT=-1;
  int maxNumberEvents = 100000;
  gammaGEANT = new double * [maxNumberEvents];
  for (int i=0; i<maxNumberEvents; i++) {
    gammaGEANT[i] = new double [20]; // 20 = max. g-ray multiplicity for 1 fragment
    std::fill_n (gammaGEANT[i], 20, 0.0);
  }
#endif
  
}

/*******************************************************************************
 Initialize the neutron and gamma-ray structures (emittedParticleType).
 ******************************************************************************/

void FissionEvents::initEmittedParticles (emittedParticleType *particles) {
  
  particles->counts   = 0;
  particles->lfCounts = 0;
  particles->hfCounts = 0;
  
  particles->Ecm   = 0.0;
  particles->lfEcm = 0.0;
  particles->hfEcm = 0.0;
  
  particles->Elab   = 0.0;
  particles->lfElab = 0.0;
  particles->hfElab = 0.0;
  
  particles->EcmA    = new double [NUMA];
  particles->ElabA   = new double [NUMA];
  particles->nuA     = new double [NUMA];
  
  particles->EcmAcut = new double [NUMA];
  particles->nuAcut  = new double [NUMA];
  
  for (int i=0; i<NUMA; i++) {
    particles->EcmA[i] = 0.0;
    particles->ElabA[i] = 0.0;
    particles->nuA[i] = 0.0;
    particles->EcmAcut[i] = 0.0;
    particles->nuAcut[i] = 0.0;
  }
  
  particles->EcmTKE    = new double [NUMBER_EXCITATION_ENERGY_GRID];
  particles->ElabTKE   = new double [NUMBER_EXCITATION_ENERGY_GRID];
  particles->nuTKE     = new double [NUMBER_EXCITATION_ENERGY_GRID];
  
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    particles->EcmTKE[i] = 0.0;
    particles->ElabTKE[i] = 0.0;
    particles->nuTKE[i] = 0.0;
  }
  
  particles->Pnu   = new double [NUMMULT];
  particles->lfPnu = new double [NUMMULT];
  particles->hfPnu = new double [NUMMULT];
  
  for (int i=0; i<NUMMULT; i++) {
    particles->Pnu[i]   = 0.0;
    particles->lfPnu[i] = 0.0;
    particles->hfPnu[i] = 0.0;
  }
  
  particles->cmSpectrum  = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  particles->labSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  
  particles->ccLabSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  particles->cdLabSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  particles->ddLabSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  
  particles->cmLightFragmentSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  particles->cmHeavyFragmentSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  
  particles->labLightFragmentSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  particles->labHeavyFragmentSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  
  for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) {
    particles->cmSpectrum[i] = 0.0;
    particles->labSpectrum[i] = 0.0;
    particles->ccLabSpectrum[i] = 0.0;
    particles->cdLabSpectrum[i] = 0.0;
    particles->ddLabSpectrum[i] = 0.0;
    particles->cmLightFragmentSpectrum[i] = 0.0;
    particles->cmHeavyFragmentSpectrum[i] = 0.0;
    particles->labLightFragmentSpectrum[i] = 0.0;
    particles->labHeavyFragmentSpectrum[i] = 0.0;
  }
  
  
  
  particles->ddLabSpectrumA = new double * [NUMA];
  for (int i=0; i<NUMA; i++) { particles->ddLabSpectrumA[i] = new double [NUMBER_SPECTRUM_ENERGY_GRID]; }
  particles->cdLabSpectrumA = new double * [NUMA];
  for (int i=0; i<NUMA; i++) { particles->cdLabSpectrumA[i] = new double [NUMBER_SPECTRUM_ENERGY_GRID]; }
  particles->ccLabSpectrumA = new double * [NUMA];
  for (int i=0; i<NUMA; i++) { particles->ccLabSpectrumA[i] = new double [NUMBER_SPECTRUM_ENERGY_GRID]; }
  
  particles->cmSpectrumA = new double * [NUMA];
  for (int i=0; i<NUMA; i++) {
    particles->cmSpectrumA[i] = new double [NUMBER_SPECTRUM_ENERGY_GRID];
    for (int j=0; j<NUMBER_SPECTRUM_ENERGY_GRID; j++) {
      particles->cmSpectrumA[i][j] = 0.0;
    }
  }
  
  // energy grid for histograms from 0.0 to 15.0 MeV by steps of 0.05 MeV (301 points); half as Alf :) Hambsch, 2014
  particles->cmHistogramA = new double * [NUMA];
  for (int i=0; i<NUMA; i++) {
    particles->cmHistogramA[i] = new double [301];
    for (int j=0; j<301; j++) {
      particles->cmHistogramA[i][j] = 0.0;
    }
  }
  
  // histogram(nu,A,Ecm)
  particles->cmHistogramNuA = new double ** [5];
  for (int i=0; i<5; i++) {
    particles->cmHistogramNuA[i] = new double * [NUMA];
    for (int j=0; j<NUMA; j++) {
      particles->cmHistogramNuA[i][j] = new double [301];
      for (int k=0; k<301; k++) {
        particles->cmHistogramNuA[i][j][k] = 0.0;
      }
    }
  }
  
  
  
  particles->cmExclusiveSpectra = new double * [10];
  for (int i=0; i<10; i++) {
    particles->cmExclusiveSpectra[i] = new double [NUMBER_SPECTRUM_ENERGY_GRID];
    for (int j=0; j<NUMBER_SPECTRUM_ENERGY_GRID; j++) {
      particles->cmExclusiveSpectra[i][j] = 0.0;
    }
  }
  
  // spectrum(nu,A,Ecm)
  particles->cmSpectrumNuA = new double ** [5];
  for (int i=0; i<5; i++) {
    particles->cmSpectrumNuA[i] = new double * [NUMA];
    for (int j=0; j<NUMA; j++) {
      particles->cmSpectrumNuA[i][j] = new double [NUMBER_SPECTRUM_ENERGY_GRID];
      for (int k=0; k<NUMBER_SPECTRUM_ENERGY_GRID; k++) {
        particles->cmSpectrumNuA[i][j][k] = 0.0;
      }
    }
  }
  
  
  particles->energyCorrelations = new double * [NUMBER_OF_PRIVATE_GRID];
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) {
    particles->energyCorrelations[i] = new double [NUMBER_OF_PRIVATE_GRID];
    for (int j=0; j<NUMBER_OF_PRIVATE_GRID; j++) {
      particles->energyCorrelations[i][j] = 0.0;
    }
  }
  
  particles->angularCorrelations = new double * [NUMANGLES];
  for (int i=0; i<NUMANGLES; i++) {
    particles->angularCorrelations[i] = new double [NUMANGLES];
    for (int j=0; j<NUMANGLES; j++) {
      particles->angularCorrelations[i][j] = 0.0;
    }
  }
  
  particles->nLF_angles = new double [NUMANGLES];
  particles->nn_angles  = new double [NUMANGLES];
  
  for (int i=0; i<NUMANGLES; i++) {
    particles->nLF_angles[i] = 0.0;
    particles->nn_angles[i] = 0.0;
  }
  
  // To compare with Vorobyev's data <Ecm>(A,TKE)
  particles->EcmATKE    = new double * [NUMA];
  particles->nubarATKE  = new double * [NUMA];
  particles->countsATKE = new int    * [NUMA];
  for (int i=0; i<NUMA; i++) {
    particles->EcmATKE[i]    = new double [NUMBER_EXCITATION_ENERGY_GRID];
    particles->nubarATKE[i]  = new double [NUMBER_EXCITATION_ENERGY_GRID];
    particles->countsATKE[i] = new int    [NUMBER_EXCITATION_ENERGY_GRID];
    
    for (int j=0; j<NUMBER_EXCITATION_ENERGY_GRID; j++) {
      particles->EcmATKE[i][j]    = 0.0;
      particles->nubarATKE[i][j]  = 0.0;
      particles->countsATKE[i][j] = 0;
    }
    
  }
  
}



/*!
 
 Add light and heavy fragment characteristics to a FissionEvents
 
 */
void FissionEvents::addFragments (fissionFragmentType lf, fissionFragmentType hf) {
  
  numberFissionEvents++; // increment number of fission events
  
  int k, kl, kh;
  
  Al  = lf.A;
  Zl  = lf.Z;
  KEl = lf.KE;
  Ul  = lf.U;
  Ml  = lf.mass;
  
  Ah  = hf.A;
  Zh  = hf.Z;
  KEh = hf.KE;
  Uh  = hf.U;
  Mh  = hf.mass;
  
  YZ[Zl]++;
  YZ[Zh]++;
  YA[Al]++;
  YA[Ah]++;
  
  kl = findEnergyIndex (KEl, excitationEnergyGrid);    YKEl[kl]++;
  kh = findEnergyIndex (KEh, excitationEnergyGrid);    YKEh[kh]++;
  k  = findEnergyIndex (KEl+KEh, excitationEnergyGrid); YTKE[k]++;
  
  kl = findEnergyIndex (Ul, excitationEnergyGrid);      YUl[kl]++;
  kh = findEnergyIndex (Uh, excitationEnergyGrid);      YUh[kh]++;
  k  = findEnergyIndex (Ul+Uh, excitationEnergyGrid);  YUtot[k]++;
  
  YUiA [Al][kl]++;
  YUiA [Ah][kh]++;
  
  UiA[Al] += Ul;
  UiA[Ah] += Uh;
  
  initialSpinA[Al] += lf.spin;
  initialSpinA[Ah] += hf.spin;
  
  lfPspin[(int) (2*lf.spin)]++;
  hfPspin[(int) (2*hf.spin)]++;
  
  for (int i=0; i<3; i++) {
    momentumLF[i] = lf.pc[i];
    momentumHF[i] = hf.pc[i];
  }

}

void FissionEvents::addFragments (int Zl, int Al, double Ul, double Jl, int Pil, double KEl,
                                  int Zh, int Ah, double Uh, double Jh, int Pih, double KEh) {
  
  int kl, kh, k;
  YZ[Zl]++;
  YZ[Zh]++;
  YA[Al]++;
  YA[Ah]++;
  
  kl = findEnergyIndex (KEl, excitationEnergyGrid);    YKEl[kl]++;
  kh = findEnergyIndex (KEh, excitationEnergyGrid);    YKEh[kh]++;
  k  = findEnergyIndex (KEl+KEh, excitationEnergyGrid); YTKE[k]++;
  
  initialSpinTKE[k] += (Jl+Jh);
  
  kl = findEnergyIndex (Ul, excitationEnergyGrid);      YUl[kl]++;
  kh = findEnergyIndex (Uh, excitationEnergyGrid);      YUh[kh]++;
  k  = findEnergyIndex (Ul+Uh, excitationEnergyGrid);  YUtot[k]++;
  
  YUiA [Al][kl]++;
  YUiA [Ah][kh]++;
  
  UiA[Al] += Ul;
  UiA[Ah] += Uh;
  
  initialSpinA[Al] += Jl;
  initialSpinA[Ah] += Jh;
  
  lfPspin[(int) (2*Jl)]++;
  hfPspin[(int) (2*Jh)]++;
  
}




/*******************************************************************************
 Save the data related to a particular fission fragment, light (index=0) or
 heavy (index=1). Those data are the neutron multiplicity and energies, and
 gamma multiplicity and energies.
 ******************************************************************************/
void FissionEvents::saveFragmentEvent (int Af, int Zf, double initialEnergy, double initialSpin, int initialParity,
                                       int neutronMultiplicity, double *neutronEnergies,
                                       int gammaMultiplicity, double *gammaEnergies, int *gammaTypes,
                                       int icMultiplicity, double *icEnergies)
{
  
  int index;
  if (Af<=Asym && lfCounter==hfCounter) { index=0; lfCounter++; } else { index=1; hfCounter++; }
  
  fragmentEvents[index].A  = Af;
  fragmentEvents[index].Z  = Zf;
  fragmentEvents[index].Ui = initialEnergy;
  fragmentEvents[index].Ji = initialSpin;
  fragmentEvents[index].Pi = initialParity;
  
  fragmentEvents[index].emissions[0].multiplicity = neutronMultiplicity;
  for (int i=0; i<neutronMultiplicity; i++) { fragmentEvents[index].emissions[0].cmEnergies[i] = neutronEnergies[i]; }
  
  fragmentEvents[index].emissions[1].multiplicity = gammaMultiplicity;
  for (int i=0; i<gammaMultiplicity; i++) {
    fragmentEvents[index].emissions[1].cmEnergies[i]      = gammaEnergies[i];
    fragmentEvents[index].emissions[1].transitionTypes[i] = gammaTypes[i];
  }
  
  fragmentEvents[index].emissions[2].multiplicity = icMultiplicity;
  for (int i=0; i<icMultiplicity; i++) { fragmentEvents[index].emissions[2].cmEnergies[i] = icEnergies[i]; }
  
  double phi;
  double cosCmTheta, cmTheta;
  double Ecm, Elab, labTheta;
  double KEf;
  
  fragmentEventType *event;
  
  if (index==0) {
    event = &lightFragments[lfCounter];
    event->KE=KEl;
    event->mass=Ml;
    for (int i=0; i<3; i++) event->preMomentum[i] = momentumLF[i]; // classical
  } else {
    event = &heavyFragments[hfCounter];
    event->KE=KEh;
    event->mass=Mh;
    for (int i=0; i<3; i++) event->preMomentum[i] = momentumHF[i];
  }
	
  event->A  = Af;
  event->Z  = Zf;
  event->Ui = initialEnergy;
  event->Ji = initialSpin;
  event->Pi = initialParity;
  
  event->emissions[0].multiplicity = neutronMultiplicity;
  KEf = event->KE;
  
  // classical momentum components (x,y,z)
  double pfx = event->preMomentum[0];
  double pfy = event->preMomentum[1];
  double pfz = event->preMomentum[2];
  
  double pf = sqrt (pfx*pfx+pfy*pfy+pfz*pfz);
  
  double Mf = Af*amuMeV+mass_excess(Zf,Af);
  double Mn = neutronMass*amuMeV;
  
  // Fragment direction
  double costhetaF = pfz/pf;
  double sinthetaF = sqrt(1.0-costhetaF*costhetaF); // abs(pfx/pf);
  double phiF;
	
//  cout << "Pf: " << pfx << " " << pfy << " " << pfz << " -- " << pf << " " << pf*pf/(2.0*Mf) << " " << KEf << endl;
  
  if (pfx<0) {
    phiF = PI-atan(-pfy/pfx);
  } else if (pfy<0) {
    phiF = PI2+atan(pfy/pfx);
  } else {
    phiF = atan(pfy/pfx);
  }
  
  // in reference frame of fragment
  double pfx2 = sinthetaF*cos(phiF)*pfx+sinthetaF*sin(phiF)*pfy+costhetaF*pfz; // along e_r
  double pfy2 = costhetaF*cos(phiF)*pfx + sin(phiF)*costhetaF*pfy - sinthetaF*pfz; // along e_theta
  double pfz2 = -sin(phiF)*pfx+cos(phiF)*pfy; // along e_phi
  
  // recover pf in original (x,y,z) lab frame
  //	cout << sinthetaF*cos(phiF)*pfx2+costhetaF*cos(phiF)*pfy2-sin(phiF)*pfz2 << endl;
  //	cout << sinthetaF*sin(phiF)*pfx2+costhetaF*sin(phiF)*pfy2+cos(phiF)*pfz2 << endl;
  //	cout << costhetaF*pfx2-sinthetaF*pfy2 << endl;
  
  //  cout << neutronMultiplicity << " " << index << endl;
  
  for (int i=0; i<neutronMultiplicity; i++) {
    
    Mf = Af*amuMeV+mass_excess(Zf,Af);
    
    Ecm = neutronEnergies[i];
    
    // direction of neutron emission chosen isotropic in c.m. of fission fragment
    phi = twopi*genrand_real1();
    cosCmTheta = 2.0*genrand_real1()-1.0;
    cmTheta = acos(cosCmTheta);
    
    // momentum neutron in c.m.
    double pn = sqrt (2.0*Ecm*Mn);
    
    // neutron in c.m. of the fragment
    double pnx2 = pn*cosCmTheta;
    double pny2 = pn*sin(cmTheta)*cos(phi);
    double pnz2 = pn*sin(cmTheta)*sin(phi);
    
    // record neutron event in c.m. of fragment
    event->emissions[0].Vcm[i][0] = pnx2/Mn;
    event->emissions[0].Vcm[i][1] = pny2/Mn;
    event->emissions[0].Vcm[i][2] = pnz2/Mn;
    
    // vector addition in c.m. of fragment --> lab.
    double vlx2 = pnx2/Mn + pfx2/Mf;
    double vly2 = pny2/Mn + pfy2/Mf;
    double vlz2 = pnz2/Mn + pfz2/Mf;
    
//    double vl = sqrt(vlx2*vlx2+vly2*vly2+vlz2*vlz2);
    
//    cout << "v: " << vlx2/vl << " " << vly2/vl << " " << vlz2/vl << " -- " << index << " " << pfx2 << endl;
//		cout << "energies: " << Ecm << " " << 0.5*Mn*vl*vl << " " << cmTheta << " " << acos(vlx/vl) << " >>> ";

    // revert neutron vector to lab. frame
    double vlx = sinthetaF*cos(phiF)*vlx2+costhetaF*cos(phiF)*vly2-sin(phiF)*vlz2;
    double vly = sinthetaF*sin(phiF)*vlx2+costhetaF*sin(phiF)*vly2+cos(phiF)*vlz2;
    double vlz = costhetaF*vlx2-sinthetaF*vly2;

//    cout << vl << " " << sqrt(vlx*vlx+vly*vly+vlz*vlz) << endl;
//    cout << "vlab: " << vlx/vl << " " << vly/vl << " " << vlz/vl << endl;
    
    
    //    Elab = Ecm + KEf/Af + 2.0*sqrt(KEf/Af*Ecm)*cosCmTheta; // Eq.(1) Terrell, Phys. Rev. 113, 527 (1959)
    Elab = Ecm + KEf*Mn/Mf + 2.0*sqrt(KEf*Mn/Mf*Ecm)*cosCmTheta; // Eq.(1) Terrell, Phys. Rev. 113, 527 (1959)
    
    //    labTheta = acos ( (sqrt(Ecm)*cosCmTheta+sqrt(KEf/Af))/sqrt(Elab) ); // Eq.(34) Terrell, Id.
    //    labTheta = acos ( (sqrt(Ecm)*cosCmTheta+sqrt(KEf*Mn/Mf))/sqrt(Elab) ); // Eq.(34) Terrell, Id.
    // equivalently, labTheta = atan( sin(cmTheta)/(cosCmTheta+sqrt(KEf/Af/Ecm)))
    // cout << atan( sin(cmTheta)/(cosCmTheta+sqrt(KEf/Af/Ecm))) << endl;
    labTheta = atan( sin(cmTheta)/(cosCmTheta+sqrt(KEf*Mn/Mf/Ecm)));
    
    if (labTheta<0) labTheta = PI+labTheta; // TO REVISIT ........................................................................
    //		if (index!=0) labTheta = PI-labTheta; // emission from heavy fragment instead....................................<<<<<<<<<<<<<<<<<<<<<<<<
    
    //		cout << "old: " << Ecm << " " << Elab << " " << cmTheta << " " << labTheta << endl;
    
    // record neutron event in laboratory frame
    event->emissions[0].Vlab[i][0] = vlx;
    event->emissions[0].Vlab[i][1] = vly;
    event->emissions[0].Vlab[i][2] = vlz;

    
    // recoil of the fragment --------------------------------------------------
    
    KEf = (Mf-Mn)/Mf*KEf + Mn/(Mf-Mn)*Ecm-2.0*sqrt(Ecm*KEf/Mf*Mn)*cosCmTheta;
		
    // new fragment velocity components
    // (pf will be multiplied later by new Mf to get back to momentum)
    pfx2 = pfx2/Mf - pnx2/(Mf-Mn);
    pfy2 = pfy2/Mf - pny2/(Mf-Mn);
    pfz2 = pfz2/Mf - pnz2/(Mf-Mn);
    
    Af--; // remove 1 neutron from evaporating fission fragment
    Mf = Af*amuMeV+mass_excess(Zf,Af);
    
    pfx2 = pfx2*Mf;
    pfy2 = pfy2*Mf;
    pfz2 = pfz2*Mf;
    
    pf = sqrt (pfx2*pfx2+pfy2*pfy2+pfz2*pfz2);
    //    cout << KEf << " " << pf << " " << pf*pf/(2.0*Mf) << endl;

		// recalculate Elab based on lab. neutron vector
		Elab = (vlx*vlx+vly*vly+vlz*vlz)*0.5*Mn;
		
    event->emissions[0].cmEnergies[i]  = Ecm;
    event->emissions[0].labEnergies[i] = Elab;
    
    if (index==0) {
      event->emissions[0].cmAngles[i]    = cmTheta;
      event->emissions[0].labAngles[i]   = labTheta;
    } else {
      event->emissions[0].cmAngles[i]    = PI-cmTheta;
      event->emissions[0].labAngles[i]   = PI-labTheta;
    }
    
    
  } // end loop over neutron emissions

	
	// recover pf in original (x,y,z) lab frame
	pfx = sinthetaF*cos(phiF)*pfx2+costhetaF*cos(phiF)*pfy2-sin(phiF)*pfz2;
	pfy = sinthetaF*sin(phiF)*pfx2+costhetaF*sin(phiF)*pfy2+cos(phiF)*pfz2;
	pfz = costhetaF*pfx2-sinthetaF*pfy2;
	
	pf = sqrt(pfx*pfx+pfy*pfy+pfz*pfz);

	// record post-neutron emission fragment momentum in the lab frame
	event->postMomentum[0] = pfx;
	event->postMomentum[1] = pfy;
	event->postMomentum[2] = pfz;
	
	// gammas --------------------------------------------------------------------
	
  event->emissions[1].multiplicity = gammaMultiplicity;
  for (int i=0; i<gammaMultiplicity; i++) {
		
		// direction of gamma emission chosen isotropic in c.m. of fission fragment
		phi = twopi*genrand_real1();
		cosCmTheta = 2.0*genrand_real1()-1.0;
		cmTheta = acos(cosCmTheta);

		double Eg = gammaEnergies[i];
		
		// photon momentum in c.m. of the fragment
		double pgx2 = Eg*cosCmTheta;
		double pgy2 = Eg*sin(cmTheta)*cos(phi);
		double pgz2 = Eg*sin(cmTheta)*sin(phi);
		
		// record photon event in c.m. of fragment
		event->emissions[1].Vcm[i][0] = pgx2;
		event->emissions[1].Vcm[i][1] = pgy2;
		event->emissions[1].Vcm[i][2] = pgz2;
		
		double b = pf/Mf;
		double g = 1./sqrt(1.0-b*b);

		double Egl;

//		if (index!=0) cmTheta=PI-cmTheta;
		Egl= g*Eg*(1.0-b*cos(cmTheta)); // Doppler shift
		
    event->emissions[1].cmEnergies[i]  = gammaEnergies[i];
//		event->emissions[1].labEnergies[i] = gammaEnergies[i]; // same as center-of-mass energies FOR NOW
		event->emissions[1].labEnergies[i] = Egl;

//		// Boost photon in c.m. of fragment --> lab
//		double vlx2 = pgx2 + pfx2/Mf;
//		double vly2 = pgy2 + pfy2/Mf;
//		double vlz2 = pgz2 + pfz2/Mf;
//		
//		// revert photon vector to lab. frame
//		double vlx = sinthetaF*cos(phiF)*vlx2+costhetaF*cos(phiF)*vly2-sin(phiF)*vlz2;
//		double vly = sinthetaF*sin(phiF)*vlx2+costhetaF*sin(phiF)*vly2+cos(phiF)*vlz2;
//		double vlz = costhetaF*vlx2-sinthetaF*vly2;
//		double vl = sqrt(vlx*vlx+vly*vly+vlz*vlz);
		
		// FOR NOW, keep same orientation but change gamma-ray energy only
		double pgx = Egl*cosCmTheta;
		double pgy = Egl*sin(cmTheta)*cos(phi);
		double pgz = Egl*sin(cmTheta)*sin(phi);
		event->emissions[1].Vlab[i][0] = pgx;
		event->emissions[1].Vlab[i][1] = pgy;
		event->emissions[1].Vlab[i][2] = pgz;
		
    if (index==0) {
			event->emissions[1].cmAngles[i]    = cmTheta;
			event->emissions[1].labAngles[i]   = labTheta;
		} else {
			event->emissions[1].cmAngles[i]    = PI-cmTheta;
			event->emissions[1].labAngles[i]   = PI-cmTheta; // same as center-of-mass FOR NOW
		}
    
    event->emissions[1].transitionTypes[i] = gammaTypes[i];
  }
  
  event->emissions[2].multiplicity = icMultiplicity;
  
  for (int i=0; i<icMultiplicity; i++) { event->emissions[2].cmEnergies[i] = icEnergies[i]; }
  for (int i=0; i<icMultiplicity; i++) { event->emissions[2].labEnergies[i] = icEnergies[i]; } // same as center-of-mass energies for now
  
  //  cout << neutronMultiplicity << " " << gammaMultiplicity << endl;
  //  exit(0);
  //
  
}



/*******************************************************************************
 Write Monte Carlo histories into file.
 ******************************************************************************/
void FissionEvents::writeHistories (string outputFilename, int ZAIDt, double Einc, int nevents, double alphaI) {

  int ip=0, np=1;
  fragmentEventType *event;
  
  std::ostream* fp = &cout;
  std::ofstream fout;
  
  string outputFilename2;
  std::ostream* fp2;
  std::ofstream fout2;
  
#ifdef MPIRUN
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
#endif

  outputFilename2 = "histories-vectors.CGMF"+fileExt;
	
  if (ip==0 && outputFilename!="") {
    remove (&(WORKDIR+outputFilename)[0]);
    remove (&(WORKDIR+outputFilename2)[0]);
    fout.open(&(WORKDIR+outputFilename)[0], ios::out);
    fout << "# " << ZAIDt << " " << Einc << " " << nevents << " " << alphaI << "\n";
    fout.close();
    fout2.open(&(WORKDIR+outputFilename2)[0], ios::out);
    fout2 << "# " << ZAIDt << " " << Einc << " " << nevents << " " << alphaI << "\n";
    fout2.close();
  }
  
  
  for (int iproc=0; iproc<np; iproc++) { // loop over processors
    
    if (iproc==ip) {
      
      if (outputFilename!="") {
        fout.open(&(WORKDIR+outputFilename)[0], ios::out | ios::app );
        fp = &fout;
        
        fout2.open(&(WORKDIR+outputFilename2)[0], ios::out | ios::app);
        fp2 = &fout2;
      }
      
      *fp << fixed;
      *fp2 << fixed;
      
      if (lfCounter!=hfCounter) { cerr << "LF and HF counters mismatch!\n"; exit(-1); }
      
      for (int i=0; i<=lfCounter; i++) { // loop over all fission events
        
        numberPFN = Ac - lightFragments[i].A - heavyFragments[i].A;
        
        for (int k=0; k<2; k++) { // loop over light and heavy fragments
          
          if (k==0) { event = &lightFragments[i]; } else { event = &heavyFragments[i]; }
          
          *fp << event->A << " " << event->Z << " " << setprecision(3) << event->Ui << " " << setprecision(1) << event->Ji << " "
          << event->Pi << " " << setprecision(3) << event->KE << " "  << event->emissions[0].multiplicity << " "
          << event->emissions[1].multiplicity << " " << event->emissions[2].multiplicity;

          *fp2 << event->A << " " << event->Z << " " << setprecision(3) << event->Ui << " " << setprecision(1) << event->Ji << " "
          << event->Pi << " " << setprecision(3) << event->KE << " "  << event->emissions[0].multiplicity << " "
          << event->emissions[1].multiplicity << " " << event->emissions[2].multiplicity << endl;

	  *fp2 << event->preMomentum[0] << " " << event->preMomentum[1] << " " << event->preMomentum[2] << " ";
	  *fp2 << event->postMomentum[0] << " " << event->postMomentum[1] << " " << event->postMomentum[2] << "\n";
					
          // write out type of transitions
          for (int m=0; m<=1; m++) { // loop over emission type (n, g, or IC)
            if (event->emissions[m].multiplicity>0) {
              for (int j=0; j<event->emissions[m].multiplicity; j++) {
                *fp << " ";
                if (m==0) { // neutrons
                  *fp << "0";
                } else if (m==1) { // gammas
                  *fp << event->emissions[1].transitionTypes[j];
                } else { // IC
                  *fp << "4";
                }
              }
            }
          }
          *fp << "\n";
          
					*fp << setprecision(3);
          *fp2 << setprecision(3);
          for (int m=0; m<=1; m++) { // loop over emission type (n, g, or IC)
            if (event->emissions[m].multiplicity>0) {
              for (int j=0; j<event->emissions[m].multiplicity; j++) {
								*fp << event->emissions[m].cmEnergies[j] << " " << event->emissions[m].cmAngles[j] << " "
                  << event->emissions[m].labEnergies[j] << " " << event->emissions[m].labAngles[j] << " ";
                *fp2 <<
                  event->emissions[m].Vcm[j][0] << " " << event->emissions[m].Vcm[j][1] << " " << event->emissions[m].Vcm[j][2] << " "
									<< event->emissions[m].cmEnergies[j] << " " <<
                  event->emissions[m].Vlab[j][0] << " " << event->emissions[m].Vlab[j][1] << " " << event->emissions[m].Vlab[j][2] << " "
									<< event->emissions[m].labEnergies[j] << " ";
              }
              *fp << endl;
              *fp2 << endl;
            }
          }
          
        } // end loop over light and heavy fragments
        
        if (numberPFN>0) {
          for (int m=0; m<numberPFN; m++) {
            *fp << pfnInfo[i][m] << " ";
            *fp2 << pfnInfo[i][m] << " ";
          }
          *fp << endl;
          *fp2 << endl;
          delete [] pfnInfo[i];
        }
        
      } // end loop over fission events
      
      if (fout.is_open()) { fout.close(); }
      if (fout2.is_open()) { fout2.close(); }
    }
#ifdef MPIRUN
    MPI_Barrier( MPI_COMM_WORLD ) ;
#endif
  }
  
}


/*******************************************************************************
 Read Monte Carlo fission events from history file produced earlier. The number
 of events read is given as input as "numberEvents".
 ******************************************************************************/
void FissionEvents::readHistories (string inputFilename, int *numberEvents, int *ZAIDt, double *Einc, double *alphaI) {
  
  int ip=0;
  int Zl=0, Al=0, Pil=0, Zh=0, Ah=0, Pih=0;
  double Ul=0.0, KEl=0.0, Uh=0.0, KEh=0.0, Jl=0.0, Jh=0.0;
  fragmentEventType *event;
  double *pfnEvent = new double [5];
  
  ifstream histories;
  
  int totalMultiplicity;
  int transitionType;
  
#ifdef MPIRUN
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
#endif
  
  int ZAIDc;
  if (*Einc>0.0) {
    ZAIDc = *ZAIDt+1;
  } else {
    ZAIDc = *ZAIDt;
  }
  int Ac = ZAIDc%1000;
  
  if (ip!=0) { cerr << "Only one processor required when reading history file!\n"; exit(-1); }
  
  numberFissionEvents = *numberEvents;
  lfCounter = *numberEvents;
  hfCounter = *numberEvents;
  
  histories.open(&inputFilename[0], ios::in);
  
  if (!histories) { cerr << "[readHistories] Cannot find history file: " << inputFilename << "\n"; exit(-1); }
  
  // read first information line
  string dummy;
  int savedNumberEvents;
  histories >> dummy >> *ZAIDt >> *Einc >> savedNumberEvents >> *alphaI;
  if (savedNumberEvents<*numberEvents) *numberEvents = savedNumberEvents;
  
  int tenPercent = (int) (0.1* *numberEvents);
  
  int c=-1;
  while (!histories.eof()) {
    
    c++; // increment event counter
    if (c>=*numberEvents) {
      printf ("Total number of events read in: %i\n", c);
      numberFissionEvents = c;
      lfCounter = c;
      hfCounter = c;
      histories.close();
      return;
    }
    
    if (*numberEvents>=100)
      if ((c+1)%tenPercent==0) cout << " Reading fission events: " << double(c+1)/double(*numberEvents)*100 << "%\n";
    
    for (int f=0; f<=1; f++) { // loop over light and heavy fragments -------------------
      
      if (f==0) { event = &lightFragments[c]; } else { event = &heavyFragments[c]; pfnEvent = preFissionNeutrons[c]; }
      
      histories >> event->A >> event->Z >> event->Ui >> event->Ji >> event->Pi >> event->KE >>
      event->emissions[0].multiplicity >> event->emissions[1].multiplicity >> event->emissions[2].multiplicity;
      
      // for Ionel's history files
      //      histories >> event->Z >> event->A >> event->Ui >> event->Ji >> event->Pi >> event->KE >>
      //      event->emissions[0].multiplicity >> event->emissions[1].multiplicity >> event->emissions[2].multiplicity;
      
      if (f==0) {
        Zl=event->Z; Al=event->A; Ul=event->Ui; Jl=event->Ji; Pil=event->Pi; KEl=event->KE;
      } else {
        Zh=event->Z; Ah=event->A; Uh=event->Ui; Jh=event->Ji; Pih=event->Pi; KEh=event->KE;
        addFragments (Zl, Al, Ul, Jl, Pil, KEl, Zh, Ah, Uh, Jh, Pih, KEh);
      }
      
      // Read in type of transitions
      totalMultiplicity = event->emissions[0].multiplicity + event->emissions[1].multiplicity; // + event->emissions[2].multiplicity;
      //      totalMultiplicity = event->emissions[0].multiplicity + event->emissions[1].multiplicity + event->emissions[2].multiplicity; // Ionel's old files
      int ig=-1;
      for (int k=0; k<totalMultiplicity; k++) {
        histories >> transitionType;
        if (transitionType>=1 && transitionType<=3) event->emissions[1].transitionTypes[++ig] = transitionType;
      }
      
      for (int m=0; m<2; m++) { // loop over emission type (n, g, or IC) -- no IC for now
        //			for (int m=0; m<3; m++) { // loop over emission type (n, g, or IC) // Ionel's old files
        for (int j=0; j<event->emissions[m].multiplicity; j++) {
          //          histories >> event->emissions[m].labEnergies[j] >> event->emissions[m].labAngles[j]; // Ionel's old history files
          histories >> event->emissions[m].cmEnergies[j] >> event->emissions[m].cmAngles[j] >>
          event->emissions[m].labEnergies[j] >> event->emissions[m].labAngles[j];
        }
      }
      
      if (f>0) {
        std::fill_n (pfnEvent, 5, 0.0);
        if (Al+Ah<Ac) { // pre-fission neutrons emitted
          numberPFN = Ac-Al-Ah;
          for (int j=0; j<numberPFN; j++) histories >> pfnEvent[j];
        }
        Al=0;
        Ah=0;
      }
      
    } // end loop over light and heavy fragments ------------------------------------------
    
    //		if (c>0) { exit(0); }
    
  } // end reading file
  histories.close();
  
  if (c<*numberEvents) {
    numberFissionEvents = c;
    lfCounter = c;
    hfCounter = c;
    printf ("Total number of events read in: %i\n", c);
  }
  
}


/*******************************************************************************
 Analyze all Monte Carlo histories, and compute average quantities, distributions
 and correlations.
 ******************************************************************************/
void FissionEvents::analyzeResults (int ZAIDc) {
  
  emittedParticleType *particles;
  
  fragmentEventType *event;
  double *pfnEvent = new double [5];
  
  int A;
  int Al=0, Ah=0;
  double KEl=0.0, KEh=0.0, TKE=0.0;
  int iTKE=0, iEcm, iElab;
  double Ecm, Elab;
  double cmTheta, labTheta;
  int multiplicity;
  
  int lfNu=0, lfNuGamma=0, hfNu=0, hfNuGamma=0;
  
  struct tmpParticleType {
    int nuTot;
    double EcmTot, ElabTot;
  };
  
  tmpParticleType tmpParticles [3]; // [0] neutrons [1] gammas [2] internal conversion
  
  numberPFN = 0;
  int Acn = ZAIDc%1000; // Mass of original compound nucleus
  
  //  for (int n=0; n<lfCounter+1; n++) { // loop over all fission events ------------------
  for (int n=0; n<lfCounter; n++) { // loop over all fission events ------------------
    
    for (int i=0; i<=2; i++) {
      tmpParticles[i].nuTot=0;
      tmpParticles[i].EcmTot=0.0;
      tmpParticles[i].ElabTot=0.0;
    }
    
    // pre-fission neutron energies (max. 5)
    int pfnNumber = Acn-lightFragments[n].A-heavyFragments[n].A; // number of pre-fission neutrons for this event
    numberPFN += pfnNumber;
    
    pfnEvent = preFissionNeutrons[n]; // energies
    for (int m=0; m<pfnNumber; m++) {
      aveEnergyPFN += pfnEvent[m];
    }
    
    for (int index=0; index<=1; index++) { // loop over light and heavy fragments ---------------------
      
      if (index==0) {
        event     = &lightFragments[n];
        Al        = event->A;
        KEl       = event->KE;
        lfNu      = event->emissions[0].multiplicity;
        lfNuGamma = event->emissions[1].multiplicity;
      } else {
        event     = &heavyFragments[n];
        Ah        = event->A;
        KEh       = event->KE;
        hfNu      = event->emissions[0].multiplicity;
        hfNuGamma = event->emissions[1].multiplicity;
        
        TKE = KEl + KEh;
        iTKE = findEnergyIndex (TKE, excitationEnergyGrid);
      }
      
      // to clarify coding only...
      A = event->A;
      //      Z = event->Z;
      
      // cut ----------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      /*			if (A>Asym) {
       if (A<125 || A>140) { continue; }
       } else {
       if (A<96 || A>111) { continue; }
       }
       */
      
      YApost[event->A-event->emissions[0].multiplicity]++;
      
      for (int type=0; type<=1; type++) { // loop over neutrons, gammas (NOT FOR NOW: and IC ) ---
        
        multiplicity = event->emissions[type].multiplicity;
        
        if (type==0) {
          particles = &neutrons;
        } else if (type==1) {
          particles = &gammas;
        } else {
          particles = &ic;
        }
        
        particles->nuA[A] += multiplicity;
        
        particles->counts += multiplicity;
        particles->nubar  += multiplicity;
        
        if (index==0) {
          particles->lfCounts += multiplicity;
          particles->lfNubar  += multiplicity;
          particles->lfPnu[multiplicity]++;
        } else {
          particles->hfCounts += multiplicity;
          particles->hfNubar  += multiplicity;
          particles->hfPnu[multiplicity]++;
        }
        
        for (int i=0; i<multiplicity; i++) {
          
          Ecm  = event->emissions[type].cmEnergies[i];
          Elab = event->emissions[type].labEnergies[i];
          
          cmTheta = event->emissions[type].cmAngles[i];
          labTheta = event->emissions[type].labAngles[i];
          
          //          if ((Elab<0.1 || Elab>6.0) && (type==1)) continue; // energy cut for gamma rays only
          
          iEcm = findEnergyIndex (Ecm, spectrumEnergyGrid);
          iElab = findEnergyIndex(Elab, spectrumEnergyGrid);
          
          particles->Ecm     += Ecm;
          particles->EcmA[A] += Ecm;
          
          particles->cmSpectrum[iEcm]++;
          particles->cmSpectrumA[A][iEcm]++;
          particles->cmHistogramA[A][min(int(Ecm/0.05),300)]++;
          
          //	Changed to account for spectra for first, second, third, etc, neutron spectra instead !!!
          //				if (multiplicity<10) particles->cmExclusiveSpectra[multiplicity][iEcm]++;
          if (multiplicity<10) particles->cmExclusiveSpectra[i+1][iEcm]++;
          
          if (i<4) { particles->cmSpectrumNuA[i+1][A][iEcm]++; }
          if (i<4) { particles->cmHistogramNuA[i+1][A][min(int(Ecm/0.05),300)]++; }
          
          particles->Elab  += Elab;
          particles->labSpectrum[iElab]++;
          
          if (particles->particleType=='g') {
            if (event->emissions[1].transitionTypes[i]==1) {
              particles->ccLabSpectrum[iElab]++;
              particles->ccLabSpectrumA[A][iElab]++;
            } else if (event->emissions[1].transitionTypes[i]==2) {
              particles->cdLabSpectrum[iElab]++;
              particles->cdLabSpectrumA[A][iElab]++;
            } else if (event->emissions[1].transitionTypes[i]==3) {
              particles->ddLabSpectrum[iElab]++;
              particles->ddLabSpectrumA[A][iElab]++;
            }
          }
          
          particles->ElabA[A] += Elab;
          
          if (particles->particleType=='n') {
            YThetacm[int(cmTheta/dTheta/PI*180.0)]++;
            YThetalab[int(labTheta/dTheta/PI*180.0)]++;
          }
          
          if (index==0) {
            particles->lfEcm  += Ecm;
            particles->lfElab += Elab;
            particles->cmLightFragmentSpectrum[iEcm]++;
            particles->labLightFragmentSpectrum[iElab]++;
          } else {
            particles->hfEcm  += Ecm;
            particles->hfElab += Elab;
            particles->cmHeavyFragmentSpectrum[iEcm]++;
            particles->labHeavyFragmentSpectrum[iElab]++;
          }
          
          tmpParticles[type].nuTot++;
          tmpParticles[type].EcmTot+=Ecm;
          tmpParticles[type].ElabTot+=Elab;
          
        } // end loop over particle multiplicity
        
      } // end loop over neutrons and gammas
      
    } // end loop over light and heavy fragments --------------------------------------
    
    for (int type=0; type<=2; type++) { // loop over neutrons, gammas and IC ---
      
      if (type==0) {
        particles = &neutrons;
      } else if (type==1) {
        particles = &gammas;
      } else {
        particles = &ic;
      }
      
      //		particles->nubar += tmpParticles[type].nuTot;
      particles->Pnu[tmpParticles[type].nuTot]++;
      particles->nuTKE[iTKE]   += tmpParticles[type].nuTot;
      particles->EcmTKE[iTKE]  += tmpParticles[type].EcmTot;
      particles->ElabTKE[iTKE] += tmpParticles[type].ElabTot;
      
      if (type==0) {
        particles->nubarATKE[Al][iTKE] += lfNu;
        particles->nubarATKE[Ah][iTKE] += hfNu;
      } else if (type==1) {
        particles->nubarATKE[Al][iTKE] += lfNuGamma;
        particles->nubarATKE[Ah][iTKE] += hfNuGamma;
      }
      
      particles->countsATKE[Al][iTKE]++;
      particles->countsATKE[Ah][iTKE]++;
      
    }
    
    totalGammaEnergyTKE[iTKE] += tmpParticles[1].ElabTot;
    totalGammaEnergyA[Al] += tmpParticles[1].ElabTot;
    
    // neutron-gamma correlations
    gammaMultiplicityNu[tmpParticles[0].nuTot]+=tmpParticles[1].nuTot;
    totalGammaEnergyNu[tmpParticles[0].nuTot]+=tmpParticles[1].ElabTot;
    countsNu[tmpParticles[0].nuTot]++;
    
    //if (n>10) exit(0);
    
  } // end loop over all fission events -----------------------------------------------
  
  
#ifdef MPIRUN
  
  double dummy;
  int idummy;
  double *array;
  
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  
  MPI_Reduce (&numberPFN, &idummy, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
  numberPFN = idummy;
  
  for (int type=0; type<=1; type++) { // loop over neutrons, gammas and IC
    
    if (type==0) {
      particles=&neutrons;
    } else if (type==1) {
      particles=&gammas;
    } else {
      particles=&ic;
    }
    
    MPI_Reduce (&particles->nubar, &dummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->nubar = dummy;
    MPI_Reduce (&particles->lfNubar, &dummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->lfNubar = dummy;
    MPI_Reduce (&particles->hfNubar, &dummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->hfNubar = dummy;
    
    MPI_Reduce (&particles->counts, &idummy, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->counts = idummy;
    MPI_Reduce (&particles->lfCounts, &idummy, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->lfCounts = idummy;
    MPI_Reduce (&particles->hfCounts, &idummy, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->hfCounts = idummy;
    
    MPI_Reduce (&particles->Ecm, &dummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->Ecm = dummy;
    MPI_Reduce (&particles->lfEcm, &dummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->lfEcm = dummy;
    MPI_Reduce (&particles->hfEcm, &dummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->hfEcm = dummy;
    
    MPI_Reduce (&particles->Elab, &dummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->Elab = dummy;
    MPI_Reduce (&particles->lfElab, &dummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->lfElab = dummy;
    MPI_Reduce (&particles->hfElab, &dummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    particles->hfElab = dummy;
    
    // f(nu)
    array = new double [NUMMULT];
    MPI_Reduce (particles->Pnu, array, NUMMULT, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMMULT; i++) { particles->Pnu[i] = array[i]; }
    MPI_Reduce (particles->lfPnu, array, NUMMULT, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMMULT; i++) { particles->lfPnu[i] = array[i]; }
    MPI_Reduce (particles->hfPnu, array, NUMMULT, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMMULT; i++) { particles->hfPnu[i] = array[i]; }
    delete [] array;
    
    // f(Eout)
    array = new double [NUMBER_SPECTRUM_ENERGY_GRID];
    MPI_Reduce (particles->cmSpectrum, array, NUMBER_SPECTRUM_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) { particles->cmSpectrum[i] = array[i]; }
    MPI_Reduce (particles->cmLightFragmentSpectrum, array, NUMBER_SPECTRUM_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) { particles->cmLightFragmentSpectrum[i] = array[i]; }
    MPI_Reduce (particles->cmHeavyFragmentSpectrum, array, NUMBER_SPECTRUM_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) { particles->cmHeavyFragmentSpectrum[i] = array[i]; }
    MPI_Reduce (particles->labSpectrum, array, NUMBER_SPECTRUM_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) { particles->labSpectrum[i] = array[i]; }
    MPI_Reduce (particles->labLightFragmentSpectrum, array, NUMBER_SPECTRUM_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) { particles->labLightFragmentSpectrum[i] = array[i]; }
    MPI_Reduce (particles->labHeavyFragmentSpectrum, array, NUMBER_SPECTRUM_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) { particles->labHeavyFragmentSpectrum[i] = array[i]; }
    
    MPI_Reduce (particles->ccLabSpectrum, array, NUMBER_SPECTRUM_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) { particles->ccLabSpectrum[i] = array[i]; }
    MPI_Reduce (particles->cdLabSpectrum, array, NUMBER_SPECTRUM_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) { particles->cdLabSpectrum[i] = array[i]; }
    MPI_Reduce (particles->ddLabSpectrum, array, NUMBER_SPECTRUM_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) { particles->ddLabSpectrum[i] = array[i]; }
    
    delete [] array;
    
    // f(A)
    array = new double [NUMA];
    MPI_Reduce (particles->nuA, array, NUMA, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMA; i++) { particles->nuA[i] = array[i]; }
    MPI_Reduce (particles->EcmA, array, NUMA, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMA; i++) { particles->EcmA[i] = array[i]; }
    delete [] array;
    
    // f(TKE)
    array = new double [NUMBER_EXCITATION_ENERGY_GRID];
    MPI_Reduce (particles->nuTKE, array, NUMBER_EXCITATION_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) { particles->nuTKE[i] = array[i]; }
    MPI_Reduce (particles->EcmTKE, array, NUMBER_EXCITATION_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) { particles->EcmTKE[i] = array[i]; }
    MPI_Reduce (particles->ElabTKE, array, NUMBER_EXCITATION_ENERGY_GRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) { particles->ElabTKE[i] = array[i]; }
    delete [] array;
    
  }
  
#endif
  
  
}



/*******************************************************************************
 Compute the final results for the fission fragments, the emitted neutrons and
 emitted gamma rays.
 ******************************************************************************/
void FissionEvents::computeFinalResults () {
  
  double sum, suml, sumh;
  double de [NUMBER_EXCITATION_ENERGY_GRID];
  
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID-1; i++) {
    de[i] = excitationEnergyGrid[i+1]-excitationEnergyGrid[i];
  }
  
  computeFinalResults (&neutrons);
  computeFinalResults (&gammas);
  
  // <Ui>(A), <nu>(A), <Ecm>(A), ...
  // (IMPORTANT: to be computed before renormalizing Y(A) !)
  for (int i=0; i<NUMA; i++) {
    if (YA[i]!=0.0) {
      UiA[i]               /= YA[i];
      initialSpinA[i]      /= YA[i];
      totalGammaEnergyA[i] /= YA[i];
    }
  }
  
  // Proba (spin) in light and heavy fragments
  sum=0.0;
  for (int i=0; i<100; i++) sum += lfPspin[i];
  for (int i=0; i<100; i++) lfPspin[i] /= sum;
  sum=0.0;
  for (int i=0; i<100; i++) sum += hfPspin[i];
  for (int i=0; i<100; i++) hfPspin[i] /= sum;
  
  // Yields Y(A) (normalized to 1)
  sum=0.0;
  for (int i=0; i<NUMA; i++) sum += YA[i];
  for (int i=0; i<NUMA; i++) YA[i] /= sum;
  
  // post-neutron emission yields Y(A)
  sum=0.0;
  for (int i=0; i<NUMA; i++) sum += YApost[i];
  for (int i=0; i<NUMA; i++) YApost[i] /= sum;
  
  // Yields Y(Z) (normalized to 1)
  sum=0.0;
  for (int i=0; i<NUMZ; i++) sum += YZ[i];
  for (int i=0; i<NUMZ; i++) YZ[i] /= sum;
  
  // <Eg^tot>(TKE)
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    if (YTKE[i]!=0) {
      totalGammaEnergyTKE [i] /= YTKE[i];
      initialSpinTKE[i] /= (2.0*YTKE[i]); // factor 2.0 because both LF and HF spins are counted
    }
  }
  
  // Yields Y(KEl), Y(KEh), Y(TKE) (normalized to 1)
  sum=0.0;
  suml = 0.0;
  sumh = 0.0;
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    sum  += ( YTKE[i] * de[i] );
    suml += ( YKEl[i] * de[i] );
    sumh += ( YKEh[i] * de[i] );
  }
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    YTKE[i] /= sum;
    YKEl[i] /= suml;
    YKEh[i] /= sumh;
  }
  
  // Yields Y(Ul), Y(Uh), Y(Utot) (normalized to 1)
  sum=0.0;
  suml = 0.0;
  sumh = 0.0;
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID-1; i++) {
    sum  += YUtot[i];
    suml += YUl[i];
    sumh += YUh[i];
  }
  
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    YUtot[i] /= (sum*de[i]);
    YUl[i]   /= (suml*de[i]);
    YUh[i]   /= (sumh*de[i]);
  }
  
  for (int j=0; j<NUMA; j++) {
    sum=0.0;
    for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) { sum += YUiA[j][i]; }
    if (sum!=0.0) { for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) { YUiA[j][i] /= sum; } }
  }
  
  // Average yields
  
  aveMass=aveMassLF=aveMassHF=0.0;
  double x;
  for (int i=0; i<NUMA; i++) {
    x=YA[i]*i;
    aveMass+=x;
    if (i<Asym) { aveMassLF+=x; } else { aveMassHF+=x; }
  }
  aveMassLF*=2.0;
  aveMassHF*=2.0;
  
  aveCharge=aveChargeLF=aveChargeHF=0.0;
  for (int i=0; i<NUMZ; i++) {
    x = YZ[i]*i;
    aveCharge+=x;
    if (i<Zsym) { aveChargeLF+=x; } else { aveChargeHF+=x; }
  }
  aveChargeLF*=2.0;
  aveChargeHF*=2.0;
  
  aveTKE=0.0;
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    aveTKE += (YTKE[i]*excitationEnergyGrid[i]*de[i]);
    aveKEl += (YKEl[i]*excitationEnergyGrid[i]*de[i]);
    aveKEh += (YKEh[i]*excitationEnergyGrid[i]*de[i]);
  }
  
  // <Eg>(nu) and <Ng>(nu)
  for (int i=0; i<20; i++) {
    if (countsNu[i]!=0) {
      totalGammaEnergyNu[i]  /= countsNu[i];
      gammaMultiplicityNu[i] /= countsNu[i];
    }
  }
  
  // renormalize Angular distributions
  sum=0.0;
  for (int i=0; i<NUMANGLES; i++) { sum += (YThetacm[i]*sin(i*dTheta/180.0*PI)*dTheta); }
  for (int i=0; i<NUMANGLES; i++) { YThetacm[i] /= sum; }
  
}




/*******************************************************************************
 Compute the final results for a specific type of emitted particles, neutron
 or gamma.
 ******************************************************************************/
void FissionEvents::computeFinalResults (emittedParticleType * particles) {
  
  if (numberFissionEvents==0) {
    cerr << "No fission occurred! Stopping the program.\n";
    exit(-1);
  }
  
  int n = NUMBER_SPECTRUM_ENERGY_GRID;
  double sum;
  
  // Transform histograms into spectra
  for (int i=0; i<n; i++) {
    particles->labSpectrum[i] /= spectrumEnergySteps[i];
    particles->cmSpectrum[i]  /= spectrumEnergySteps[i];
    particles->cmLightFragmentSpectrum[i] /= spectrumEnergySteps[i];
    particles->cmHeavyFragmentSpectrum[i] /= spectrumEnergySteps[i];
    particles->ccLabSpectrum[i] /= spectrumEnergySteps[i];
    particles->cdLabSpectrum[i] /= spectrumEnergySteps[i];
    particles->ddLabSpectrum[i] /= spectrumEnergySteps[i];
    for (int j=1; j<5; j++) particles->cmExclusiveSpectra[j][i] /= spectrumEnergySteps[i];
    for (int j=1; j<NUMA; j++) particles->cmSpectrumA[j][i] /= spectrumEnergySteps[i];
    for (int k=0; k<5; k++) {
      for (int j=1; j<NUMA; j++) particles->cmSpectrumNuA[k][j][i] /= spectrumEnergySteps[i];
    }
  }
  
  // renormalize center-of-mass spectrum
  sum=0.0;
  for (int i=0; i<n; i++) sum += (particles->cmSpectrum[i]*spectrumEnergySteps[i]);
  if (sum!=0.0) { for (int i=0; i<n; i++) particles->cmSpectrum[i] /= sum; }
  
  // renormalize laboratory spectrum
  sum=0.0;
  for (int i=0; i<n; i++) sum += (particles->labSpectrum[i]*spectrumEnergySteps[i]);
  if (sum!=0.0) { for (int i=0; i<n; i++) particles->labSpectrum[i] /= sum; }
  
  // renormalize c-c, c-d, and d-d lab. spectra
  sum=0.0;
  for (int i=0; i<n; i++) sum += (particles->ccLabSpectrum[i]*spectrumEnergySteps[i]);
  if (sum!=0.0) { for (int i=0; i<n; i++) particles->ccLabSpectrum[i] /= sum; }
  sum=0.0;
  for (int i=0; i<n; i++) sum += (particles->cdLabSpectrum[i]*spectrumEnergySteps[i]);
  if (sum!=0.0) { for (int i=0; i<n; i++) particles->cdLabSpectrum[i] /= sum; }
  sum=0.0;
  for (int i=0; i<n; i++) sum += (particles->ddLabSpectrum[i]*spectrumEnergySteps[i]);
  if (sum!=0.0) { for (int i=0; i<n; i++) particles->ddLabSpectrum[i] /= sum; }
  
  
  // renormalize average c.m. and lab spectra for Light Fragments
  sum=0.0;
  for (int i=0; i<n; i++) sum += (particles->cmLightFragmentSpectrum[i]*spectrumEnergySteps[i]);
  if (sum!=0.0) { for (int i=0; i<n; i++) particles->labLightFragmentSpectrum[i] /= sum; }
  sum=0.0;
  for (int i=0; i<n; i++) sum += (particles->labLightFragmentSpectrum[i]*spectrumEnergySteps[i]);
  if (sum!=0.0) { for (int i=0; i<n; i++) particles->labLightFragmentSpectrum[i] /= sum; }
  
  // renormalize average c.m. and lab spectra for Heavy Fragments
  sum=0.0;
  for (int i=0; i<n; i++) sum += (particles->cmHeavyFragmentSpectrum[i]*spectrumEnergySteps[i]);
  if (sum!=0.0) { for (int i=0; i<n; i++) particles->labHeavyFragmentSpectrum[i] /= sum; }
  sum=0.0;
  for (int i=0; i<n; i++) sum += (particles->labHeavyFragmentSpectrum[i]*spectrumEnergySteps[i]);
  if (sum!=0.0) { for (int i=0; i<n; i++) particles->labHeavyFragmentSpectrum[i] /= sum; }
  
  
  for (int j=1; j<5; j++) {
    sum=0.0;
    for (int i=0; i<n; i++) sum += (particles->cmExclusiveSpectra[j][i]*spectrumEnergySteps[i]);
    if (sum!=0.0) { for (int i=0; i<n; i++) particles->cmExclusiveSpectra[j][i] /= sum; }
  }
  
  // renormalize mass-dependent center-of-mass spectra
  /*  for (int j=1; j<NUMA; j++) {
   sum=0.0;
   for (int i=0; i<n; i++) sum += (particles->cmSpectrumA[j][i]*spectrumEnergySteps[i]);
   if (sum!=0.0) { for (int i=0; i<n; i++) particles->cmSpectrumA[j][i] /= sum; }
   }
   */
  
  // multiplicity distributions
  
  sum=0.0;
  double suml=0.0;
  double sumh=0.0;
  for (int i=0; i<NUMMULT; i++) {
    sum  += particles->Pnu[i];
    suml += particles->lfPnu[i];
    sumh += particles->hfPnu[i];
  }
  
  for (int i=0; i<NUMMULT; i++) {
    particles->Pnu[i]   /= sum;
    particles->lfPnu[i] /= suml;
    particles->hfPnu[i] /= sumh;
  }
  
  // multiplicities
  
  particles->nubar   = 0.0;
  particles->lfNubar = 0.0;
  particles->hfNubar = 0.0;
  
  for (int i=0; i<NUMMULT; i++) {
    particles->nubar   += particles->Pnu[i]*i;
    particles->lfNubar += particles->lfPnu[i]*i;
    particles->hfNubar += particles->hfPnu[i]*i;
  }
  
  // energies
  
  particles->Ecm   /= particles->counts;
  particles->lfEcm /= particles->lfCounts;
  particles->hfEcm /= particles->hfCounts;
  
  particles->Elab   /= particles->counts;
  particles->lfElab /= particles->lfCounts;
  particles->hfElab /= particles->hfCounts;
  
  //-- <Ecm>(A), <nu>(A), ...
  for (int i=0; i<NUMA; i++) {
    if (particles->nuA[i]!=0.0) {
      particles->EcmA[i]  /= particles->nuA[i];
      particles->ElabA[i] /= particles->nuA[i];
      particles->nuA[i]   /= YA[i]; // renormalize <nu>(A) -- IMPORTANT: to be computed after <Ecm>(A) and <Elab>(A) !!!
    }
  }
  
  //-- <Ecm>(A), <nu>(A), ... w/ 1-6 MeV cut on Ecm
  for (int i=0; i<NUMA; i++) {
    if (particles->nuAcut[i]!=0.0) {
      particles->EcmAcut[i]  /= particles->nuAcut[i];
      particles->nuAcut[i]   /= YA[i]; // renormalize <nu>(A) -- IMPORTANT: to be computed after <Ecm>(A) and <Elab>(A) !!!
    }
  }
  
  //-- <Ecm>(TKE), <nu>(TKE), ...
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    if (YTKE[i]!=0.0 && particles->nuTKE[i]!=0.0) {
      particles->EcmTKE[i]  /= particles->nuTKE[i];
      particles->ElabTKE[i] /= particles->nuTKE[i];
      particles->nuTKE[i]   /= YTKE[i]; // renormalize <nu>(TKE) -- IMPORTANT: to be computed after <Ecm>(TKE) and <Elab>(TKE)
    }
  }
  
  
  
  
  // Particle-particle correlations
  
  if (particles->particleType=='n') { // WHY NOT GAMMA-GAMMA CORRELATIONS HERE?
    
    double **correlationMatrix;   correlationMatrix = new double * [NUMBER_OF_PRIVATE_GRID];
    //  if (ip==0) {
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) {
      correlationMatrix[i] = new double [NUMBER_OF_PRIVATE_GRID];
      for (int j=0; j<NUMBER_OF_PRIVATE_GRID; j++) {
        if (particles->energyCorrelations[i][i] != 0.0 && particles->energyCorrelations[j][j] != 0.0) {
          correlationMatrix[i][j] = particles->energyCorrelations[i][j] /
          sqrt ( particles->energyCorrelations[i][i] * particles->energyCorrelations[j][j] );
        } else {
          correlationMatrix[i][j] = 0.0;
        }
      }
    }
    
    for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) {
      for (int j=0; j<NUMBER_OF_PRIVATE_GRID; j++) {
        particles->energyCorrelations[i][j] = correlationMatrix[i][j];
      }
    }
    
    delete [] correlationMatrix;
    
  }
  
  // To infer <Ecm>(A) from <Ecm(A,TKE)> (following Vorobyev's question)
  if (particles->particleType=='n') {
    
    for (int i=0; i<NUMA; i++) {
      for (int j=0; j<NUMBER_EXCITATION_ENERGY_GRID; j++) {
        if (particles->nubarATKE[i][j]!=0) {
          particles->EcmATKE[i][j]   /= particles->nubarATKE[i][j];
          particles->nubarATKE[i][j] /= particles->countsATKE[i][j];
        }
      }
    }
    
    double EcmA [NUMA];
    int normA [NUMA];
    std::fill_n (normA, NUMA, 0);
    std::fill_n (EcmA, NUMA, 0.0);
    //    cout << particles->particleType << "\n";
    for (int i=0; i<NUMA; i++) {
      for (int j=0; j<NUMBER_EXCITATION_ENERGY_GRID; j++) {
        if (particles->countsATKE[i][j]!=0) {
          EcmA[i] += ( particles->EcmATKE[i][j]*particles->countsATKE[i][j] );
          normA[i] += particles->countsATKE[i][j];
        }
      }
      if (normA[i]!=0.0) { EcmA[i] /= float(normA[i]); }
      //    cout << i << " " << normA[i] << " " << EcmA[i] << "\n";
    }
    
  }
  
}


/*******************************************************************************
 Deallocate dynamical arrays in particles. (Called by destructor.)
 ******************************************************************************/
void FissionEvents::deleteParticles (emittedParticleType *particles) {
  
  delete [] particles->EcmA;
  delete [] particles->ElabA;
  delete [] particles->nuA;
  
  delete [] particles->EcmAcut;
  delete [] particles->nuAcut;
  
  delete [] particles->EcmTKE;
  delete [] particles->ElabTKE;
  delete [] particles->nuTKE;
  
  delete [] particles->cmSpectrum;
  delete [] particles->labSpectrum;
  
  delete [] particles->ccLabSpectrum;
  delete [] particles->cdLabSpectrum;
  delete [] particles->ddLabSpectrum;
  
  for (int i=0; i<NUMA; i++) {
    delete [] particles->ddLabSpectrumA[i];
    delete [] particles->cdLabSpectrumA[i];
    delete [] particles->ccLabSpectrumA[i];
  }
  delete [] particles->ddLabSpectrumA;
  delete [] particles->cdLabSpectrumA;
  delete [] particles->ccLabSpectrumA;
  
  delete [] particles->cmLightFragmentSpectrum;
  delete [] particles->cmHeavyFragmentSpectrum;
  
  delete [] particles->labLightFragmentSpectrum;
  delete [] particles->labHeavyFragmentSpectrum;
  
  for (int i=0; i<NUMA; i++) { delete [] particles->cmSpectrumA[i]; }
  delete [] particles->cmSpectrumA;
  
  for (int i=0; i<5; i++) {
    for (int j=0; j<NUMA; j++) { delete [] particles->cmSpectrumNuA[i][j]; }
    delete [] particles->cmSpectrumNuA[i];
  }
  delete [] particles->cmSpectrumNuA;
  
  for (int i=0; i<10; i++) { delete [] particles->cmExclusiveSpectra[i]; }
  delete [] particles->cmExclusiveSpectra;
  
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) { delete [] particles->energyCorrelations[i]; }
  delete [] particles->energyCorrelations;
  
}

/*******************************************************************************
 Finds the index corresponding to the value 'x0' in the array 'xarray'.
 ******************************************************************************/

// TO DO :: transform this template to accommodate double* array, not just constant double array !!!!!!!!!!!!!!!

template <int N> int FissionEvents::findEnergyIndex (double x0, double (&xarray)[N]) {
  int k0;
  int numberElements = N;
  for (int k=1; k<numberElements; k++) {
    if (xarray[k]>x0) {
      k0=k-1;
      break;
    }
  }
  if (k0<0 || k0>numberElements) { cerr << "STOP: problem in findEnergyIndex!\n"; exit(-1); }
  return k0;
}



/*******************************************************************************
 Initialize MPI communications.
 ******************************************************************************/
#ifdef MPIRUN
void FissionEvents::setMPI (void) {
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  return ;
}
#endif


/*!
 Save <nu>(A) results into file for R_T(A) analyses.
 */
void FissionEvents::printNubarA (string outputFilename) {
  ofstream out;
  out.open(&outputFilename[0]);
  for (int i=Asym; i<180; i++) {
    out << i << " " << neutrons.nuA[i] << " " << neutrons.nuA[Ac-i] << "\n";
  }
  out.close();
}

/*******************************************************************************
 Print formatted results.
 It should only be called once (for ip=0) in case of MPIRUN.
 ******************************************************************************/


void FissionEvents::printSummaryResults (string outputFilename, int ZAIDt, double Einc, int nevents, double alphaI) {
  
  emittedParticleType *particles;
  
  std::ostream* fp = &cout;
  std::ofstream fout;
  
  if (outputFilename!="") {
    fout.open(&(WORKDIR+outputFilename)[0]);
    fp = &fout;
  }
  
  *fp << setprecision(4) << setiosflags(ios::scientific);
  
  *fp << endl;
  *fp << "-------------------------------------------------------" << endl;
  *fp << "                      User Input                       " << endl;
  *fp << "-------------------------------------------------------" << endl;
  *fp << " ZAID Target                       = " << ZAIDt << endl;
  *fp << " Incident Neutron Energy (MeV)     = " << Einc << endl;
  *fp << " Spin distribution factor (alphaI) = " << alphaI << endl;
  *fp << " Time Gate (sec)                   = " << timeGate << endl;
  *fp << " Number of Monte Carlo events      = " << numberFissionEvents << endl;
  *fp << "-------------------------------------------------------" << endl << endl;
  
  
  *fp << "*******************************************************" << endl;
  *fp << "                C G M F  Summary of Results            " << endl;
  *fp << "*******************************************************" << endl;
  if (incidentNeutronEnergy!=0.0) {
    *fp << "                n ( " << incidentNeutronEnergy << " MeV ) + " << Zc*1000+Ac-1 << endl;
  } else {
    *fp << "               Spontaneous fission of " << Zc*1000+Ac << endl;
  }
  *fp << "*******************************************************" << endl << endl;
  
  *fp << fixed << setprecision(3);
  *fp << endl;
  
  *fp << std::string(96,'-') << endl;
  *fp << "|              | All Fragments | Light Fragments | Heavy Fragments | Pre-Fission |    Total    |\n";
  *fp << std::string(96,'-') << endl;
  
  *fp << "|      A       |" << setw(10) << aveMass   << "     | " << setw(10) << aveMassLF   << "      | " << setw(10) << aveMassHF   << "      |\n";
  *fp << "|      Z       |" << setw(10) << aveCharge << "     | " << setw(10) << aveChargeLF << "      | " << setw(10) << aveChargeHF << "      |\n";
  *fp << "| TKE/KE (MeV) |" << setw(10) << aveTKE    << "     | " << setw(10) << aveKEl      << "      | " << setw(10) << aveKEh      << "      |\n";
  *fp << std::string(96,'-') << endl;
  
  string title;
  for (int k=0; k<=1; k++) {
    
    if (k==0) {
      particles=&neutrons;
      title="Neutrons";
    } else if (k==1) {
      particles=&gammas;
      title="Gamma Rays";
    } else {
      particles=&ic;
      title="Internal Conversion";
    }
    
    *fp << endl << " << " << title << " >> " << endl << endl;
    *fp << std::string(96,'-') << endl;
    *fp << "|    < nu >    |" << setw(10) << particles->nubar << "     | " << setw(10) << particles->lfNubar << "      | " << setw(10) << particles->hfNubar << "      |";
    if (k==0 && numberPFN>0) {
      *fp << " " << setw(10) << float(numberPFN)/float(numberFissionEvents) << "  | " << setw(10) << particles->nubar+float(numberPFN)/float(numberFissionEvents) << "  |\n";
    } else {
      *fp << endl;
    }
    *fp << "|    < Ecm >   |" << setw(10) << particles->Ecm   << "     | " << setw(10) << particles->lfEcm   << "      | " << setw(10) << particles->hfEcm   << "      |\n";
    *fp << "|    < Elab >  |" << setw(10) << particles->Elab  << "     | " << setw(10) << particles->lfElab  << "      | " << setw(10) << particles->hfElab  << "      |";
    if (k==0 && numberPFN>0) {
      double aveEnergyTotal = ( particles->nubar*particles->Elab+float(aveEnergyPFN)/float(numberFissionEvents) ) / ( particles->nubar+float(numberPFN)/float(numberFissionEvents));
      *fp << " " << setw(10) << float(aveEnergyPFN)/float(numberPFN) << "  | " << setw(10) << aveEnergyTotal << "  |\n";
    } else {
      *fp << endl;
    }
    *fp << std::string(96,'-') << endl;
    
  }
  
  if (fout.is_open()) { fout.close(); }
  
}




void FissionEvents::printCompactSummaryResults () {
  
  cout << neutrons.nubar << " " << neutrons.Ecm << "\n";
  
}



/*******************************************************************************
 Stores all CGMF results into data files readable by GNUPLOT. It calls the
 method saveResultsToDataFile (emittedParticleType, string) for each neutron
 and gamma-ray results.
 ******************************************************************************/

void FissionEvents::saveResultsToGnuplot () {
  
  saveGeneralResults ("results.CGMF"+fileExt);
  
  //  saveMonteCarloYieldsToFile ("yields.CGMF"); // redundant -- to remove
  
  saveParticleResultsToDataFile (neutrons, "neutrons.CGMF"+fileExt);
  saveParticleResultsToDataFile (gammas, "gammas.CGMF"+fileExt);
  //	saveParticleResultsToDataFile (ic, "ic.CGMF");
  
}


/******************************************************************************
 Save results common to both neutrons and photons into a file.
 ******************************************************************************/
void FissionEvents::saveGeneralResults (string filename) {
  
  ofstream out;
  out.open(&(WORKDIR+filename)[0]);
  out << setprecision(4) << setiosflags(ios::scientific);
  out << "#\n";
  out << "# CGMF General Results\n";
  out << "#\n\n";
  
  int gnuplotIndex=-1;
  
  //-- Mass-dependent quantities: Y(A), <Ui>(A), <Ji>(A), YApost(A)
  out << "\n# [gnuplot #" << ++gnuplotIndex << "]\n";
  out << "# A   Ypre(A)     <Ui>(A)    <Ji>(A)   Ypost(A)\n\n";
  for (int i=0; i<NUMA; i++)
    out << setw(3) << i << " " << YA[i] << " " << " " << UiA[i] << " " << initialSpinA[i] << " " << YApost[i] << "\n";
  
  //-- Charge-dependent quantities: Y(Z), ...
  out << "\n#  [gnuplot #" << ++gnuplotIndex << "]\n";
  out << "#  Z   Y(Z) \n\n";
  for (int i=0; i<NUMZ; i++)
    out << setw(2) << i << " " << YZ[i] << "\n";
  
  //-- Total Kinetic Energy Yields Y(TKE)
  out << "\n# [gnuplot #" << ++gnuplotIndex << "]\n";
  out << "#    TKE     Y(TKE)   <Ji>(TKE)\n";
  out << "#   (MeV)    (n/f)     (hbar)\n\n";
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++)
    out << excitationEnergyGrid[i] << " " << YTKE[i] << " " << initialSpinTKE[i] << "\n";
  
  //-- Initial Excitation Energy Distributions Y(Ul), Y(Uh) and Y(Utot)
  out << "\n# [gnuplot #" << ++gnuplotIndex << "]\n";
  out << "# U (MeV)     Y(Ul)     Y(Uh)     Y(Utot)\n\n";
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++)
    out << excitationEnergyGrid[i] << " " << YUl[i] << " " << YUh[i] << " " << YUtot[i] << "\n";
  
  //-- Spin probabilities in light and heavy fragments
  out << "\n#  [gnuplot #" << ++gnuplotIndex << "]\n";
  out << "#  Spin   LF Proba    HF Proba\n";
  out << "#  (hbar)\n\n";
  for (int i=0; i<50; i++)
    out << setw(6) << i << "   " << lfPspin[2*i] << "   " << hfPspin[2*i] << "\n";
  
  //-- Total gamma-ray energy and multiplicity as a function of neutron multiplicity
  out << "\n# [gnuplot #" << ++gnuplotIndex << "]\n";
  out << "# nu   <Eg^tot>    <Ng>\n\n";
  for (int i=0; i<20; i++) {
    if (countsNu[i]!=0) {
      out << setw(4) << i << " " << totalGammaEnergyNu[i] << " " << gammaMultiplicityNu[i] << "\n";
    }
  }
  
  //-- Total gamma-ray energy as a function of light fragment mass number
  out << "\n# [gnuplot #" << ++gnuplotIndex << "]\n";
  out << "#  A_LF   <Eg^tot>\n\n";
  for (int i=0; i<=Asym; i++) {
    if (YA[i]!=0) {
      out << setw(3) << i << " " << totalGammaEnergyA[i] << "\n";
    }
  }
  
  //-- Total gamma-ray energy as a function of TKE
  out << "\n# [gnuplot #" << ++gnuplotIndex << "]\n";
  out << "#   TKE   <Eg^tot>\n";
  out << "#  (MeV)     (MeV)\n\n";
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    if (YTKE[i]!=0) {
      out << excitationEnergyGrid[i] << " " << totalGammaEnergyTKE[i] << "\n";
    }
  }
  
  out << "\n# [gnuplot #" << ++gnuplotIndex << "]\n";
  out << "# Angular distributions n-LF\n";
  out << "# Theta   Y(Theta_cm)   Y(Theta_lab)\n";
  out << "#  (deg)\n\n";
  
  for (int i=0; i<NUMANGLES; i++) {
    out << i*dTheta << " " << YThetacm[i] << " " << YThetalab[i] << " " << "\n";
  }
  
  out.close();
  
}


/******************************************************************************
 Save CGMF results into a data file (readable by gnuplot) for a specific type
 of emitted particle (neutron or gamma). It also creates an additional
 output file to store particle-particle energy correlations.
 ******************************************************************************/
void FissionEvents::saveParticleResultsToDataFile (emittedParticleType particles, string filename) {
  
  ofstream out;
  out.open(&(WORKDIR+filename)[0]);
  out << setprecision(4) << setiosflags(ios::scientific);
  out << "#\n# CGMF Results\n#\n";
  
  int gnuplotIndex=-1;
  
  //-- Laboratory PFNS
  out << "# [gnuplot #" << ++gnuplotIndex << "] lab spec, cm spec, LF lab spec, LF cm spec, HF lab spec, HF cm spec\n\n";
  for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++)
    out << spectrumEnergyGrid[i] << " " << particles.labSpectrum[i] << " " << particles.cmSpectrum[i] << " " <<
    particles.labLightFragmentSpectrum[i] << " " << particles.cmLightFragmentSpectrum[i] << " " <<
    particles.labHeavyFragmentSpectrum[i] << " " << particles.cmHeavyFragmentSpectrum[i] << "\n";
  
  //-- cm PFNS and multiplicity-dependent exclusive spectra
  out << "\n# [gnuplot #" << ++gnuplotIndex << "] c.m. spectrum\n";
  out << "# Energy    Total       nu=1      nu=2      nu=3       nu=4      nu=5\n";
  out << "# (MeV)     (1/MeV)\n\n";
  for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++)
    out << spectrumEnergyGrid[i] << " " << particles.cmSpectrum[i] << " " <<
    particles.cmExclusiveSpectra[1][i] << " " << particles.cmExclusiveSpectra[2][i] << " " <<
    particles.cmExclusiveSpectra[3][i] << " " << particles.cmExclusiveSpectra[4][i] << " " <<
    particles.cmExclusiveSpectra[5][i] << "\n";
  
  //-- Multiplicity distributions P(nu), and same for light and heavy fragments resp.
  out << "\n# [gnuplot #" << ++gnuplotIndex << "] P(nu) P_LF(nu) P_HF(nu)\n\n";
  for (int i=0; i<NUMMULT; i++)
    out << i << " " << particles.Pnu[i] << " " << particles.lfPnu[i] << " " << particles.hfPnu[i] << "\n";
  
  //-- Fragment mass-dependent results: mass Yields A Y(A), <nu>(A), <Ecm>(A), <Ecm>(A) w/ cut, etc.
  out << "\n# [gnuplot #" << ++gnuplotIndex << "] A  Y(A)   <nu>(A)   <Ecm>(A)  <Ecm>(A)|cut  <Elab>(A)\n\n";
  for (int i=0; i<NUMA; i++)
    out << i << " " << YA[i] << " " << particles.nuA[i] << " " << particles.EcmA[i] <<
    " " << particles.EcmAcut[i] << " "  << particles.ElabA[i] << "\n";
  
  //-- Z-dependent results: Y(Z), ...
  /*out << "\n#  [gnuplot #" << ++gnuplotIndex << "] Z   Y(Z) \n\n";
   for (int i=0; i<NUMZ; i++) {
   out << i << " " << YZ[i] << "\n";
   }*/
  
  
  //-- Angular distributions
  out << "\n# [gnuplot #" << ++gnuplotIndex << "] Theta   n-LF    n-n\n";
  out << "#       (deg)   (arb. u.)    (arb. u.)\n\n";
  if (particles.particleType=='n') {
    for (int i=0; i<NUMANGLES; i++) {
      out  << i*dTheta << " " << particles.nLF_angles[i] << " " << particles.nn_angles[i] << " " << sin(i*dTheta*PI/180.0) << "\n";
    }
  } else {
    out << "0 0 0\n";
  }
  
  //-- Y(TKE), <nu>(TKE), <Ecm>(TKE), <Elab>(TKE)
  out << "\n# [gnuplot #" << ++gnuplotIndex << "] TKE   Y(TKE)   <nu>  <Ecm>  <Elab>\n";
  out << "#             (MeV)   (n/f)  (MeV)   (MeV)   (MeV) \n\n";
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    out << excitationEnergyGrid[i] << " " << YTKE[i] << " " << particles.nuTKE[i] <<
    " " << particles.EcmTKE[i] << " " << particles.ElabTKE[i] << "\n";
  }
  
  //-- Einc, <nu>, <nu_l>, <nu_h>, <Ecm>, <Ecm_l>, <Ecm_h>, <Elab>, <Elab_l>, <Elab_h>
  out << "\n# [gnuplot #" << ++gnuplotIndex << "] Einc, <nu>, <nu_l>, <nu_h>, <Ecm>, <Ecm_l>, <Ecm_h>, <Elab>, <Elab_l>, <Elab_h>\n\n";
  out << incidentNeutronEnergy << " " << particles.nubar << " " << particles.lfNubar
  << " " << particles.hfNubar << " " << particles.Ecm << " " << particles.lfEcm
  << " " << particles.hfEcm << " " << particles.Elab << " " << particles.lfElab
  << " " << particles.hfElab << "\n";
  
  
  //-- Laboratory PFNS *** for gammas only ***
  //-- continuum-to-continuum, continuum-to-discrete, discrete-to-discrete
  if (particles.particleType=='g') {
    out << "\n# [gnuplot #" << ++gnuplotIndex << "] lab spec. c-c, c-d, d-d\n\n";
    for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++)
      out << spectrumEnergyGrid[i] << " " << particles.ccLabSpectrum[i] << " " <<
      particles.cdLabSpectrum[i] << " " << particles.ddLabSpectrum[i] << "\n";
  }
  
  out.close();
  
  /*
   //-- continuum-to-continuum, continuum-to-discrete, discrete-to-discrete
   //-- per fission fragment mass!
   if (particles.particleType=='g') {
   gnuplotIndex=-1;
   string fn;
   fn = WORKDIR + "exclusiveGammaSpectra.CGMF";
   out.open(&fn[0]);
   for (int iA=0; iA<NUMA; iA++) {
   if (YA[iA]!=0.0) {
   out << "\n# [gnuplot #" << ++gnuplotIndex << "] A, lab spec. c-c, c-d, d-d\n";
   out << "# " << iA << "\n\n";
   for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++)
   out << spectrumEnergyGrid[i] << " " << particles.ccLabSpectrumA[iA][i] << " " <<
   particles.cdLabSpectrumA[iA][i] << " " << particles.ddLabSpectrumA[iA][i] << "\n";
   }
   }
   out.close();
   }
   */
  
  /*
   // particle-particle Energy Correlations -------------------------------------
   //  if (ip==0) {
   string correlationFilename;
   correlationFilename = WORKDIR + filename + ".corr";
   out.open(&correlationFilename[0]);
   out << "#\n# particle-particle energy correlations\n#\n";
   for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) {
   for (int j=0; j<NUMBER_OF_PRIVATE_GRID; j++) {
   out << .5 * ( custom_energy_grid[i] + custom_energy_grid[i+1] ) << " "
   << .5 * ( custom_energy_grid[j] + custom_energy_grid[j+1] ) << " "
   << particles.energyCorrelations[i][j] << " " << particles.energyCorrelations[i][j] << "\n";
   }
   out << endl ;
   }
   out << endl ;
   out.close();
   //  }
   
   if (particles.particleType=='n') {
   string correlationFilename;
   correlationFilename = WORKDIR + filename + ".ad.corr";
   out.open(&correlationFilename[0]);
   out << "#\n# particle-particle angular correlations\n#\n";
   for (int i=0; i<NUMANGLES; i++) {
   for (int j=0; j<NUMANGLES; j++) {
   out << 0.5*i << " " << 0.5*j << " " << particles.angularCorrelations[i][j] << "\n";
   }
   out << endl ;
   }
   out << endl ;
   out.close();
   }
   */
  
  // end particle-particle correlations ----------------------------------------
  
  /*
   // save c.m. spectrum as a function of fragment mass -------------------------
   string spectrumFilename;
   spectrumFilename = WORKDIR + filename + ".cmSpectrumA";
   
   //	cout << spectrumFilename << endl;
   
   out.open(&spectrumFilename[0]);
   out << setprecision(4) << setiosflags(ios::scientific);
   
   out << "#\n# c.m. spectrum as a function of primary fission fragment mass\n#";
   
   
   cout << "HERE?\n";
   cout << particles.counts << " " << particles.nubar << endl;
   
   double sum=0.0;
   for (int i=0; i<NUMA; i++) {
   for (int j=0; j<NUMBER_SPECTRUM_ENERGY_GRID; j++) {
			sum+=(particles.cmSpectrumA[i][j]*spectrumEnergySteps[j]);
   }
   }
   cout << sum << endl;
   
   sum=0.0;
   for (int i=0; i<NUMA; i++) {
   for (int j=0; j<301; j++) {
			sum+=particles.cmHistogramA[i][j];
   }
   }
   cout << sum << endl;
   
   gnuplotIndex=-1;
   for (int i=0; i<NUMA; i++) {
   if (YA[i]!=0) {
   out << "\n# [gnuplot #" << ++gnuplotIndex << "] A = " << i << "\n\n";
   for (int j=0; j<NUMBER_SPECTRUM_ENERGY_GRID; j++) {
   out << spectrumEnergyGrid[j] << " " << particles.cmSpectrumA[i][j] << " "
   << particles.cmSpectrumNuA[1][i][j] << " " << particles.cmSpectrumNuA[2][i][j] << " "
   << particles.cmSpectrumNuA[3][i][j] << " " << particles.cmSpectrumNuA[4][i][j] << " "
   << YA[i]*particles.counts*200 << "\n";
   }
   }
   }
   out.close();
   // end c.m. spectrum as a function of fragment mass --------------------------
  	
   spectrumFilename = WORKDIR + filename + ".cmHistogramA";
   out.open(&spectrumFilename[0]);
   out << setprecision(4) << setiosflags(ios::scientific);
   
   out << "#\n# c.m. histograms as a function of primary fission fragment mass\n#";
   
   gnuplotIndex=-1;
   for (int i=0; i<NUMA; i++) {
   if (YA[i]!=0) {
   out << "\n# [gnuplot #" << ++gnuplotIndex << "] A = " << i << "\n\n";
   for (int j=0; j<301; j++) {
   out << j*0.05 << " " << particles.cmHistogramA[i][j]/0.05 << " "
   << particles.cmHistogramNuA[1][i][j]/0.05 << " " << particles.cmHistogramNuA[2][i][j]/0.05 << " "
   << particles.cmHistogramNuA[3][i][j]/0.05 << " " << particles.cmHistogramNuA[4][i][j]/0.05 << " "
   << YA[i]*particles.counts*200 << "\n";
			}
   }
   }
   out.close();
   */
  
}

/*!
 
 Save primary fission fragment yields to an output file (readable by gnuplot).
 Y(A), Y(Z), Y(TKE), ...
 
 */
void FissionEvents::saveMonteCarloYieldsToFile (string filename) {
  
  int gnuplotIndex = 0;
  ofstream out;
  out.open(&(WORKDIR+filename)[0]);
  out << setprecision(4) << setiosflags(ios::scientific);
  out << "#\n# CGMF Monte Carlo Generated Fragment Yields\n#\n";
  
  // Y(A)
  out << "#\n# [gnuplot #" << gnuplotIndex++ << "] Mass Yields Y(A)\n\n";
  for (int i=0; i<NUMA; i++) out << i << " " << YA[i] << " " << YApost[i] << endl;
  
  // Y(Z)
  out << "\n# [gnuplot #" << gnuplotIndex++ << "] Charge Yields Y(Z)\n\n";
  for (int i=0; i<NUMZ; i++) out << i << " " << YZ[i] << endl;
  
  // Y(KE)
  out << "\n# [gnuplot #" << gnuplotIndex++ << "] Fragment Kinetic Energy Y(KE)\n\n";
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++)
    out << excitationEnergyGrid[i] << " " << YKEl[i] << " " << YKEh[i] << endl;
  
  // Y(TKE)
  out << "\n# [gnuplot #" << gnuplotIndex++ << "] Total Kinetic Energy Y(TKE)\n\n";
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++)
    out << excitationEnergyGrid[i] << " " << YTKE[i] << endl;
  
  
}

/*!
 
 Set the (Zc,Ac) of the compound (fissioning) nucleus, as well as the symmetric
 fission fragment (Zsym,Asym).
 
 */
void FissionEvents::setZAIDcn (int ZAIDt, double Einc) {
  
  int ZAIDcn;
  
  if (Einc!=0.0) { // neutron-induced fission
    ZAIDcn = ZAIDt+1;
    incidentNeutronEnergy = Einc;
  } else { // spontaneous fission
    ZAIDcn = ZAIDt;
    incidentNeutronEnergy=0.0;
  }
  
  Ac    = ZAIDcn%1000;
  Zc    = int(ZAIDcn/1000);
  
  // symmetric fission fragment
  Asym = Ac/2;
  Zsym = Zc/2;
  
}

/*!
 Save results of gamma rays for GEANT simulations (J.Ullmann, LANSCE-NS & M.Jandel, CNR).
 */
#if defined(GEANT) && defined(MPIRUN)
void FissionEvents::saveResultsGEANT (void) {
  
  ofstream OUT;
  int i;
  
  if (ip==0) remove ("gammaGEANT.cgmf");
  
  for (int n=0; n<np; n++) {
    
    if (n==ip) {
      
      OUT.open ("gammaGEANT.cgmf", ios::app);
      OUT << setprecision(4) << setiosflags(ios::fixed);
      
      i=0;
      do {
        for (int j=0; j<20; j++) if (gammaGEANT[i][j]!=0.0) OUT << gammaGEANT[i][j] << " "; // light fragment
        i++;
        for (int j=0; j<20; j++) if (gammaGEANT[i][j]!=0.0) OUT << gammaGEANT[i][j] << " "; // heavy fragment
        i++;
        OUT << "\n";
      } while (i<indexGEANT);
      
      
      OUT.close ();
    }
    
    MPI_Barrier (comm);
    
  }
}
#endif





