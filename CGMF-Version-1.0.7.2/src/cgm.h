// $Id: cgm.h 5 2011-08-09 02:08:58Z kawano $
/*
   cgm.h :
        prototype of functions for the CGM code
 */

#ifndef __STRUCTUR_H__
#define __STRUCTUR_H__
#include "structur.h"
#endif

#ifndef CGM_TOPLEVEL
//extern  Calculation  ctl;
extern  Nucleus      *ncl;
extern  double       *spc[SPECTRA_OUTPUT];
#endif


// to pass as argument to MCNP fission event
struct fissionEventType {
	int nu;  // number of neutrons
	int nug; // number of photons
	double neutronEnergies [10];
	double photonEnergies [30];
};


/*** spin int / half-int distinction */
#define halfint(x)       ( (((int)(x*2.0+0.001) % 2)==1) ? 0.5 : 0.0 )

/*** parity calculator */
#define parity(a)        ((((a)/2)%2==1) ? -1 : 1)


/**************************************/
/*      Gamma-ray multipolality       */
/**************************************/
typedef enum {SL, GL}           GammaProfile;
typedef enum {E1=0, M1=1, E2=2} GammaMultipol;

#define E1_MULTIPOL         1
#define M1_MULTIPOL         1
#define E2_MULTIPOL         2



/**************************************/
/*      omcalc.cpp                    */
/**************************************/
int     omCalc                          (double, Pdata *, ZAnumber *, double,
                                         double *, CrossSection *);

/**************************************/
/*      setup.cpp                     */
/**************************************/
void    statSetupInitSystem             (int, Pdata *);
int     statSetupEnergyBin              (Nucleus *);
void    statSetupLevelDensityParameter  (Nucleus *, LevelDensity *);
void    statSetupLevelDensity           (Nucleus *, LevelDensity *);
void    statGenerateLevels              (int, Nucleus *, LevelDensity *);
void    statFixDiscreteLevels           (Nucleus *);
void    statSetupGdrParameter           (Nucleus *, double, GDR *, double);
int     statTotalNeutrons               (double, ZAnumber *);


/**************************************/
/*      levden.cpp                    */
/**************************************/
double  ldLevelDensity                  (double, double, LevelDensity *);
double  ldDensityParameter              (double, double, LevelDensity *);
double  ldSpinDistribution              (double, double, double);
double  ldParityDistribution            (int);
double  ldLevelSpinCutoff               (int, Level *);
void    ldTGconnect                     (double, double, double, LevelDensity *);
void    ldTextrapolate                  (double, double, LevelDensity *);
void    ldGextrapolate                  (double, LevelDensity *);


/**************************************/
/*      popinit.cpp                   */
/**************************************/
void    statSetupInitialPopulationSpec  (double, int, int, double, Transmission *);
void    statSetupInitialPopulationBeta  (double, int, int);
int     statStoreLevelExcite            (Nucleus *);
void    statClearPopulation             (Nucleus *);


/**************************************/
/*      stattrans.cpp                 */
/**************************************/
void    statStoreContinuumTransmission  (int, int, Pdata *, Transmission *);
void    statStoreDiscreteTransmission   (int, int, Pdata *, Transmission *);
void    statStoreGammaTransmission      (int, double **, Nucleus *);


/**************************************/
/*      gtrans.cpp                    */
/**************************************/
void    gdrM1norm                       (double, double, GDR *);
void    gdrParameterSave                (int, GDR *);
double  gdrGammaTransmission            (GammaProfile,int,double, GammaMultipol, double, double);
void    gdrRenormE1StrengthFunction     (int, int, int, double, Nucleus *, double **);


/**************************************/
/*      specmain.cpp                  */
/**************************************/
void    specMain                        (double, int, double, double **);
void    spectra                         (int, Pdata *, double **);
void    spectraBinary                   (int, Pdata *, double **);
int     specFindEnergyBin               (double, double);
int     specFindCustomGrid              (int, double, double *);


/**************************************/
/*      specmc.cpp                    */
/**************************************/
void    specMCMain                      (double, int, double, double, unsigned long, double **);
void    specMCMain2                     (int, int, double, double, int, double, double, double **, fissionEventType *);

/**************************************/
/*      specevap.cpp                  */
/**************************************/
void    specEvaporation                 (int, Pdata *, double **);
void    specEvaporationBin              (int, Pdata *, double **);

/**************************************/
/*      betamain.cpp                  */
/**************************************/
void    betaMain                        (int, double **);


/**************************************/
/*      betaprofile.cpp               */
/**************************************/
void    betaStrengthProfile             (ZAnumber *, double *, double *, Nucleus *);
void    betaInitialPopulation           (double *, double *, Nucleus *);


/**************************************/
/*      elepop.cpp                    */
/**************************************/
void    betaElectronSpectrum            (double, double *, double *, Nucleus *, double **);


/**************************************/
/*      gampop.cpp                    */
/**************************************/
double  specTransitionGamma             (Statcalcmode, int, int, int,
                                         double **, double, Nucleus *, double *);


/**************************************/
/*      parpop.cpp                    */
/**************************************/
double  specTransitionParticle          (Statcalcmode, int, int, int,
                                         Transmission *, Transmission *,
                                         double, Nucleus *, Nucleus *, double *);

/**************************************/
/*      gamcas.cpp                    */
/**************************************/
int     specGammaCascade                (int, double *, Nucleus *);
int     specGammaDiscrete               (int, double *, Nucleus *);


/**************************************/
/*      output.cpp                    */
/**************************************/
void    cgmPrintSpectra                 (bool, double, double **);
void    cgmPrintBroaden                 (double, double **);
void    cgmPrintCustomGrid              (int, double *, double *, double **);
void    cgmPrintPopulation              (Nucleus *, Nucleus *);
int     cgmZeroCut                      (int, double **);


/**************************************/
/*      cmslab.cpp                    */
/**************************************/
void    cgmLabSpectrum                  (double, double *, double **);
