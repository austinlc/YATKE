// $Id: config-BETA.h 5 2011-08-09 02:08:58Z kawano $
/* 
   config.h : 
        customize CGM, define data directory, system dependent parameters
*/


/*----------------------------------------------------------
Top directory in which CGM can find data files - 
     RIPL3 level data   /dir/to/levels/z001.dat
     kcksyst.dat        /dir/to/kcksyst.dat
     beta strength      /dir/to/bstrength/strint053137.dat 
*/

#define  DATADIR   "/usr/local/share/coh"


/*----------------------------------------------------------
Neutron emission channel
     This value allows / suppresses competing neutron emission
*/

const bool INCLUDE_NEUTRON_EMISSION = false;

// if excitation energy is Bn + this, gamma-ray emission ignored
// undef to include all gammas

const double IGNORE_CAPTURE_GAMMA_RAY = 0.0;


/*----------------------------------------------------------
Energy bin width in the continuum
*/

//const double ENERGY_BIN = 0.10;  // 100 keV
const double ENERGY_BIN = 0.01;  // 10 keV
//const double ENERGY_BIN = 0.05;  // 50 keV


/*----------------------------------------------------------
Beta strength function Gaussian broadening resolution
*/

const double PROFILE_BROADENING_WIDTH = 0.1;  // 100 keV


/*----------------------------------------------------------
Binary reaction calculation
     suppress multiple step reactions,
     calculate decay from the initial state only
*/

const bool BINARY_REACTION_ONLY = false;


/*----------------------------------------------------------
Custom energy grid for spectra
*/

#undef HAVE_PRIVATE_ENERGY_GRID

#ifdef   HAVE_PRIVATE_ENERGY_GRID
const int NUMBER_OF_PRIVATE_GRID = 565;
#define   GRID_STRUCTURE_FILE    "privategrid1.h"
//const int NUMBER_OF_PRIVATE_GRID = 272;
//#define   GRID_STRUCTURE_FILE    "privategrid2.h"
#endif


/*----------------------------------------------------------
Monte Carlo control
   RANDOM_SEED_BY_TIME
       change random number seed each time by using
       internal clock

   PERTURB_EXCITATION ENERGY
       random noize is added to the energies for
       the continuum to continuum transitions

   EVENT_OUTPUT_FORMAT
       0: output MC history in the CGM standard format
       1: output MC history in the GEANT4 input format
       2: output MC history in the FFD input format
*/

const bool RANDOM_SEED_BY_TIME = false;
const bool PERTURB_EXCITATON_ENERGY = true;
const int  EVENT_OUTPUT_FORMAT = 0;

