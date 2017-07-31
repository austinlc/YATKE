// $Id: config-ORIG.h 5 2011-08-09 02:08:58Z kawano $
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

//#define DATADIR   "/usr/local/share/coh/"
#define DATADIR "/home/austinlc/CGMF/CGMF-Version-1.0.7.2/data/"

#define CGMF
//#undef CGMF


/*----------------------------------------------------------
 Current version of CGMF code.
 */
const std::string CGMF_VERSION = "1.0.7";


/*----------------------------------------------------------
 Neutron emission channel
 This value allows / suppresses competing neutron emission
 */

const bool INCLUDE_NEUTRON_EMISSION = true;

// if excitation energy is Bn + this, gamma-ray emission ignored
// set zero to include all gammas

const double IGNORE_CAPTURE_GAMMA_RAY = 0.0;


/*----------------------------------------------------------
 Energy bin width in the continuum
 */

//const double ENERGY_BIN = 0.10;  // 100 keV
//const double ENERGY_BIN = 0.01;  // 10 keV
const double ENERGY_BIN = 0.05;  // 50 keV

/*----------------------------------------------------------
 Experimental time coincidence window (in seconds)
 By default, set to -99.0 to let all levels decay to the 
 ground.state.
 */

// const double EXPERIMENTAL_TIME_WINDOW = -99.0; // [default]
const double EXPERIMENTAL_TIME_WINDOW = 1e-8; // [default] but can replaced by timeGate from user input with option '-t'


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
 Neutron spectrum calculated by Evaporation model
 no gamma-ray emission
 */
const bool EVAPORATION_SPECTRUM = false;


/*----------------------------------------------------------
 Custom energy grid for spectra
 */

#undef HAVE_PRIVATE_ENERGY_GRID
//#define HAVE_PRIVATE_ENERGY_GRID

#ifdef   HAVE_PRIVATE_ENERGY_GRID
//const int NUMBER_OF_PRIVATE_GRID = 565;
//#define   GRID_STRUCTURE_FILE    "privategrid1.h"
const int NUMBER_OF_PRIVATE_GRID = 272;
#define   GRID_STRUCTURE_FILE    "privategrid2.h"
#endif


// no gamma decay if excitation energy of continuum is less than this
const double CONTINUUM_LOWER_CUT = 0.02;


/*----------------------------------------------------------
 Discrete gamma-ray Internal Conversion control
 ICC considered in discrete transisions
 */

const bool INCLUDE_INTERNAL_CONVERSION = true;


/*----------------------------------------------------------
 Monte Carlo control
 RANDOM_SEED_BY_TIME
 change random number seed each time by using
 internal clock
 
 PERTURB_EXCITATION ENERGY
 random noise is added to the energies for
 the continuum to continuum transitions
 
 EVENT_OUTPUT_FORMAT
 0: output MC history in the CGM standard format
 1: output MC history in the GEANT4 input format
 2: output MC history in the FFD input format
 */

const bool RANDOM_SEED_BY_TIME = true;
const bool PERTURB_EXCITATON_ENERGY = true;
const int  EVENT_OUTPUT_FORMAT = 0;


