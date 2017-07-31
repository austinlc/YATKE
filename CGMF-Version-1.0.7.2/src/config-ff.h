/*!
 Additional configuration file for CGM-F
 */

/*----------------------------------------------------------
 Turn on/off MPI calls.
 */

#define MPIRUN // to turn on MPI calls
//#undef MPIRUN    // to turn off MPI calls

#define HISTORIES // to write out Monte Carlo history file
//#undef HISTORIES

/*----------------------------------------------------------
 Prompt gamma-ray results are saved into a format usable for
 GEANT simulations of the DANCE g-ray detector.
 */

//#define GEANT
#undef GEANT

const std::string WORKDIR = "./";
//const std::string WORKDIR = "/Users/talou/tmp/";

const std::string EXT = ".10ns"; // extension to result filenames

const int    NUMA   =  300; // number of masses A
const int    NUMZ   =  100; // number of charges Z
const int    NUMTKE =  300; // number of Total Kinetic Energy values

const int      NUME =  401; // number of energies in level density tables; dE=0.25 MeV; up to Emax=100 MeV
const double deltaE = 0.25; // energy-bin size used in level density tables
const int     NUME2 = 32;

const int    NUMdZ  =   21; // [-dZ:+dZ] if dZ=10 for charge distribution around most probable Zp[A] 

const int   NUMMULT =   50; // number of multiplicities

const int NUMANGLES =   73; // number of angles in angular distribution
const double dTheta =  2.5; // angular bins (degrees)

const int MAX_NUMBER_PARTICLES = 50; // max. number of particles (n or g) emitted per fragment in a fission event

const int NUMBER_SPECTRUM_ENERGY_GRID = 641; //551;

