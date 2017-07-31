// $Id: constant.h 5 2011-08-09 02:08:58Z kawano $
/****************************/
/*      ARRAY SIZE          */
/****************************/

const int MAX_ENERGY_BIN   = 2000 ;  /* maximum energy bins                   */
const int MAX_J            =   70 ;  /* maximum J-value                       */
const int MAX_LEVELS       =  200 ;  /* maximum discrete levels               */
const int MAX_GAMMA_BRANCH =  100 ;  /* maximum gamma-ray branches            */
const int MAX_MULTIPOL     =    3 ;  /* multipolarity, E1, E1(def), M1, E2    */
const int MAX_GDR          =    6 ;  /* maximum GDRs and pygmy resonance      */
const int MAX_CHANNEL      =    2 ;  /* maximum decay channel                 */
const int MAX_GAMMA_LINES  = 1000 ;  /* maximum line gamma-rays               */
const int MAX_COMPOUND     =    5 ;  /* maximum number of compound nucleus    */
const int SPECTRA_OUTPUT   =    4 ;  /* number of spectra to be calculated    */


/****************************/
/*      GENERAL             */
/****************************/

const double EXP1       =  2.71828182845904523536  ;  /* Napier's constant    */
const double PI         =  3.14159265358979323846  ;  /* circular constant    */
const double PI2        =  6.28318530717958647692  ;  /* PI*2                 */
const double PI4        = 12.56637061435917295384  ;  /* PI*4                 */
const double PI_2       =  1.57079632679489661923  ;  /* PI/2                 */
const double PI_4       =  0.785398163397448309616 ;  /* PI/4                 */
const double PI4_3      =  4.18879020478639053     ;  /* PI*4/3               */
const double SQRTPI     =  1.77245385090552        ;  /* sqrt(PI)             */
const double SQRT_1_PI  =  0.564189583547756286948 ;  /* sqrt(1/PI)           */
const double SQRT_2_PI  =  0.797884560802865406    ;  /* sqrt(2/PI)           */
const double LOG_2PI    =  1.83787706640934548     ;  /* log(2PI)             */
const double LOG_2      =  0.6931471805599453      ;  /* log(2.0)             */


/****************************/
/*      DEFAULT PARAMETER   */
/****************************/

const double NORM_FACT      = 10.0 ; /* conversion factor to [mb]             */
const double SPIN_CUTOFF    = 1e-8 ; /* cut-off of total anuglar momentum     */


/****************************/
/*      PHYSICAL CONSTANT   */
/****************************/

const double HBAR       =  6.58212196e-22 ; /* Planck's constant/2pi [MeV sec]*/
const double HBARSQ     =  4.33243296e-43 ; /* HBAR*HBAR                      */
const double COULOMB    =  1.60217733e-19 ; /* J = 1eV                        */
const double COULOMBSQ  =  2.56697220e-38 ; /* COULOMB*COULOMB                */
const double PERMITTIV  =  5.60958617e+37 ; /* permittivity [MeV fm /C^2]     */

const double AMUNIT     =  931.4943335    ; /* MeV = 1amu                     */
const double VLIGHT     =  2.99792458e+23 ; /* light velocty [fm/sec]         */
const double VLIGHTSQ   =  8.98755179e+46 ; /* VLIGHT*VLIGHT                  */

const double MNEUTRON   =  1.008664891    ; /* Neutron Weight  [amu]          */
const double MPROTON    =  1.007276487    ; /* Proton Weight   [amu]          */
const double MDEUTERON  =  2.0141018      ; /* Deuteron Weight [amu]          */
const double MTRITON    =  3.0160494      ; /* Triton Weight   [amu]          */
const double MHELIUM3   =  3.0260294      ; /* Helium-3 Weight [amu]          */
const double MALPHA     =  4.0026032      ; /* Alpha Weight    [amu]          */

const double ENEUTRON   =   8.071431      ; /* Neutron Mass Excess  [MeV]     */
const double EPROTON    =   7.289034      ; /* Proton Mass Excess   [MeV]     */
const double EDEUTERON  =  13.13584       ; /* Deuteron Mass Excess [MeV]     */
const double ETRITON    =  14.94994       ; /* Triron Mass Excess   [MeV]     */
const double EHELIUM3   =  14.93132       ; /* Helium-3 Mass Excess [MeV]     */
const double EALPHA     =   2.42492       ; /* Alpha Mass Excess    [MeV]     */
const double EELECTRON  =  0.510998910    ; /* Electron Mass        [MeV]     */

const double CSPO       =   2.04553       ; /* (HBAR/Mpi_meson/VLIGHT)^2      */
const double ALPHA      = 7.2973525376e-03; /* fine structure constant        */
