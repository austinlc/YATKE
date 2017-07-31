/*
 *  CGM-F
 *  [ physics.h ]
 * 
 *  Version: 1.0.1
 *  April 18, 2012
 *
 *  P.Talou, talou@lanl.gov
 *
 */

const double pi                           = 3.14159265359;
const double twopi                        = 6.28318530718;   // 2*pi
const double elementaryCharge             = 1.602176462e-19; // [C]
const double inverseFineStructureConstant = 137.03599976;
const double amu                          = 931.494013e+6;   // atomic mass unit [eV]
const double amuMeV                       = 931.494013;      // [MeV]
const double PlanckConstant               = 4.13566727e-15;  // [eV.s]
const double PlanckConstantOverTwoPi      = 6.58211889e-16;  // Planck constant / (2 pi) [eV.s]
const double BoltzmannConstant            = 8.617342e-5;     // [eV/K]
const double speedOfLight                 = 299792458.0;     // [m/s]
const double AvogadroNumber               = 6.02214199e23;   // [/mol]

const double hbarc = PlanckConstantOverTwoPi*speedOfLight*1e9;  // Hbar.c [MeV.fm]

//--- Masses [amu] ---

const double neutronMass  = 1.00866491578;
const double electronMass = 5.48579911e-4;
const double protonMass   = 1.00727646688;	
const double deuteronMass = 2.01355321271;
const double tritonMass   = 3.016049268;
const double he3Mass      = 3.01493223469;
const double alphaMass    = 4.0015061747;

//--- Energies needed to break particles into their constituent nucleons [MeV] ---
const double deuteronBreakupEnergy = 2.22;
const double tritonBreakupEnergy   = 8.48;
const double he3BreakupEnergy      = 7.72;
const double alphaBreakupEnergy    = 28.3;

//--- Conversion Factors ---
const double amu2ev = 9.31494013e+8; //-- amu -> eV
