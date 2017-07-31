/*
 *  CGMF
 *  Version: 1.0.6
 *
 *  [ FissionFragments.cpp ]
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <sstream>

// To include with XML I/O
//#include <libxml/xpath.h>
//#include <libxml/xpathInternals.h>
//#include "ffInput.hxx"

#include "mt19937ar.h"
#include "physics.h"
#include "masstable.h"

using namespace std;

#include "FissionFragments.h"

// << 1.0.6 >> pre-fission neutron emission

//#include "Yields.h"
#include "excinterface.h"
#include "evapinterface.h"
#include "terminate.h"


// from CGM
#include "config.h"
#include "cgm.h"
#include "ripl2levels.h"

#include "global-mcnp.h"

double **pfnInfo;

FissionFragments::FissionFragments (void) { } // default constructor


/*******************************************************************************
 * Constructor
 *------------------------------------------------------------------------------
 * Define an input filename (.dat) that corresponds to a ZAID and incident
 * energy. Then calls routine init(inputFilename);
 ******************************************************************************/
FissionFragments::FissionFragments (int ZAID, double Einc, double *alphaSpin) {
  
  incidentEnergy  = Einc;
  
  // Target nucleus
  ZAIDt=ZAID;
  Zt=int(ZAIDt/1000.0);
  At=ZAIDt-1000*Zt;
  
  //-- compound (parent, fissioning) nucleus
  if (incidentEnergy==0.0) { // spontaneous fission
    Zc=Zt; Ac=At;
  } else { // neutron-induced fission
    Zc=Zt; Ac=At+1;
  }
  ZAIDc = 1000*Zc+Ac;
  
  //-- symmetric fragment
  Asym=Ac/2; Zsym=Zc/2;
  
  // << 1.0.6 >> set original (before any pre-fission neutron emission) compound fissioning nucleus
  incidentEnergy0 = Einc;
  ZAIDc0=ZAIDc;
  Zc0=Zc;
  Ac0=Ac;
  
  setOptions();
  init();
  
  if (*alphaSpin==0.0) {
    *alphaSpin = alphaI;
  } else {
    alphaI = *alphaSpin;
  }
  
  
  
}

/*******************************************************************************
 * Constructor
 *------------------------------------------------------------------------------
 * Define an input filename (.dat) that corresponds to a ZAID and incident
 * energy. Then calls routine init(inputFilename);
 ******************************************************************************/
FissionFragments::FissionFragments (int ZAID, double Einc, double *alphaSpin, double RTin) {
  
  incidentEnergy  = Einc;
  
  // Target nucleus
  ZAIDt=ZAID;
  Zt=int(ZAIDt/1000.0);
  At=ZAIDt-1000*Zt;
  
  //-- compound (parent, fissioning) nucleus
  if (incidentEnergy==0.0) { // spontaneous fission
    Zc=Zt; Ac=At;
  } else { // neutron-induced fission
    Zc=Zt; Ac=At+1;
  }
  ZAIDc = 1000*Zc+Ac;
  
  //-- symmetric fragment
  Asym=Ac/2; Zsym=Zc/2;
  
  // << 1.0.6 >> set original (before any pre-fission neutron emission) compound fissioning nucleus
  incidentEnergy0 = Einc;
  ZAIDc0=ZAIDc;
  Zc0=Zc;
  Ac0=Ac;
  
  setOptions();

  // override setoptions for the RT parameter
  if (RTin >= 0.5 && RTin <= 1.5){
    energySortingOption="fixedRT";
    RT=RTin;
  }

  init();
  
  if (*alphaSpin==0.0) {
    *alphaSpin = alphaI;
  } else {
    alphaI = *alphaSpin;
  }
  
  
  
}

/*******************************************************************************
 * Constructor
 *------------------------------------------------------------------------------
 * Constructs an instance of FissionFragments based on user input:
 * ZAID:       ZAID number of target nucleus
 * Einc:       incident neutron energy in MeV
 * alphaSpin:  spin distribution factor
 * yieldsFile: file providing the initial fission fragment yields in (A,Z,TKE,Y)
 ******************************************************************************/
FissionFragments::FissionFragments (int ZAID, double Einc, double *alphaSpin, string yieldsFile) {
  
  incidentEnergy=Einc;
  
  // Target nucleus
  ZAIDt=ZAID;
  Zt=int(ZAIDt/1000.0);
  At=ZAIDt-1000*Zt;
  
  //-- compound (parent, fissioning) nucleus
  if (incidentEnergy==0.0) { // spontaneous fission
    Zc=Zt; Ac=At;
  } else { // neutron-induced fission
    Zc=Zt; Ac=At+1;
  }
  ZAIDc = 1000*Zc+Ac;
  
  //-- symmetric fragment
  Asym=Ac/2; Zsym=Zc/2;
  
  // setting options
  yieldsFilename = yieldsFile;
  energySortingOption="RTA";
  alphaI=*alphaSpin;
  
  if(ZAIDt == 98252){
    if (incidentEnergy!=0.0) {
      cerr << "ERROR: only spontaneous fission is treated for Cf-252." << endl;
      exit(-1);
    }
    
    // yieldsOption="YATKE";
    // yieldsFilename += "yieldsATKE.Cf252.sf.dat3";

    yieldsOption = "BrosaYields";
    //yieldsFilename += "BrosaYields.dat"
    dZ=3;
    sigmaZ=0.15;
    if (alphaI==0.0) alphaI=1.7;
    energySortingOption="RTA";
  }else{
    yieldsOption="YAZTKEfile";
    dZ=1;
    sigmaZ=0;
  }
   
  init();
  
}




/*******************************************************************************
 * Constructor
 *------------------------------------------------------------------------------
 * Define an input filename that corresponds to a ZAID and incidentEnergy. Then
 * reads the corresponding XML input file and builds the FF yields.
 ******************************************************************************/
FissionFragments::FissionFragments (int ZAIDf, double excitationEnergy, string option)
{
  
  if (option.compare("single")!=0) { cerr << "Wrong option for FissionFragments\n"; }
  
  // no need to do anything
  
}


/*******************************************************************************
 setOptions for different isotopes/energies
 
 - yieldsOption:        format in which fission fragment yields are given
 - yieldsFilename:      file containing fission fragment yields
 - dZ:                  number of nuclei to consider around most probable Zp
 - sigmaZ:              width of charge distribution in systematics
 - alphaI:              parameter to determine spin cut-off parameter in FF
 - energySortingOption: used to sort out TXE between the two fragments
 
 ******************************************************************************/
void FissionFragments::setOptions () {
  
  yieldsFilename = DATADIR;

  
  switch (ZAIDt)
  {
    
    case 92235: //-- NEUTRON-INDUCED FISSION OF U-235 --------------------------
    
    if (abs(incidentEnergy-2.53e-8)<1e-7) {
      yieldsOption="YATKE"; // keep thermal point outside of systematics
      yieldsFilename += "yieldsATKE.U235.Thermal.dat";
    } else {
      yieldsOption = "Systematics2"; // << 1.0.6 >>
    }
    
    dZ=4;
    sigmaZ=0.15;
    if (alphaI==0.0) alphaI=1.7;
    energySortingOption="RTA";
    
    //energySortingOption="fixedRT";
    //RT=1.0;
   

    //    yieldsFilename += "yieldsATKE.U235.Thermal.dat";
    //    yieldsFilename += "yieldsATKE.U235.0.5MeV.dat";
    //    yieldsFilename += "yieldsATKE.U235.5.5MeV.dat";
    
    break;
    
    case 92238: //-- SPONTANEOUS FISSION OF U-238 ------------------------------
    
    if (incidentEnergy!=0.0) {
      cerr << "ERROR: Only U-238 Spontaneous Fission is treated." << endl;
      exit(-1);
    }
    
    yieldsOption="Systematics2";
    
    // default -- should be reviewed
    dZ=2;
    sigmaZ=0.0;
    if (alphaI==0.0) alphaI=1.0;
    energySortingOption="fixedRT";
    RT=1.0;
    
    break;
    
    cerr << "TO BE IMPLEMENTED IN FIRST OPEN-SOURCE VERSION OF THE CODE!\n";
    exit(0);
    
    
    case 94239: //-- NEUTRON-INDUCED FISSION OF PU-239 -------------------------
    
    if (abs(incidentEnergy-2.53e-8)<1e-7) {
      yieldsOption="YATKE"; // keep thermal data out of systematics
      yieldsFilename += "yieldsATKE.Pu239.Thermal.dat";
    } else {
      yieldsOption="Systematics2"; // << 1.0.6 >>
    }
    
    dZ=2;
    sigmaZ=0.0;
    if (alphaI==0.0) alphaI=1.5;
    energySortingOption="RTA";
    
    break;
    
    
    case 94240: //-- SPONTANEOUS FISSION OF PU-240 -----------------------------
    
    if (incidentEnergy!=0.0) {
      cerr << "ERROR: Only Pu-240 Spontaneous Fission is treated." << endl;
      exit(-1);
    }
    
    //    yieldsOption="Systematics2";
    //
    //    // default -- should be reviewed
    //    dZ=2;
    //    sigmaZ=0.0;
    //    if (alphaI==0.0) alphaI=1.0;
    //    energySortingOption="fixedRT";
    //    RT=1.0;
    
    yieldsOption = "YATKE";
    yieldsFilename += "yieldsATKE.Pu240sf.dat";
    dZ=7;
    sigmaZ=0.15;
    if (alphaI==0.0) alphaI=1.0;
    energySortingOption="fixedRT";
    RT=1.2;
    break;
    
    
    case 98252: //-- CF-252 SPONTANEOUS FISSION --------------------------------
    
    if (incidentEnergy!=0.0) {
      cerr << "ERROR: only spontaneous fission is treated for Cf-252." << endl;
      exit(-1);
    }
    
    // yieldsOption="YATKE";
    // yieldsFilename += "yieldsATKE.Cf252.sf.dat3";

    yieldsOption = "BrosaYields";
    //yieldsFilename += "BrosaYields.dat"
    dZ=7;
    sigmaZ=0.15;
    if (alphaI==0.0) alphaI=1.7;
    energySortingOption="RTA";
    //energySortingOption="fixedRT";
    //RT=1.0;
    break;
    
    //-- IF ISOTOPE/REACTION NOT DEFINED ABOVE, QUIT WITH ERROR MESSAGE --------
    default:
    cerr << "ERROR: Cannot use CGMF (yet) for this isotope or/and incident energy." << endl;
    exit(-1);
    
  }
  
}


/*******************************************************************************
 * Init
 *------------------------------------------------------------------------------
 * First main routine to be called after constructor
 ******************************************************************************/
void FissionFragments::init () {
  
  try {
    fragments = new Nucleus[2];
  }
  catch (bad_alloc) {
    cerr << "Memory allocation error trying to create CGM-type fragments";
    exit(-1);
  }
  
  // define energy grid
  for (int k=0; k<NUME; k++) ldEnergyGrid[k] = k*deltaE;
  
  // define energy grid #2
  for (int k=0; k<NUME2; k++) {
    if (k<11) ldEnergyGrid2[k] = k*1.0;
    if (k>=11 && k<=15) ldEnergyGrid2[k] = 10.0+(k-10)*2.0;
    if (k>15) ldEnergyGrid2[k] = 20.0+(k-15)*5.0;
  }
  
  resetYields ();
  buildYields ();
  
  //readNeutronBindingEnergies();
  readMasses();
  computeNeutronSeparationEnergies();
  
  // To save time in the fission fragment yield sampling loops
  if (energySortingOption == "Pfactors") computePfactors ();
  
  if (energySortingOption == "RTA") readRTAh (Zt,At);
  //  if (energySortingOption == "RTA") readRTAh (Zt,At,"RTA.dat");
  //  if (energySortingOption == "RTA") readRTAh (Zt,At,"RTA.Becker.U235th.dat");
  
  
  //  computePfactors();
  readDeformations();
  
  //  readMollerDeformationEnergies ("deformationEnergies.Moller2012.".ZAIDc.".dat"); /// HERE ////////////////////
  
  //  readLevelDensityParameterTables("ldparamtables.dat");
  //	readLevelDensityTables("ldtables.dat");
  
  
  readTemperatures ();
  readLevelDensityParameters();
  
  /*
   computeLevelDensityParameterTables ();
   
   computeLevelDensityTables(); // rho(Z,A,U)
   computeTemperatureTables();  // T(Z,A,U) // needed for beta parameter in spin-dependent level density
   
   dZ=3;
   for (int i=Amin; i<=Amax; i++) {
   for (int j=-dZ; j<=+dZ; j++) {
			computeTemperatures((int) Z0[i]+j,i);
   }
   }
   
   exit(0);
   */
  
}

/*
 readTemperatures
 */
void FissionFragments::readTemperatures () {
  
  string tempFilename = DATADIR;
  tempFilename += "temperatures.dat";
	  
  ifstream tempFile;
  
  tempFile.open(&tempFilename[0]);
  if (!tempFile) { cerr << "ERROR: Cannot find nuclear temperatures file: " << tempFilename << endl; exit(-1); }
  
  int c=-1;
  while (!tempFile.eof()) {
    c++;
    tempFile >> ZAIDlist[c];
    for (int i=0; i<32; i++) tempFile >> temp[c][i];
  }
  
}

/*
 readLevelDensityParameters
 */
void FissionFragments::readLevelDensityParameters () {
  
  string filename = DATADIR;
  filename += "ldp.dat";
  
  ifstream dataFile;
  
  dataFile.open(&filename[0]);
  if (!dataFile) { cerr << "ERROR: Cannot find level density parameters file: " << filename << endl; exit(-1); }
  
  int c=-1;
  while (!dataFile.eof()) {
    c++;
    dataFile >> ZAIDlist[c];
    for (int i=0; i<32; i++) dataFile >> ldp[c][i];
  }
  
}



/*******************************************************************************
 * Set yields YA[], YZ[], ... to zero.
 ******************************************************************************/
void FissionFragments::resetYields() {
  
  for (int i=0; i<NUMA; i++) {
    YA[i]=0.0;
    maxYieldZ[i]=0.0;
    QfA[i]=0.0;
  }
  
  for (int i=0; i<NUMZ; i++) YZ[i]=0.0;
  
  for (int i=0; i<NUMTKE; i++) {
    YTKE[i]=0.0;
    QfTKE[i]=0.0;
  }
  
  for (int i=0; i<NUMA; i++) {
    for (int j=0; j<NUMTKE; j++) {
      YATKE[i][j] = 0.0;
    }
  }
  
  for (int i=0; i<NUMdZ; i++) {
    for (int j=0; j<NUMA; j++) {
      YZA[i][j] = 0.0;
    }
  }
  
  for (int i=0; i<NUMZ; i++) {
    for (int j=0; j<NUMA; j++) {
      YZA2[i][j] = 0.0;
    }
  }
  
  for (int i=0; i<NUMZ; i++) {
    for (int j=0; j<NUMA; j++) {
      for (int k=0; k<NUMTKE; k++) {
        YZATKE[i][j][k] = 0.0;
      }
    }
  }
  
  
  
}

/*******************************************************************************
 * Reads a main input file without XML formatting.
 ******************************************************************************/
void FissionFragments::readMainInputFile (string inputFilename) {
  
  ifstream inputFile;
  string line, key, value;
  size_t found;
  
  inputFile.open(&inputFilename[0]);
  if (!inputFile) { cerr << "[readMainInputFile] Cannot find main input file!\n"; cerr << inputFilename << "\n"; exit(-1); }
  
  while (!inputFile.eof()) {
    
    inputFile >> line;
    
    found = line.find("=");
    if (found!=string::npos) {
      
      key = line.substr(0, int(found));
      value = line.substr(int(found)+1,80);
      //cout << key << " : " << value << "\n";
      
      if (key == "At") At = atoi(value.c_str());
      if (key == "Zt") Zt = atoi(value.c_str());
      if (key == "Ap") Ap = atoi(value.c_str());
      if (key == "Zp") Zp = atoi(value.c_str());
      
      if (key == "incidentEnergy") incidentEnergy = atof(value.c_str());
      
      if (key == "yieldsOption") yieldsOption = value.c_str();
      if (key == "yieldsFilename") yieldsFilename = value.c_str();
      
      //      if (key == "TKEmin") TKEmin = atoi(value.c_str());
      //      if (key == "TKEmax") TKEmax = atoi(value.c_str());
      
      if (key == "Amin") Amin = atoi(value.c_str());
      if (key == "Amax") Amax = atoi(value.c_str());
      if (key == "Zmin") Zmin = atoi(value.c_str());
      if (key == "Zmax") Zmax = atoi(value.c_str());
      
      if (key == "dZ") dZ = atoi(value.c_str());
      if (key == "sigmaZ") sigmaZ = atof(value.c_str());
      
      if (key == "energySortingOption") energySortingOption = value;
      
      if (key == "RT") RT = atof(value.c_str());
      
      if (key == "alphaI") alphaI=atof(value.c_str());
      
    }
    
  }
  
  
  
  //-- compound (parent, fissioning) nucleus
  if (incidentEnergy==0.0) { // spontaneous fission
    Zc=Zt; Ac=At;
  } else { // neutron-induced fission
    Zc=Zt; Ac=At+1;
  }
  
  //-- symmetric fragment
  Asym=Ac/2; Zsym=Zc/2;
  
  inputFile.close();
  
}

// TO INCLUDE WITH XML I/O
/*******************************************************************************
 * First reads an XML input data file, and creates an instance of FissionFragments.
 * It then builds the fission fragment yields Y(A,Z,TKE).
 *
 * << XCode note >>
 *
 * It requires XERCES and XSD C++ libraries installed and loaded as "external
 * frameworks and libraries". Also, paths to the libraries are to be included
 * in the project header files path search.
 *
 ******************************************************************************/
/* void FissionFragments::readXMLInputFile (string inputFilename) {
 
 try {
 
 auto_ptr<ffInputType> mainInput (ffInput(inputFilename)); // parses XML data file using XSD C++ routines
 
 //-- target nucleus
 Zt = mainInput->target().Z();
 At = mainInput->target().A();
 
 incidentEnergy = mainInput->incidentEnergy();
 
 //-- compound (parent, fissioning) nucleus
 if (incidentEnergy==0.0) { // spontaneous fission
 Zc=Zt; Ac=At;
 } else { // neutron-induced fission
 Zc=Zt; Ac=At+1;
 }
 
 //-- symmetric fragment
 Asym=Ac/2; Zsym=Zc/2;
 
 TKEmin = mainInput->yields().TKEmin();
 TKEmax = mainInput->yields().TKEmax();
 
 Amin = mainInput->yields().Amin();
 Amax = mainInput->yields().Amax();
 
 Zmin = mainInput->yields().Zmin();
 Zmax = mainInput->yields().Zmax();
 
 dZ     = mainInput->yields().dZ();
 sigmaZ = mainInput->yields().sigmaZ();
 // dZFile = mainInput->yields().dZFile();
 
 energySortingOption = mainInput->energySorting().option();
 RT = mainInput->energySorting().RT();
 
 yieldsOption   = mainInput->yields().yieldsOption();
 yieldsFilename = mainInput->yields().yieldsFilename();
 
 }
 
 catch (const xml_schema::exception& e)
 {
 cerr << e << endl;
 exit(1);
 }
 
 
 }
 */


// Destructor
FissionFragments::~FissionFragments () {}



/*!
 Reads deformation energies at scission calculated by P.Moller, 2012.
 */

void FissionFragments::readMollerDeformationEnergies (string filename) {
  
  int A;
  double Edef;
  
  std::fill_n(EdefA, NUMA, 0.0);
  
  ifstream data;
  data.open (&filename[0]);
  
  while (!data.eof()) {
    data >> A >> Edef;
    EdefA[A] = Edef;
  }
  data.close();
  
}

/*******************************************************************************
 * readLevelDensityTables
 ******************************************************************************/
void FissionFragments::readLevelDensityTables (string datafile) {
  
  cout << "[readLevelDensityTables] Running... ";
  
  
  ifstream data;
  string str;
  unsigned pos;
  int A;
  //  int Z, deltaZ;
  
  levelDensities = new double **[NUMA];
  for (int i=Amin; i<=Amax; i++) {
    levelDensities[i] = new double *[NUMdZ];
    for (int j=-dZ; j<=dZ; j++) {
      levelDensities[i][j+dZ] = new double [NUME];
      std::fill_n(levelDensities[i][j+dZ], NUME, 0.0);
    }
  }
  
  // read data file
  data.open (&datafile[0]);
  while (getline(data,str)) {
    
    if (str.find("#")==0) {
      
      pos = (int) str.find("A=");  A  = atoi(str.substr(pos+2,3).c_str());
      //      pos = str.find("Zp="); Z  = atoi(str.substr(pos+3,3).c_str());
      //      pos = str.find("dZ="); deltaZ = atoi(str.substr(pos+3,2).c_str());
      
      if (A>Amax || A<Amin) { continue; }
      
      //			cout << A << " " << Z << " " << deltaZ << "\n";
      
      getline(data,str); getline(data,str); // skip two lines
      
      getline(data,str);
      while (!str.empty()) {
        for (int k=0; k<NUME; k++) {
          //					data >> levelDensityParameters[A][0][k];
          //					cout << k << " " << levelDensityParameters[A][0][k] << "\n";
          for (int j=-dZ; j<=dZ; j++) {
            levelDensities[A][j+dZ][k] = atof(str.substr((j+dZ+1)*10,10).c_str());
            //						cout << setprecision(6) << atof(str.substr((j+dZ+1)*10,11).c_str()) << "\n";
          }
          getline(data,str);
        }
      }
      
      // TO BE CONTINUED... NEED TO SOLVE PROBLEMS OF DIFFERENT ENERGY GRID
      // AND DIFFERENT FRAGMENTS PRODUCED IN DIFFERENT REACTIONS
      
    }
  }
  data.close();
  
  cout << "OK\n";
  
}


/*******************************************************************************
 * computeLevelDensityTables
 *------------------------------------------------------------------------------
 * Computes level density tables for (A,Z) on a given energy grid (with NUME
 * number of energy points), which are later used to compute energy sorting at
 * scission.
 ******************************************************************************/
void FissionFragments::computeLevelDensityTables () {
  
  Nucleus nucleus;
  
  ofstream out;
  string fn;
  fn=WORKDIR+"ldtables.dat";
  stringstream line;
  out.open(&fn[0], ios::out);
  out << setprecision(3) << setiosflags(ios::scientific);
  
  cout << "[computeLevelDensityTables] Running... ";
  
  levelDensities = new double **[NUMA];
  for (int i=Amin; i<=Amax; i++) {
    levelDensities[i] = new double *[NUMdZ];
    for (int j=-dZ; j<=dZ; j++) {
      levelDensities[i][j+dZ] = new double [NUME];
      //    	if (YZA[j+dZ][i]!=0.0) {
      nucleus.za.setZA((int) floor(Z0[i]+j+0.5),i);
      nucleus.ndisc = riplReadDiscreteLevels(&nucleus.za,nucleus.lev,reassign);
      //        nucleus.ndisc = riplReadDiscreteLevels(&nucleus.za,nucleus.lev,normal);
      statFixDiscreteLevels(&nucleus);
      statSetupEnergyBin(&nucleus);
      statSetupLevelDensityParameter(&nucleus,&nucleus.ldp);
      for (int k=0; k<NUME; k++) {
        levelDensities[i][j+dZ][k] = ldLevelDensity(ldEnergyGrid[k], i, &nucleus.ldp);
      }
      
      //      }
    }
    
    out << "\n# A=" << i << " ; Zp=" << int(Z0[i]+0.5) << " ; dZ=" << dZ << "\n";
    out << "# U / Z = ";
    out << setw(10);
    for (int j=-dZ; j<=dZ; j++) { out << setw(10) << int(Z0[i]+0.5)+j; }
    out << endl << endl;
    for (int k=0; k<NUME; k++) {
      out << setw(10) << ldEnergyGrid[k];
      for (int j=-dZ; j<=dZ; j++) {
        out << setw(10) << levelDensities[i][j+dZ][k];
      }
      out << "\n";
    }
    
  }
  
  cout << "OK\n";
  
}


/*******************************************************************************
 * readLevelDensityParameterTables
 ******************************************************************************/
void FissionFragments::readLevelDensityParameterTables (string datafile) {
  
  cout << "[readLevelDensityParameterTables] Running... ";
  
  
  ifstream data;
  string str;
  unsigned pos;
  int A;
  //  int Z, deltaZ;
  
  levelDensityParameters = new double **[NUMA];
  for (int i=Amin; i<=Amax; i++) {
    levelDensityParameters[i] = new double *[NUMdZ];
    for (int j=-dZ; j<=dZ; j++) {
      levelDensityParameters[i][j+dZ] = new double [NUME];
      std::fill_n(levelDensityParameters[i][j+dZ], NUME, 0.0);
    }
  }
  
  // read data file
  data.open (&datafile[0]);
  while (getline(data,str)) {
    
    if (str.find("#")==0) {
      
      pos = (int) str.find("A=");  A  = atoi(str.substr(pos+2,3).c_str());
      //      pos = (int) str.find("Zp="); Z  = atoi(str.substr(pos+3,3).c_str());
      //      pos = (int) str.find("dZ="); deltaZ = atoi(str.substr(pos+3,2).c_str());
      
      if (A>Amax || A<Amin) { continue; }
      
      getline(data,str); getline(data,str); // skip two lines
      
      getline(data,str);
      while (!str.empty()) {
        for (int k=0; k<NUME; k++) {
          //					data >> levelDensityParameters[A][0][k];
          //					cout << k << " " << levelDensityParameters[A][0][k] << "\n";
          for (int j=-dZ; j<=dZ; j++) {
            levelDensityParameters[A][j+dZ][k] = atof(str.substr((j+dZ+1)*10,10).c_str());
            //						cout << setprecision(6) << atof(str.substr((j+dZ+1)*10,11).c_str()) << "\n";
          }
          getline(data,str);
          //					cout << str << "\n";
        }
      }
      
      // TO BE CONTINUED... NEED TO SOLVE PROBLEMS OF DIFFERENT ENERGY GRID
      // AND DIFFERENT FRAGMENTS PRODUCED IN DIFFERENT REACTIONS
      
    }
  }
  data.close();
  
  cout << "OK\n";
  
}


/*******************************************************************************
 * computeLevelDensityParameterTables
 *------------------------------------------------------------------------------
 * Computes level density parameter tables f(A,Z,U), which are later used to
 * compute energy sorting at scission, if energySortingOption = "fixedRT".
 ******************************************************************************/
void FissionFragments::computeLevelDensityParameterTables () {
  
  Nucleus nucleus;
  
  ofstream out;
  string fn;
  fn = WORKDIR + "ldparamtables.dat";
  stringstream line;
  out.open(&fn[0], ios::out);
  out << setprecision(6) << setiosflags(ios::fixed);
  
  cout << "[computeLevelDensityParameterTables] Running... ";
  
  levelDensityParameters = new double **[NUMA];
  for (int i=Amin; i<=Amax; i++) {
    
    levelDensityParameters[i] = new double *[NUMdZ];
    for (int j=-dZ; j<=dZ; j++) {
      
      levelDensityParameters[i][j+dZ] = new double [NUME];
      nucleus.za.setZA((int) Z0[i]+j,i);
      nucleus.ndisc = riplReadDiscreteLevels(&nucleus.za,nucleus.lev,reassign);
      //nucleus.ndisc = riplReadDiscreteLevels(&nucleus.za,nucleus.lev,normal);
      statFixDiscreteLevels(&nucleus);
      statSetupEnergyBin(&nucleus);
      statSetupLevelDensityParameter(&nucleus,&nucleus.ldp);
      for (int k=0; k<NUME; k++) {
        levelDensityParameters[i][j+dZ][k] = ldDensityParameter(ldEnergyGrid[k], i, &nucleus.ldp);
      }
    }
    
    out << "\n# A=" << i << " ; Zp=" << int(Z0[i]+0.5) << " ; dZ=" << dZ << "\n";
    out << "# U / Z = ";
    out << setw(10);
    for (int j=-dZ; j<=dZ; j++) { out << setw(10) << int(Z0[i]+0.5)+j; }
    out << endl << endl;
    for (int k=0; k<NUME; k++) {
      out << setw(10) << ldEnergyGrid[k];
      for (int j=-dZ; j<=dZ; j++) {
        out << setw(10) << levelDensityParameters[i][j+dZ][k];
      }
      out << "\n";
    }
    
  }
  
  out.close();
  
  cout << "OK\n";
  
}

/*******************************************************************************
 * computeEnergyFromMaxEntropy
 *------------------------------------------------------------------------------
 * Compute mean light fragment <Ul> for a given total intrinsic excitation
 * energy, following KHS formula [Phys. Rev. C83, 061601(R) (2011)].
 ******************************************************************************/
double FissionFragments::computeEnergyFromMaxEntropy (int Zl, int Al, double TIXE) {
  // TIXE: Total Intrinsic eXcitation Energy (MeV)
  // (Zl,Al): mass and charge of light fragment
  
  double Ul;
  int Ah, Zh;
  
  Ah=Ac-Al;
  Zh=Zc-Zl;
  
  int dZl, dZh;
  
  dZl = int ( Zl - floor(Z0[Al]+0.5) + dZ );
  dZh = int ( Zh - floor(Z0[Ah]+0.5) + dZ );
  
  int kl, kh;
  double numerator, denominator;
  
  int maxIndex = (int) (TIXE/deltaE)+1;
  
  numerator=0.0;
  denominator=0.0;
  for (kl=1; kl<=maxIndex; kl++) {
    Ul = (kl-1)*deltaE; // deltaE=0.25 (see config-ff.h)
    kh =  (int) (floor((TIXE-Ul)/deltaE) + 1);
    numerator += Ul * levelDensities[Al][dZl][kl] * levelDensities[Ah][dZh][kh];
    denominator += levelDensities[Al][dZl][kl] * levelDensities[Ah][dZh][kh];
  }
  if (denominator!=0) {
    return numerator/denominator;
  } else {
    cout << "[computeEnergyFromMaxEntropy] denominator is null!\n";
    cout << Al << " " << Zl << "\n";
    return 0;
  }
}

/*******************************************************************************
 * studyEnergySorting()
 *------------------------------------------------------------------------------
 ******************************************************************************/
void FissionFragments::studyEnergySorting () {
  int i, count;
  string ldFilename, tempFilename, ctFilename;
  
  ldFilename   = WORKDIR + "levelDensities.CGMF.out";
  tempFilename = WORKDIR + "temperatures.CGMF.out";
  ctFilename   = WORKDIR + "constantTemperatures.CGMF.out";
		
  // LEVEL DENSITIES
  ofstream ldFile;
  ldFile.open (&ldFilename[0]);
  count=-1;
  for (i=Amin; i<Amax; i++) {
    for (int j=-dZ; j<=dZ; j++) {
      ldFile << "# [gnuplot #" << ++count << "] A=" << i << " ; Z=" << (int) floor(Z0[i]+j+0.5) << " (j=" << j << ")\n\n";
      for (int k=0; k<NUME; k++) {
        ldFile << ldEnergyGrid[k] << " " << levelDensities[i][j+dZ][k] << " " << log(levelDensities[i][j+dZ][k]) << "\n";
      }
      ldFile << "\n";
    }
  }
  ldFile.close();
  
  
  // TEMPERATURES
  ofstream tempFile;
  tempFile.open (&tempFilename[0]);
  count=-1;
  for (int i=Amin; i<=Amax; i++) {
    for (int j=0; j<=2*dZ; j++) {
      tempFile << "# [gnuplot #" << ++count << "] A=" << i << " ; Z=" << (int) floor(Z0[i]+j-dZ+0.5) << "\n\n";
      tempFile << setprecision(3) << setiosflags(ios::scientific);
      for (int k=0; k<NUME; k++) {
        tempFile << ldEnergyGrid[k] << " " << temperatures[i][j][k] << " " << levelDensities[i][j][k] << "\n";
      }
      tempFile << "\n";
    }
  }
  tempFile.close();
  
  // LOW-ENERGY TEMPERATURES AS A FUNCTION OF MASS
  ofstream ctFile;
  ctFile.open (&ctFilename[0])  ;
  double t, sum;
  for (int i=Amin; i<=Amax; i++) {
    t=0.0;
    sum=0.0;
    for (int j=0; j<=2*dZ; j++) {
      t += temperatures[i][j][0] * YZA[j][i];
      sum += YZA[j][i];
    }
    if (sum!=0.0) t /= sum;
    ctFile << i << " " << t << "\n";
  }
  ctFile.close();
  
  ofstream maxEntropyFile;
  maxEntropyFile.open ("maxEntropy.CGMF.out");
  // Compute <Ul> from max. entropy principle
  int A1;
  int Z1;
  
  double x, TXE;
  count=-1;
  for (int i=70; i<=Asym; i++) {
    A1 = i;
    Z1 = (int) floor(Z0[i]+0.5);
    maxEntropyFile << "# [gnuplot #" << ++count << "] A=" << i << " ; Z=" << Z1 << "\n\n";
    //	  cout << "# [gnuplot #" << ++count << "] A=" << i << " ; Z=" << Z1 << "\n\n";
    for (int i=1; i<=200; i++) {
      TXE = (i-1)*deltaE;
      x = computeEnergyFromMaxEntropy(Z1, A1, TXE);
      maxEntropyFile << TXE << " " << x << " " << TXE/(1+float(Ac-A1)/float(A1)) << "\n";
    }
    maxEntropyFile << "\n";
  }
  
  maxEntropyFile.close();
  
  cout << "TMP\n";
  exit(0);
  
}

/*******************************************************************************
 * generateInitialFissionFragmentHistories (array out)
 *------------------------------------------------------------------------------
 * Generates Monte Carlo histories of primary fission fragment with their
 * initial excitation energy and intial spin. The number of generated histories
 * is a user input (numberEvents).
 ******************************************************************************/
void FissionFragments::generateInitialFissionFragmentHistories
(fissionFragmentType* lightFragment, fissionFragmentType* heavyFragment, const int numberEvents) {
  
  int i=0 ;
  
  if (yieldsOption=="Systematics2") { // << 1.0.6 >>
    
    Yields ffy (ZAIDc);
    
    while (i<numberEvents) {
      index_save=i;
      constructYields (&ffy); // SUPER SLOWDOWN HERE
      sampleFissionFragments (lightFragment+i, heavyFragment+i, &ffy);
      if (lightFragment[i].U!=0.0 && heavyFragment[i].U!=0.0) i++;
    }
    
  } else if (yieldsOption == "BrosaYields") {

    BrosaYields BROSA(yieldsFilename, ZAIDc, Amin);
    while (i<numberEvents) {
      index_save=i;

      sampleFissionFragmentsBROSA (lightFragment+i, heavyFragment+i, &BROSA);
      if (lightFragment[i].U!=0.0 && heavyFragment[i].U!=0.0) i++;
    }
  }

  else {
    
    while (i<numberEvents) {
      sampleFissionFragments (lightFragment+i, heavyFragment+i);
      if (!(lightFragment[i].U==0.0 && heavyFragment[i].U==0.0)) i++;
    }
    
  }
  
}

/*******************************************************************************
 * generateSingleFissionFragments (array out)
 *------------------------------------------------------------------------------
 * Used to study the decay of only one fission fragment, so it generates a
 * certain number of events, each corresponding to the same fragment (A,Z) and
 * excitation energy. [TESTING]
 ******************************************************************************/
void FissionFragments::generateSingleFissionFragments
(fissionFragmentType* lightFragment, fissionFragmentType* heavyFragment, int ZAIDf, double Eexc, const int numberEvents) {
  
  int Zf = int(ZAIDf/1000);
  int Af = ZAIDf%1000;
  
  for (int i=0; i<numberEvents; i++) {
    lightFragment[i].A = Af;
    lightFragment[i].Z = Zf;
    lightFragment[i].U = Eexc;
    lightFragment[i].KE = 0.0;
    lightFragment[i].spin = 0;
    lightFragment[i].parity = 0;
    heavyFragment[i].A = Af;
    heavyFragment[i].Z = Zf;
    heavyFragment[i].U = Eexc;
    heavyFragment[i].KE = 0.0;
    heavyFragment[i].spin = 0;
    heavyFragment[i].parity = 0;
  }
  
}


/*******************************************************************************
 * generateInitialFissionFragmentHistories (file out)
 *------------------------------------------------------------------------------
 * Generates Monte Carlo histories of primary fission fragment with their
 * initial excitation energy and intial spin. The number of generated histories
 * is a user input (numberEvents).
 ******************************************************************************/
void FissionFragments::generateInitialFissionFragmentHistories (string outputFilename, int numberEvents) {
  fissionFragmentType lightFragment, heavyFragment;
  
  /*
   
   // special format for IAEA CRP, 10/04/2013
   
   cout << "\nIAEA Yields...\n\n";
   
   double ***Yields;
   Yields = new double **[NUMA];
   for (int i=Amin; i<=Amax; i++) {
   Yields[i] = new double *[NUMZ];
   for (int j=0; j<=NUMZ; j++) {
   Yields[i][j] = new double [NUMTKE];
   for (int k=0; k<NUMTKE; k++) Yields[i][j][k]=0.0;
   }
   }
   
   cout << "Number of fission events: " << numberEvents << endl;
   
   for (int i=1; i<=numberEvents; i++) {
   sampleFissionFragments(&lightFragment, &heavyFragment);
   if (lightFragment.U!=0.0 && heavyFragment.U!=0.0) {
   Yields[lightFragment.A][int(lightFragment.Z+0.5)][int(lightFragment.KE+heavyFragment.KE+0.5)]+=1;
   Yields[heavyFragment.A][int(heavyFragment.Z+0.5)][int(lightFragment.KE+heavyFragment.KE+0.5)]+=1;
   //      cout << lightFragment.A << " " << int(lightFragment.Z-Z0[lightFragment.A])+dZ << " " << int(lightFragment.KE+heavyFragment.KE+0.5) << " " <<  Yields[lightFragment.A][int(lightFragment.Z-Z0[lightFragment.A])+dZ][int(lightFragment.KE+heavyFragment.KE+0.5)] << "\n";
   }
   }
   
   ofstream out;
   out.open("YieldsAZTKE.IAEA.dat");
   out.precision(2);
   out << fixed;
   cout << "\n Saving yields... \n\n";
   for (int i=Amin; i<Amax; i++) {
   for (int j=0; j<NUMZ; j++) {
   for (int k=0; k<NUMTKE; k++) {
   if (Yields[i][j][k]!=0.0) out << i << " " << j << " " << k << " " << Yields[i][j][k] << "\n";
   }
   }
   }
   out.close();
   
   return; ////////////////////////////////////////////////////////////////////// <<< TO DELETE
   */
  
  ofstream outputFile;
  outputFilename = WORKDIR+outputFilename;
  outputFile.open(&outputFilename[0]);
  outputFile.precision(2);
  outputFile << fixed;
  
  if (yieldsOption == "Systematics2") { // Use of fission fragment yield systematics
    
    Yields ffy (ZAIDc);
    
    for (int i=1; i<=numberEvents; i++) { // << 1.0.6 >>
      
      index_save=i-1;
      constructYields (&ffy);
      sampleFissionFragments (&lightFragment, &heavyFragment, &ffy);
      
      if (lightFragment.U>0.0 && heavyFragment.U>0.0) {
        
        if (i>1) outputFile << endl;
        
        outputFile << setw(4) << lightFragment.Z << setw(4) << lightFragment.A <<
        setw(10) << lightFragment.KE << setw(8) << lightFragment.U << setw(8) <<
        lightFragment.spin << setw(3) << lightFragment.parity << endl;
        
        outputFile << setw(4) << heavyFragment.Z << setw(4) << heavyFragment.A <<
        setw(10) << heavyFragment.KE << setw(8) << heavyFragment.U << setw(8) <<
        heavyFragment.spin << setw(3) << heavyFragment.parity;
      }
      
    }
    
  } else if (yieldsOption == "BrosaYields"){

      BrosaYields BROSA(yieldsFilename, ZAIDc, Amin);
      
      
      for (int i=1; i<=numberEvents; i++) { // << 1.0.6 >>
      
      index_save=i-1;

      sampleFissionFragmentsBROSA (&lightFragment, &heavyFragment, &BROSA);
      
      if (lightFragment.U>0.0 && heavyFragment.U>0.0) {
        
        if (i>1) outputFile << endl;
        
        outputFile << setw(4) << lightFragment.Z << setw(4) << lightFragment.A <<
        setw(10) << lightFragment.KE << setw(8) << lightFragment.U << setw(8) <<
        lightFragment.spin << setw(3) << lightFragment.parity << endl;
        
        outputFile << setw(4) << heavyFragment.Z << setw(4) << heavyFragment.A <<
        setw(10) << heavyFragment.KE << setw(8) << heavyFragment.U << setw(8) <<
        heavyFragment.spin << setw(3) << heavyFragment.parity;
      }
      
    }
  }

   else { // use of specific fission fragment yields
    
    for (int i=1; i<=numberEvents; i++) {
      
      sampleFissionFragments (&lightFragment, &heavyFragment);
      
      if (lightFragment.U>0.0 && heavyFragment.U>0.0) {
        
        if (i>1) outputFile << endl;
        
        /*				outputFile << setw(4) << lightFragment.Z << setw(4) << lightFragment.A <<
         setw(10) << lightFragment.KE << setw(8) << lightFragment.U << setw(8) <<
         lightFragment.spin << setw(3) << lightFragment.parity << endl;
         
         outputFile << setw(4) << heavyFragment.Z << setw(4) << heavyFragment.A <<
         setw(10) << heavyFragment.KE << setw(8) << heavyFragment.U << setw(8) <<
         heavyFragment.spin << setw(3) << heavyFragment.parity;
         */
        
        outputFile << setw(4) << lightFragment.Z << setw(4) << lightFragment.A <<
        setw(10) << lightFragment.KE << setw(8) << lightFragment.U << setw(8) <<
        lightFragment.spin << setw(3) << lightFragment.parity << endl;
        
        //setw(8) <<
        //				ExpMassMeV[lightFragment.Z][lightFragment.A] << " " << S1n[lightFragment.Z][lightFragment.A] << endl;
        
        outputFile << setw(4) << heavyFragment.Z << setw(4) << heavyFragment.A <<
        setw(10) << heavyFragment.KE << setw(8) << heavyFragment.U << setw(8) <<
        heavyFragment.spin << setw(3) << heavyFragment.parity;
        
        //              setw(8) <<
        //				ExpMassMeV[heavyFragment.Z][heavyFragment.A] << " " << S1n[heavyFragment.Z][heavyFragment.A];
        
      }
      
    }
    
  }
  
  outputFile.close();
  
}

/*******************************************************************************
 * sampleFissionFragments
 *------------------------------------------------------------------------------
 * Performs a Monte Carlo sampling the fission fragment yields Y(A,Z,TKE) to
 * choose a specific triplet {A,Z,TKE}. It then calls the method
 * 'computeFragmentInitialConditions' to return the excitation energy, spin and
 * parity of each fission fragment.
 ******************************************************************************/
void FissionFragments::sampleFissionFragments
(fissionFragmentType *lightFragment, fissionFragmentType *heavyFragment) {
  
  double TKE;
  double maxYz;
  double KEl, KEh;
  int iTKE, iA, iZ=0;
  int j;
  int Al, Zl;
  int Ah, Zh;
  
  double Ul, Uh;       //-- excitation energies
  double spinl, spinh; //-- spins
  int parl, parh;      //-- parities
  
  double y, yz;
  
  double energyRelease;
  
  double *momentum;
  momentum = new double [3]; // (x,y,z) components of fragment momentum
  
  int i=0;
  do {
    
    i++;
    
    //-- first, select A and TKE...
    
    TKE = genrand_real3()*(TKEmax-TKEmin)+TKEmin;
    iTKE = (int) floor(TKE+0.5);
    iA = (int) floor(genrand_real3()*(Amax-Amin)+Amin+0.5);
    
    y = genrand_real3()*maxYieldATKE;
    
    if (y <= YATKE[iA][iTKE]) {
      
      maxYz = 1.01*maxYieldZ[iA];
      j=0;
      do {
        j++;
        //				iZ = int(3*genrand_real3())%3-1;
        iZ=rand()%3-1;
        //          iZ = 0; //-- only 1 Z per fragment mass!!!
        yz = genrand_real3()*maxYz;
        if (j>=10000) { cerr << "Charge sampling exceeded!" << endl; }
        //      } while (yz > YZA[iZ+dZ][iA]);
      } while (yz > YZA[iZ+dZ][iA]);
      
    }
    
    if (i>=10000) { cerr << "[sampleFissionFragments] Sampling did not converge!" << endl; }
    
  } while (y > YATKE[iA][iTKE]);
  
  if (iA<Asym) {
    Al = iA;
    Zl = iZ+int(floor(Z0[iA]+0.5));
  } else {
    Al = Ac-iA;
    Zl = Zc-iZ-int(floor(Z0[iA]+0.5));
  }
  Zh = Zc-Zl;
  Ah = Ac-Al;
  

  
  /*	double maxYield=1.0;
   
   int i=0;
   do {
   i++;
   
   TKE = genrand_real3()*(TKEmax-TKEmin)+TKEmin;
   iTKE = (int) floor(TKE+0.5);
   iA = (int) floor(genrand_real3()*(Amax-Amin)+Amin+0.5);
   iZ = (int) floor(genrand_real3()*2*dZ+0.5)+Z0[iA]-dZ;
   y = maxYield*genrand_real3();
   
   } while (y>YZATKE[iZ][iA][iTKE]);
   
   if (iA<Asym) {
   Al = iA;
   Zl = iZ;
   } else {
   Al = Ac-iA;
   Zl = Zc-iZ;
   }
   Zh = Zc-Zl;
   Ah = Ac-Al;
   */
  
  computeFragmentInitialConditions (TKE, Zl, Al, Zh, Ah, &Ul, &Uh, &spinl,
                                    &spinh, &parl, &parh, &energyRelease, momentum);

  //-- momentum light and heavy fragments

  for (int i=0; i<3; i++) {
    lightFragment->pc[i] = momentum[i];
    heavyFragment->pc[i] = -momentum[i];
  }
  
  //-- total energy of the fragments
  
  lightFragment->p[0] = momentum[0] + Al*amuMeV + Ul;
  heavyFragment->p[0] = momentum[0] + Ah*amuMeV + Uh;
  
  //-- ground-state masses of the fragments
  
  double Ml = Al*amuMeV+mass_excess(Zl,Al);
  double Mh = Ah*amuMeV+mass_excess(Zh,Ah);
  
  lightFragment->mass = Ml;
  heavyFragment->mass = Mh;
  
  //-- fission fragment kinetic energies per nucleon
  //KEl = TKE*float(Ah)/float(Ac);
  //KEh = TKE*float(Al)/float(Ac);
  KEl = TKE*Mh/(Ml+Mh);
  KEh = TKE*Ml/(Ml+Mh);
  
  lightFragment->A = Al;
  lightFragment->Z = Zl;
  lightFragment->N = Al-Zl;
  lightFragment->U = Ul;
  lightFragment->spin = spinl;
  lightFragment->parity = parl;
  lightFragment->KE = KEl;
  lightFragment->TKE = TKE;
  lightFragment->energyRelease = energyRelease;
  
  heavyFragment->A = Ah;
  heavyFragment->Z = Zh;
  heavyFragment->N = Ah-Zh;
  heavyFragment->U = Uh;
  heavyFragment->spin = spinh;
  heavyFragment->parity = parh;
  heavyFragment->KE = KEh;
  heavyFragment->TKE = TKE;
  heavyFragment->energyRelease = energyRelease;
  
		
  // Is there any more compact way of doing this?
  //  lightFragment = {Al, Zl, Al-Zl, Ul, spinl, parl, KEl, TKE, energyRelease};
  //  heavyFragment = fissionFragmentType (Ah, Zh, Ah-Zh, Uh, spinh, parh, KEh, TKE, energyRelease);
  
}

/*******************************************************************************
 * sampleFissionFragments - Version 2 << 1.0.6 >> from Ionel, June 2015
 ******************************************************************************/
void FissionFragments::sampleFissionFragments
(fissionFragmentType *lightFragment, fissionFragmentType *heavyFragment, Yields *ffy) {
  
  double TKE;
  double maxYz;
  double KEl, KEh;
  int iA, iZ=0;
  int j;
  int Al, Zl;
  int Ah, Zh;
  
  double Ul, Uh;       //-- excitation energies
  double spinl, spinh; //-- spins
  int parl, parh;      //-- parities
  
  double y, yz;
  
  double energyRelease;
  
  double *momentum;
  momentum = new double [4]; // 4-vector relativistic energy-momentum
  
  int i=0;
  //-- first, select A and Z...
  do {
    
    i++;
    iA = (int) floor(genrand_real3()*(Amax-Amin)+Amin+0.5);
    y = genrand_real3()*maxYieldA;
    
    if (y <= YA[iA]) {
      
      maxYz = maxYieldZ[iA];
      j=0;
      do {
        j++;
        //iZ = dZ*int((genrand_real3()-1));
        iZ =genrand_int32()%(dZ+1);
        // iZ = 0; //-- only 1 Z per fragment mass!!!
        yz = genrand_real3()*maxYz;
        if (j>=10000) { cerr << "Charge sampling exceeded!" << endl; }
      } while (yz > YZA[iZ+dZ][iA]);
      
    }
    
    if (i>=10000) { cerr << "[sampleFissionFragments] Sampling did not converge!" << endl; }
    
  } while (y > YA[iA]);
  
  
  if (iA<=Ac/2) {
    Al = iA;
    Zl = iZ+int(floor(Z0[iA]+0.5));
  } else {
    Al = Ac-iA;
    Zl = Zc-iZ-int(floor(Z0[iA]+0.5));
  }
  Zh = Zc-Zl;
  Ah = Ac-Al;
  
  do{
    TKE=ffy->sampleTKE(Ah);
    computeFragmentInitialConditions (TKE, Zl, Al, Zh, Ah, &Ul, &Uh, &spinl,
                                      &spinh, &parl, &parh, &energyRelease, momentum);
  }while (Ul==0. && Uh==0.);
	
	//-- momentum light and heavy fragments
	for (int i=0; i<3; i++) {
		lightFragment->pc[i] = momentum[i];
		heavyFragment->pc[i] = -momentum[i];
	}
	
  //-- fission fragment kinetic energies per nucleon
  KEl = TKE*float(Ah)/float(Ac);
  KEh = TKE*float(Al)/float(Ac);
  
  lightFragment->A = Al;
  lightFragment->Z = Zl;
  lightFragment->N = Al-Zl;
  lightFragment->U = Ul;
  lightFragment->spin = spinl;
  lightFragment->parity = parl;
  lightFragment->KE = KEl;
  lightFragment->TKE = TKE;
  lightFragment->energyRelease = energyRelease;
  
  heavyFragment->A = Ah;
  heavyFragment->Z = Zh;
  heavyFragment->N = Ah-Zh;
  heavyFragment->U = Uh;
  heavyFragment->spin = spinh;
  heavyFragment->parity = parh;
  heavyFragment->KE = KEh;
  heavyFragment->TKE = TKE;
  heavyFragment->energyRelease = energyRelease;
  
  // Is there any more compact way of doing this?
  //  lightFragment = {Al, Zl, Al-Zl, Ul, spinl, parl, KEl, TKE, energyRelease};
  //  heavyFragment = fissionFragmentType (Ah, Zh, Ah-Zh, Uh, spinh, parh, KEh, TKE, energyRelease);
  
}


/*******************************************************************************
 * sampleFissionFragments - Version 2 << 1.0.6 >> from Ionel, June 2015
 ******************************************************************************/
void FissionFragments::sampleFissionFragmentsBROSA
(fissionFragmentType *lightFragment, fissionFragmentType *heavyFragment, BrosaYields *BY) {
  
  double TKE;
  double maxYz;
  double KEl, KEh;
  int iA, iZ=0;
  int j;
  int Al, Zl;
  int Ah, Zh;
  
  double Ul, Uh;       //-- excitation energies
  double spinl, spinh; //-- spins
  int parl, parh;      //-- parities
  
  double y, yz;
  
  double energyRelease;
  
  double *momentum;
  momentum = new double [4]; // 4-vector relativistic energy-momentum
  
  Al = 0;
  Ah = 0;
  Zl = 0;
  Zh = 0;
  //-- first, select A and Z...

  int i_mode = BY -> samplemode(); 

do{
  int A_sampled = BY -> sampleA(i_mode);
  if (A_sampled >= Ac/2 ){
    Ah = A_sampled;
    Al = Ac - Ah;
  } else{
    Al = A_sampled;
    Ah = Ac - Al;
  }
  
  Zh = BY -> sampleZA(Ah);
  Zl = Zc - Zh;

} while(Al > 0 && Al < Amax && Ah > 0 && Ah < Amax && Zl > 0 && Zl > Zmax && Zh > 0 && Zh < Zmax);


  if(Zl == 0 || Zh == 0 || Al == 0 || Ah == 0){
    cerr << "Bad Sampling in sampleFissionFragmentsBROSA" << endl;
      exit(-1);
  }
  
  
  do{
    TKE=(double)BY->sampleTKEA(Ah, i_mode);
    //cout << "TKE = " << TKE << endl;
    computeFragmentInitialConditions (TKE, Zl, Al, Zh, Ah, &Ul, &Uh, &spinl,
                                      &spinh, &parl, &parh, &energyRelease, momentum);
  }while (Ul==0. && Uh==0.);
  
  //-- momentum light and heavy fragments
  for (int i=0; i<3; i++) {
    lightFragment->pc[i] = momentum[i];
    heavyFragment->pc[i] = -momentum[i];
  }
  
  //-- fission fragment kinetic energies per nucleon
  KEl = TKE*float(Ah)/float(Ac);
  KEh = TKE*float(Al)/float(Ac);
  
  lightFragment->A = Al;
  lightFragment->Z = Zl;
  lightFragment->N = Al-Zl;
  lightFragment->U = Ul;
  lightFragment->spin = spinl;
  lightFragment->parity = parl;
  lightFragment->KE = KEl;
  lightFragment->TKE = TKE;
  lightFragment->energyRelease = energyRelease;
  
  heavyFragment->A = Ah;
  heavyFragment->Z = Zh;
  heavyFragment->N = Ah-Zh;
  heavyFragment->U = Uh;
  heavyFragment->spin = spinh;
  heavyFragment->parity = parh;
  heavyFragment->KE = KEh;
  heavyFragment->TKE = TKE;
  heavyFragment->energyRelease = energyRelease;
  
  // Is there any more compact way of doing this?
  //  lightFragment = {Al, Zl, Al-Zl, Ul, spinl, parl, KEl, TKE, energyRelease};
  //  heavyFragment = fissionFragmentType (Ah, Zh, Ah-Zh, Uh, spinh, parh, KEh, TKE, energyRelease);
  
}


/*******************************************************************************
 * computeFragmentInitialConditions
 *------------------------------------------------------------------------------
 * Given a total kinetic energy (TKE) and specific fragmentation (Zl,Al,Zh,Ah),
 * this subroutine returns the initial conditions in energy and spin in both
 * fission fragments.
 * Pf(3): (x,y,z) components of one fission fragment momentum
 ******************************************************************************/
void FissionFragments::computeFragmentInitialConditions (
  double TKE, int Zl, int Al, int Zh, int Ah, double *Ul, double *Uh,
  double *spinl, double *spinh, int *parl, int *parh, double *energyRelease,
  double *pf) {
  
  double addToEnergyRelease, totalExcitationEnergy;
  
  //-- compute total excitation energy
  addToEnergyRelease = incidentEnergy + SnCompound;
  *energyRelease = massExcessCompound - ExpMassMeV[Zl][Al] - ExpMassMeV[Zh][Ah];
  totalExcitationEnergy = *energyRelease + addToEnergyRelease - TKE;
  
  if (totalExcitationEnergy<0) {
    
    *Ul=0.0;
    *Uh=0.0;
    
  } else {
    
    //-- excitation energies
    computeFragmentExcitationEnergy (totalExcitationEnergy, Zl, Al, Zh, Ah, Ul, Uh, energySortingOption);
    
    //-- spins
    getInitialSpin (Zl, Al, Ul, spinl);
    getInitialSpin (Zh, Ah, Uh, spinh);
    
    //    cout << "E*: " << Al << " " << Ah << " : " << totalExcitationEnergy << " " << *Ul << " " << *Uh << " " << *spinl << " " << *spinh << "\n";
    
    //-- parities
    *parl=1; if (genrand_real3()<0.5) { *parl=-1; }
    *parh=1; if (genrand_real3()<0.5) { *parh=-1; }
    
  }

  //-- fragment masses
  double Ml = Al*amuMeV+mass_excess(Zl,Al); // + *Ul;
  double Mh = Ah*amuMeV+mass_excess(Zh,Ah); // + *Uh;
  
  // compute momentum vector (px,py,pz) ----------------------------------------
  double Pf0 = sqrt( (2.0*TKE)/(1.0/Ml+1.0/Mh) );
    
  // choose a (random) direction for the light fragment first (isotropic emission)
  double phi = genrand_real2()*twopi;
  double costheta = 2.0*genrand_real2()-1.0;
  double sintheta = sqrt(1.0-costheta*costheta);

  pf[0] = Pf0*sintheta*cos(phi);
  pf[1] = Pf0*sintheta*sin(phi);
  pf[2] = Pf0*costheta;
		
}


void FissionFragments::sortExcitationEnergy (double Utot, int Al, int Ah, double *Ul, double *Uh, double RT) {
  
  double alf, ahf;
  double epsilon;
  
  double El = Utot/2.0;
  double Eh = Utot/2.0;
  
  int il=0; do { il++; } while (ldEnergyGrid[il]<El);
  //    int ih=0; do { ih++; } while (ldEnergyGrid[ih]<Eh);
  int ih = il;
  
  //-- iterative procedure to find (Ul,Uh) -- converges very quickly!
  int count=0;
  do {
    
    alf = levelDensityParameters [Al][dZ][il]; // consider only the most probable Z for a given A
    ahf = levelDensityParameters [Ah][dZ][ih];
    
    Eh = ahf/(alf*RT*RT+ahf)*Utot;
    El = Utot-Eh;
    
    il=0; do { il++; } while (ldEnergyGrid[il]<El);
    ih=0; do { ih++; } while (ldEnergyGrid[ih]<Eh);
    
    epsilon = abs(alf-levelDensityParameters[Al][dZ][il]);
    
    if (++count>10) { cerr << "[sortExcitationEnergy] did not converge!\n"; exit(-1); }
    
  } while (epsilon>0.01);
  
  //    alf = float(Al)/11.;
  //    ahf = float(Ah)/11.;
  
  *Uh = Eh;
  *Ul = El;
  
  return;
}




/*******************************************************************************
 * computeFragmentExcitationEnergy
 *------------------------------------------------------------------------------
 * Given a total excitation energy (TXE) and specific fragmentation (Zl,Al,Zh,Ah),
 * this subroutine computes the initial excitation energies in the light and
 * heavy fragments, respectively. The particular mechanism used for this energy
 * sorting process is controlled by the choice of the keyword 'sortingOption'.
 ******************************************************************************/
void FissionFragments::computeFragmentExcitationEnergy
(double TXE, int Zl, int Al, int Zh, int Ah, double *Ul, double *Uh, string sortingOption) {
  
  double Uint;  //-- total intrinsic energy to be shared between the two fragments
  double Ucoll; //-- total energy stored in collective degrees of freedom at scission
  
  int iZl = int( Zl-floor(Z0[Al]+0.5)+dZ );
  int iZh = int( Zh-floor(Z0[Ah]+0.5)+dZ );
  
  /*
   cout << "TEMP3\n";
   *Ul = TXE / 2.0;
   *Uh = TXE / 2.0;
   return;
   */
  
  if (TXE>ldEnergyGrid[NUME-1]) {
    cerr << "[computeFragmentExcitationEnergy] SHOULD INCREASE NUME --> TXE = ldEnergyGrid[NUME] !\n";
    cout << TXE << "\n";
    TXE=ldEnergyGrid[NUME-1];
  }
  
  if (sortingOption=="maxEntropy") {
    
    double f;
    double Utot;
    
    // 1. Read Moller's deformation energies
    // 2. Uint = Utot - EdefLF - EdefHF
    // 3. Share Uint according to maximum entropy ==> UintLF, UintHF
    // 4. Ul = UintLF + EdefLF and Uh = UintHF + EdefHF
    
    f=0.5;
    Utot = TXE - incidentEnergy - SnCompound;
    Ucoll = f * Utot;
    Uint = (1-f) * Utot + incidentEnergy + SnCompound;
    
    cout << TXE << " " << Utot << " " << Ucoll << " " << Uint << "\n";
    exit(0);
    
    //-- Cf252sf
    if (ZAIDc==98252) {
      if (abs(Ah-Asym)>=25) {
        Ucoll = 0.3*TXE; //-- assume that x % of TXE is stored in collective energy
      } else {
        Ucoll = 0.6*TXE;
      }
    }
    
    //-- n+U235
    if  (ZAIDc==92236) {
      if (abs(Ah-Asym)>=25) {
        Ucoll = 0.2*TXE; //-- assume that x % of TXE is stored in collective energy
      } else {
        Ucoll = 0.75*TXE;
      }
    }
    
    //-- n+Pu239
    if  (ZAIDc==94240) {
      if (abs(Ah-Asym)>=20) {
        Ucoll = 0.7*TXE; //-- assume that x % of TXE is stored in collective energy
      } else {
        Ucoll = 0.7*TXE;
      }
    }
    
    Uint=TXE-Ucoll;
    
    //    Uint = TXE; // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    /*	  fragments[0].za.setZA(Zl,Al);
     fragments[1].za.setZA(Zh,Ah);
     
     fragments[0].max_energy = TXE;
     fragments[1].max_energy = TXE;
     
     for (int i=0; i<2; i++) {
     fragments[i].ndisc = riplReadDiscreteLevels(&fragments[i].za,fragments[i].lev,reassign);
     statFixDiscreteLevels(&fragments[i]);
     statSetupEnergyBin(&fragments[i]);
     statSetupLevelDensityParameter(&fragments[i],&fragments[i].ldp);
     }
     */
    // Compute from max. entropy principle
    double nominator = 0.0;
    double denominator = 0.0;
    double El, Eh;
    double x;
    
    for (int k=1; k<=NUME; k++){
      
      if (ldEnergyGrid[k]>TXE) { break; }
      
      //    	El = fragments[0].excitation[k];
      El = ldEnergyGrid[k];
      Eh = TXE - El;
      
      // linear interpolation
      double yl, yh;
      double y1, y2, x1, x2;
      
      int il=0;
      do {
        il++;
      } while (ldEnergyGrid[il]<El);
      
      y1 = log(levelDensities[Al][iZl][il-1]);
      y2 = log(levelDensities[Al][iZl][il]);
      x1 = ldEnergyGrid[il-1];
      x2 = ldEnergyGrid[il];
      yl = (y2-y1)/(x2-x1)*(El-x1)+y1;
      yl = exp(yl);
      
      int ih=0;
      do {
        ih++;
      } while (ldEnergyGrid[ih]<Eh);
      
      y1 = log(levelDensities[Ah][iZh][ih-1]);
      y2 = log(levelDensities[Ah][iZh][ih]);
      x1 = ldEnergyGrid[ih-1];
      x2 = ldEnergyGrid[ih];
      yh = (y2-y1)/(x2-x1)*(Eh-x1)+y1;
      yh = exp(yh);
      
      //    	x = ldLevelDensity(El, Al, &fragments[0].ldp)*ldLevelDensity(Eh, Ah, &fragments[1].ldp);
      
      x = yl*yh;
      nominator += El * x;
      denominator += x;
      
    }
    
    
    if (denominator!=0) {
      *Ul = nominator/denominator;
      *Uh = Uint - *Ul;
    } else {
      *Ul = 0.0;
      *Uh = 0.0;
    }
    
    return;
    
  } else if (sortingOption=="fixedRT" || sortingOption=="RTfile" ||
             sortingOption=="Pfactors" || sortingOption=="RTA") {
    
    double alf, ahf;
    
    /*    fragments[0].za.setZA(Zl,Al);
     fragments[1].za.setZA(Zh,Ah);
     
     fragments[0].max_energy = TXE;
     fragments[1].max_energy = TXE;
     
     for (int i=0; i<2; i++) {
     fragments[i].ndisc = riplReadDiscreteLevels(&fragments[i].za,fragments[i].lev,reassign);
     statFixDiscreteLevels(&fragments[i]);
     statSetupEnergyBin(&fragments[i]);
     statSetupLevelDensityParameter(&fragments[i],&fragments[i].ldp);
     }
     */
    
    if (sortingOption == "Pfactors") {
      if (Pfactors[Ah][Zh]!=0) {
        RT = 1+0.05*(Pfactors[Al][Zl]/Pfactors[Ah][Zh]-1.0);
      } else {
        RT = 1.0;
      }
    }
    
    if (sortingOption == "RTA") RT = RTAh [Ah];
    
    double El = TXE/2.0;
    double Eh = TXE/2.0;
    
    //    int il=0; do { il++; } while (ldEnergyGrid[il]<El);
    //    int ih=0; do { ih++; } while (ldEnergyGrid[ih]<Eh);
    //    int ih = il;
    
    double epsilon, alf0;
    
    int il, ih;
    int ZAID;
    int iel, ieh;
    
    ZAID = 1000*Zl+Al;
    il=0; while (ZAID!=ZAIDlist[il]) { il++; }
    ZAID = 1000*Zh+Ah;
    ih=0; while (ZAID!=ZAIDlist[ih]) { ih++; }
    
    iel = findEnergyIndex(El, ldEnergyGrid2);
    ieh = findEnergyIndex(Eh, ldEnergyGrid2);
    
    alf = (ldp[il][iel+1]-ldp[il][iel])/(ldEnergyGrid2[iel+1]-ldEnergyGrid2[iel])*(El-ldEnergyGrid2[iel])+ldp[il][iel];
    ahf = (ldp[ih][ieh+1]-ldp[ih][ieh])/(ldEnergyGrid2[ieh+1]-ldEnergyGrid2[ieh])*(Eh-ldEnergyGrid2[ieh])+ldp[ih][ieh];
    
    //-- iterative procedure to find (Ul,Uh) -- converges very quickly!
    int count=0;
    do {
      
      //      alf = levelDensityParameters [Al][iZl][il];
      //      ahf = levelDensityParameters [Ah][iZh][ih];
      
      Eh = ahf/(alf*RT*RT+ahf)*TXE;
      El = TXE-Eh;
      
      iel = findEnergyIndex(El, ldEnergyGrid2);
      ieh = findEnergyIndex(Eh, ldEnergyGrid2);
      
      alf0=alf;
      
      alf = (ldp[il][iel+1]-ldp[il][iel])/(ldEnergyGrid2[iel+1]-ldEnergyGrid2[iel])*(El-ldEnergyGrid2[iel])+ldp[il][iel];
      ahf = (ldp[ih][ieh+1]-ldp[ih][ieh])/(ldEnergyGrid2[ieh+1]-ldEnergyGrid2[ieh])*(Eh-ldEnergyGrid2[ieh])+ldp[ih][ieh];
      
      epsilon=alf-alf0;
      //      epsilon = abs(alf-levelDensityParameters[Al][iZl][il]);
      
      if (++count>10) { cerr << "[computeFragmentExcitationEnergy]did not converge!\n"; exit(-1); }
      
    } while (epsilon>0.01);
    
    //    alf = float(Al)/11.;
    //    ahf = float(Ah)/11.;
    
    *Uh = Eh;
    *Ul = El;
    
  } else {
    
    cerr << "sorting option not implemented yet!" << endl;
    
  }
  
  //write(77,'(2(i3,2x),2f12.5)') Ah, Al, Uh, Ul
  
}

/*******************************************************************************
 * readRTAh
 *------------------------------------------------------------------------------
 * Set tables R_T(Ah) for a given reaction.
 * Numbers obtained from Bjorn Becker (email on 5/9/2012), adjusted to CGM+FFD,
 * not CGMF!
 ******************************************************************************/
void FissionFragments::readRTAh (int Zt, int At) {
  
  int ZAIDt;
  ZAIDt = 1000*Zt+At;
  
  int j;
  
  std::fill_n (RTAh, NUMA, 1.0);
  
  switch (ZAIDt) {
      
    case 92235:
      
      // from Ah=118 (sym) to Ah=166
      static double U235_RTAh [49] = {
        1.000, 1.096, 1.191, 1.287, 1.382, 1.478, 1.574, 1.669, 1.588, 1.508,
        1.427, 1.346, 1.265, 1.184, 1.192, 1.200, 1.207, 1.215, 1.223, 1.231,
        1.238, 1.235, 1.231, 1.228, 1.225, 1.221, 1.218, 1.214, 1.211, 1.208,
        1.204, 1.201, 1.197, 1.194, 1.191, 1.187, 1.184, 1.180, 1.177, 1.174,
        1.170, 1.167, 1.163, 1.160, 1.157, 1.153, 1.150, 1.146, 1.143 };
      
      j=-1;
      for (int i=118; i<=166; i++) {
        RTAh[i] = U235_RTAh[++j];
      }
      
      break;
      
    case 94239:
      
      // from Ah=120 (sym) to Ah=170
      /*			static double Pu239_RTAh [51] = {
       1.000, 1.044, 1.116, 1.205, 1.298, 1.384, 1.456, 1.506, 1.530, 1.525,
       1.488, 1.421, 1.326, 1.205, 1.105, 1.105, 1.110, 1.113, 1.112, 1.110,
       1.105, 1.098, 1.088, 1.075, 1.061, 1.043, 1.024, 1.002, 0.977, 0.950,
       0.920, 0.888, 0.854, 0.817, 0.778, 0.736, 0.736, 0.736, 0.736, 0.736,
       0.736, 0.736, 0.736, 0.736, 0.736, 0.736, 0.736, 0.736, 0.736, 0.736,
       0.736 };
       */
      
      // modified on 8/22/2012
      static double Pu239_RTAh [51] = {
        1.000, 1.044, 1.116, 1.205, 1.298, 1.384, 1.456, 1.506, 1.530, 1.525,
        1.458, 1.401, 1.306, 1.225, 1.165, 1.135, 1.120, 1.105, 1.100, 1.095,
        1.090, 1.085, 1.080, 1.075, 1.061, 1.043, 1.024, 1.002, 0.977, 0.950,
        0.920, 0.888, 0.854, 0.817, 0.778, 0.736, 0.736, 0.736, 0.736, 0.736,
        0.736, 0.736, 0.736, 0.736, 0.736, 0.736, 0.736, 0.736, 0.736, 0.736,
        0.736 };
      
      j=-1;
      for (int i=120; i<=170; i++) {
        RTAh[i] = Pu239_RTAh[++j];
      }
      
      break;
      
    case 98252:
      
      // from Ah=126 (sym) to Ah=169
      static double Cf252_RTAh [44] = {
        1.000, 1.173, 1.346, 1.519, 1.692, 1.613, 1.534, 1.455, 1.376, 1.297,
        1.307, 1.317, 1.328, 1.338, 1.348, 1.359, 1.328, 1.298, 1.267, 1.237,
        1.206, 1.176, 1.146, 1.115, 1.085, 1.054, 1.051, 1.018, 0.986, 0.953,
        0.967, 0.981, 0.996, 1.010, 1.024, 1.038, 1.029, 1.020, 1.011, 1.002,
        0.993, 0.984, 0.975, 0.496 };
      
      j=-1;
      for (int i=126; i<=169; i++) {
        RTAh[i] = Cf252_RTAh[++j];
      }
      
      break;
      
    default:
      cerr << "ERROR: [readRTAh] does not have data for this isotope\n";
      exit(-1);
  }
  
}

/*********************************************************************************
 * Read R_T(Ah) from file.
 ********************************************************************************/
void FissionFragments::readRTAh (int Zt, int At, string filename) {
  
  //int ZAIDt;
  //ZAIDt = 1000*Zt+At;
  
  int Ah;
  double RT;
  
  std::fill_n (RTAh, NUMA, 1.0);
  
  ifstream dataFile;
  dataFile.open(&filename[0]);
  
  while (!dataFile.eof()) {
    dataFile >> Ah >> RT;
    RTAh[Ah] = RT;
  }
  
  dataFile.close();
  
}

/*******************************************************************************
 * getInitialSpin
 *------------------------------------------------------------------------------
 * Given a nucleus (Z,A), an initial excitation energy U, and a spin
 * parameter B, returns a spin value sampled from a discrete spin distribution
 * P(J) propto (2J+1) exp [ -(J(J+1))/(2*B^2) ].
 ******************************************************************************/
void FissionFragments::getInitialSpin (int Z, int A, double *U, double *spin) {
  
  int const numberSpins = 50;
  double spinProba[numberSpins];
  double oddSpin;
  int j;
  double xj;
  
  double B = spinParameter(*U, A, Z);
  
  /*  if (A<Asym) {
   B = 6.0;  // light fragment B parameter
   } else {
   B = 7.0;  // heavy fragment B parameter
   }
   */
  
  oddSpin=0.0;
  if (A%2!=0) { oddSpin=0.5; }
  
  double sum=0.0;
  for (j=0; j<numberSpins; j++) {
    xj = float(j)+oddSpin;
    spinProba[j] = (2*xj+1)*exp(-(xj*(xj+1))/(2*B*B));
    sum += spinProba[j];
  }
  
  for (j=0; j<numberSpins; j++) spinProba[j] /= sum;
  
  do  {
    j = int(genrand_real3()*numberSpins);
  } while (genrand_real3()>spinProba[j]);
  
  *spin = float(j)+oddSpin;
  
}

/*******************************************************************************
 * writeYieldsInFile
 *------------------------------------------------------------------------------
 * Writes computed fission fragment yields in output file.
 ******************************************************************************/
void FissionFragments::writeYieldsInFile (string outputFilename) {
  
  ofstream outputFile;
  
  outputFile.open(&outputFilename[0]);
  
  outputFile << "# Ac = " << Ac << " ; Zc = " <<Zc << endl;
  
  outputFile << endl << "# Mass Distribution " << endl << endl;
  for (int i=Amin; i<=Amax; i++) {
    outputFile << i << " " << YA[i] << endl;
  }
  outputFile << endl << "# Charge Distribution " << endl << endl;
  for (int i=Zmin; i<=Zmax; i++) {
    outputFile << i << " " << YZ[i] << endl;
  }
  outputFile << endl << "# TKE Distribution " << endl << endl;
  for (int i=TKEmin; i<=TKEmax; i++) {
    outputFile << i << " " << YTKE[i] << endl;
  }
  
  outputFile.close();
  
}

/*******************************************************************************
 * buildYields
 *------------------------------------------------------------------------------
 * Reconstruct fission fragment yields Y(A,Z,TKE) according to 'buildOption',
 * which for now is only 'YATKE'.
 ******************************************************************************/
void FissionFragments::buildYields () {
  
  int i,j;

  
  // << 1.0.6 >>
  SnCompound = -mass_excess(Zc,Ac)+mass_excess(Zc,Ac-1)+mass_excess(0,1);
  if (incidentEnergy==0.0) SnCompound=0.0;
  SnCompound0=SnCompound;
  
  if (yieldsOption == "YATKE") {
    
    readYieldsATKE (yieldsFilename);
    buildYA ();
    computeWahlParametersForChargeDistribution(); // replaces buildZp(string)
    buildYZA ();
    buildYZ ();
    buildYTKE ();
    
  } else if (yieldsOption == "Systematics") {
    
    //    cout << "Yields built from Systematics\n";
    
    computeWahlParametersForChargeDistribution(); // replaces buildZp(string)
    buildWahlYA();
    buildYATKEfromSystematics();
    buildYZA();
    buildYZ();
    buildYTKE();
    
  } else if (yieldsOption == "Systematics2") { // << 1.0.6 >> -----------------------------------------------
    
    //    cout << "Yields built from Systematics-2" << endl;
    
    Amin=20;
    Amax=210;
    computeWahlParametersForChargeDistribution(); // replaces buildZp(string)]
    
    /*** setup discrete levels, energy bins and level densities */
    ncl[0].za.setZA(Zc, Ac);
    ncl[0].max_energy=incidentEnergy+SnCompound;
    nemit = statTotalNeutrons(ncl[0].max_energy,&ncl[0].za);
    //cout<<"Ex="<<ncl[0].max_energy<< endl;
    
    // read the pre-fission neutron emission probability
    // TODO: should be called only for energies above Sn. For now, just remove it
    // for spontaneous fission reactions
    if (incidentEnergy!=0.0) readMultichanceFissionData();
    
    
    try {
      tc=new Transmission[MAX_ENERGY_BIN];
      energy_grid1=new double[MAX_ENERGY_BIN];
      for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tc[k].tran = new double [3*MAX_J];
    }
    catch(bad_alloc){
      cgmTerminateCode("memory allocation error in excinterface");
    }
    
    /*** initialize system parameters */
    statSetupInitSystem(nemit,pdt);
    
    for(int i=0 ; i<=nemit ; i++){
      //ncl[i].ndisc = riplReadDiscreteLevels(&ncl[i].za,ncl[i].lev,normal);
      ncl[i].ndisc = riplReadDiscreteLevels(&ncl[i].za,ncl[i].lev,reassign);
      statFixDiscreteLevels(&ncl[i]);
      statSetupEnergyBin(&ncl[i]);
      statSetupLevelDensityParameter(&ncl[i],&ncl[i].ldp);
      statClearPopulation(&ncl[i]);
    }
    
    sp1 = new double[MAX_ENERGY_BIN]; // evaporation spectrum from first compound nucleus
    sp  = new double[MAX_ENERGY_BIN];
    sp_pfn = new double [MAX_ENERGY_BIN];
    fill_n(sp_pfn, MAX_ENERGY_BIN, 0.0);
    
    evapInterface (pdt+neutron,tc,0,0,sp1);
    
    for(int i=0;i<MAX_ENERGY_BIN;i++){
      energy_grid1[i]=tc[i].ecms;
    }
    
    // pre-equilibrium probability (calculated using CoH-3?) << 1.0.6 >>
    // simplify digits!!
    double f_pe=1.0/(1.0+ exp((12.4924628980014-incidentEnergy)/10.2096843664751)) -0.00422967345976505*incidentEnergy-0.249024419940462;
    if(f_pe<0.)f_pe=0.;
    
    //    cout << "f_pe=" << f_pe << endl;
    
    ofstream sp_out;
    sp_out.open(&(WORKDIR+"sp1.dat")[0]);
    sp_out << scientific;
    
    if(f_pe>0.){
      double *spc[MAX_CHANNEL];
      for (int j=0; j<MAX_CHANNEL; j++) {
        spc[j]=new double[MAX_ENERGY_BIN];
      }
      excitonInterface(Zc, Ac-1, incidentEnergy,spc);
      for (int i=0; i<MAX_ENERGY_BIN; i++) {
        sp_out << energy_grid1[i] << " " << f_pe*spc[neutron][i] <<" " << (1.-f_pe)*sp1[i] << endl;
        sp1[i]=(1.-f_pe)*sp1[i]+f_pe*spc[neutron][i];
      }
      for (int j=0; j<MAX_CHANNEL; j++) {
        delete [] spc[j];
      }
    }else{
      for (int i=0; i<MAX_ENERGY_BIN; i++) {
        sp_out << energy_grid1[i] << " " << sp1[i] << endl;
      }
    }
    
    sp_out.close();
    
    // Wasn't this done already? ............... << 1.0.6 >>
    
    ncl[0].za.setZA(Zc, Ac);
    ncl[0].max_energy=incidentEnergy+SnCompound;
    /*** initialize system parameters */
    statSetupInitSystem(nemit,pdt);
    
    for(int i=0 ; i<=nemit ; i++){
      //ncl[i].ndisc = riplReadDiscreteLevels(&ncl[i].za,ncl[i].lev,normal);
      //cout << "prob["<<i<<"]="<<emissprob[i]<<endl;
      ncl[i].ndisc = riplReadDiscreteLevels(&ncl[i].za,ncl[i].lev,reassign);
      statFixDiscreteLevels(&ncl[i]);
      statSetupEnergyBin(&ncl[i]);
      statSetupLevelDensityParameter(&ncl[i],&ncl[i].ldp);
      statClearPopulation(&ncl[i]);
    }
    
  } else if (yieldsOption == "YAZTKEfile") {
    
    readYieldsAZTKE (yieldsFilename);
    buildYA ();
    buildYZ ();
    buildYTKE ();
    
  } 
  else if (yieldsOption == "BrosaYields") {

    Amin=20;
    Amax=210;
    computeWahlParametersForChargeDistribution(); // replaces buildZp(string)]
    // do nothing
    return;
    
  } 
  else {
    
    cerr << "Fission fragment build option not recognized!" << endl;
    exit(-1);
    
  }
  
  // In the case of "Systematics2", the yields are constructed "on-the-fly"
  // for each fission event because of the sampling of the pre-fission neutrons
  if( yieldsOption == "Systematics2") return;
  
  
  // Symmetrize
  for (i=Amin; i<Asym; i++) {
    YA[Ac-i] = YA[i];
    for (j=Zmin; j<=Zmax; j++) {
      YZA2[Zc-j][Ac-i] = YZA2[j][i];
    }
    for (j=0; j<=2*dZ; j++) {
      YZA[2*dZ-j][Ac-i] = YZA[j][i];
    }
  }
  
  //... and renormalize...?
  
  // determine maxYieldZ[A]
  for (i=0; i<NUMA; i++) {
    maxYieldZ[i]=0.0;
    for (j=0; j<NUMZ; j++) {
      if (YZA2[j][i]>maxYieldZ[i]) { maxYieldZ[i] = YZA2[j][i]; }
    }
  }
  
  //-- determine maximum of Y(A,TKE)
  maxYieldATKE=0.0;
  for (int i=Amin; i<=Amax; i++) {
    for (int j=TKEmin; j<=TKEmax; j++) {
      if (YATKE[i][j]>maxYieldATKE) { maxYieldATKE = YATKE[i][j]; }
    }
  }
  
}

/*******************************************************************************
 * buildYA
 *------------------------------------------------------------------------------
 * Build fission fragments mass-dependent yields Y(A) from Y(A,TKE).
 ******************************************************************************/
void FissionFragments::buildYA (void) {
  float sum;
  sum=0.0;
  for (int i=0; i<NUMA; i++) {
    YA[i] = 0.0;
    for (int j=0; j<NUMTKE; j++) {
      YA[i] = YA[i] + YATKE[i][j];
    }
    sum=sum+YA[i];
  }
  if (sum!=0) {
    for (int i=0; i<NUMA; i++) {
      YA[i] /= sum;
    }
  }
}

/*******************************************************************************
 * buildYATKE << 1.0.6 >>
 *******************************************************************************
 * Builds the Y(A,TKE) distribution
 *******************************************************************************/
void FissionFragments::buildYATKE (Yields *ffy) {
  
  TKEmax=-1; TKEmin=500;
  for (int a=Amin; a<Amax; ++a) {
    for (int t=0; t<NUMTKE; ++t) {
      double y=ffy->yieldATKE(a, (double) t);
      if(y>1.e-7){
        YATKE[a][t]=y;
        if (t<TKEmin) {
          TKEmin=t;
        }
        if (t>TKEmax) {
          TKEmax=t;
        }
      }else{
        YATKE[a][t]=0.;
      }
    }
  }
}


/*******************************************************************************
 * buildYZA
 *------------------------------------------------------------------------------
 * Builds the fission fragment yields YZA[Z][A]. The charge distribution uses
 * Wahl's systematics with parameters obtained from computeWahlParametersForChargeDistribution().
 ******************************************************************************/
void FissionFragments::buildYZA (void) {
  double sigma, c=0.0, Zp0;
  int iZp0, iZ;
  int i, j;
  double oddeven=1.0, V, W;
  
  sigma = 0.4;
  //  c=2.0*(sigma*sigma+1.0/12.0);
  
  for (i=Amin; i<=Amax; i++) {
    c = sZ0[i];
    Zp0 = Z0[i];
    iZp0 = (int) floor(Zp0+0.5);
    
    for (j=0; j<=2*dZ; j++) {
      
      //      Zi = Zp0+j-dZ;
      iZ = iZp0+j-dZ;
      YZA[j][i] = 0.0;
      YZA2[iZ][i] = 0.0;
      
      if (c>0.0) {
        //        Zi = double(int(Zp0+0.5)+j);
        if (iZ%2==0) {
          if (i%2==0) {
            oddeven = FZZ0[i]*FNZ0[i]; // even Z - even N
          } else {
            oddeven = FZZ0[i]/FNZ0[i]; // even Z - odd N
          }
        } else {
          if ((i+1)%2==0) {
            oddeven = FNZ0[i]/FZZ0[i]; // odd Z - even N
          } else {
            oddeven = 1.0/(FZZ0[i]*FNZ0[i]); // odd Z - odd N
          }
        }
        
        V = (iZ-Zp0+0.5)/(sqrt(2.0)*c);
        W = (iZ-Zp0-0.5)/(sqrt(2.0)*c);
        
        YZA[j][i] = oddeven * 0.5 * ( erf(V) - erf(W) ); // YA[i]
        YZA2[iZ][i] = YZA[j][i];
        
        //   ExpYields%Yza(j,i) = oddeven * 1.0_rk/sqrt(c*pi)*exp(-1.0_rk* ( Zi-Zp )**2/c) * 1.0_rk !ExpYields%Ya(i)
        
      } else {
        
        if (j==dZ) {
          //          YZA[j][i] = 1.0; // yields%Ya(i)
          //        	YZA2[iZ][i] = 1.0;
          YZA[j][i] = YA[i];
          YZA2[iZ][i] = YA[i];
        }
      }
      
    } // end loop over Z
    
    
    // TEMPORARY ========================================================================================================
    //		exit(0);
    // TEMPORARY ========================================================================================================
    
    double x=0;
    for (j=0; j<=2*dZ; j++) { x += YZA[j][i]; }
    if (x!=0) {
      for (j=0; j<=2*dZ; j++) {
        YZA[j][i] = YZA[j][i] * YA[i] / x;
      }
    }
    
    double y=0;
    for (j=0; j<NUMZ; j++) { y += YZA2[j][i]; }
    if (y!=0) {
      for (j=0; j<NUMZ; j++) {
        YZA2[j][i] = YZA2[j][i] * YA[i] / y;
      }
    }
    
  } // end loop over A
  
  /*  ofstream out;
   out.open("yza");
   for (i=Amin; i<=Amax; i++) {
   out << "\n#\n\n";
   int z = int(Z0[i]+0.5)-dZ;
   for (j=0; j<=2*dZ; j++) {
			out << i << " " << z+j << " " << YZA[j][i] << "\n";
   }
   }
   out.close();
   */
  
}


/*******************************************************************************
 * buildYZ
 *------------------------------------------------------------------------------
 * Builds the 1D distribution Y(Z) from the Y(Z,A) distribution. It is assumed
 * that Y(Z,A) has been calculated before entering this subroutine.
 ******************************************************************************/
void FissionFragments::buildYZ (void) {
  int i,j;
  double sum;
  
  Zpmin=100;
  Zpmax=0;
  
  //-- Determine min and max values of Zp(A) for any A
  for (i=Amin; i<=Amax; i++) {
    if (YA[i]!=0.0 && Z0[i]<Zpmin) { Zpmin = (int) floor(Z0[i]+0.5); }
    if (YA[i]!=0.0 && Z0[i]>Zpmax) { Zpmax = (int) floor(Z0[i]+0.5); }
  }
  
  //	cout << Amin << " " << Zmin << " " << Z0[Amin] << " " << dZ << " " << Zmin-Z0[Amin]+dZ << endl;
  //	cout << Amax << " " << Zmax << " " << Z0[Amax] << " " << dZ << " " << Zmax-Z0[Amax]-dZ << endl;
  
  if (Zmin>Zpmin-dZ) { cout << "ERROR [builYZ] Zmin > Zp(Amin)-dZ"; exit(-1); }
  if (Zmax<Zpmax+dZ) { cout << "ERROR [builYZ] Zmax > Zp(Amax)+dZ"; exit(-1); }
  
  sum=0.0;
  for (i=Zmin; i<=Zmax; i++) {
    YZ[i]=0.0;
    sum=0;
    for (j=Amin; j<=Amax; j++) {
      YZ[i] += YZA2[i][j];
    }
    sum += YZ[i];
  }
  if (sum!=0) {
    for (i=Zmin; i<=Zmax; i++) {
      YZ[i] /= sum;
    }
  }
  
}

/*******************************************************************************
 * buildYTKE
 *------------------------------------------------------------------------------
 * Build the 1D distribution yields%Ytke(TKEmin:TKEmax) from the distribution
 * Y(A,TKE), assumed to have been calculated earlier.
 ******************************************************************************/
void FissionFragments::buildYTKE (void) {
  int i,j;
  double sum = 0.0;
  for (i=TKEmin; i<=TKEmax; i++) {
    YTKE[i] = 0.0;
    for (j=Amin; j<=Amax; j++) {
      YTKE[i] += YATKE[j][i];
    }
    sum += YTKE[i];
  }
  for (i=TKEmin; i<=TKEmax; i++) { // normalize to unity
    YTKE[i] /= sum;
  }
}

/*******************************************************************************
 * buildYATKEfromSystematics
 *------------------------------------------------------------------------------
 * Build the yields Y[A][TKE] using systematics for <TKE>(A) and <sigma_TKE>(A).
 * It follows the form: Y(A,TKE) = Y(A)xP(TKE|A) where P(TKE) is a Gaussian
 * function whose parameters are mass-dependent.
 ******************************************************************************/
void FissionFragments::buildYATKEfromSystematics (void) {
  
  std::fill_n(meanTKE,NUMA,187.0); // default is constant
  std::fill_n(sigTKE,NUMA,10.0);   // default is constant
  
  buildSystematicsTKEvsMass (); // Zeke Blaine, 4/9/2012
  
  for (int i=Amin; i<=Amax; i++) {
    if (YA[i]!=0.0) {
      for (int j=TKEmin; j<=TKEmax; j++) {
        YATKE[i][j]= YA[i]*exp(-pow(j-meanTKE[i],2)/(2*pow(sigTKE[i],2)));
      }
    }
  }
  
}

/*******************************************************************************
 * buildYATKEfromSystematics << 1.0.6 >>
 *------------------------------------------------------------------------------
 * Build the yields Y[A][TKE] using systematics for <TKE>(A) and <sigma_TKE>(A).
 * It follows the form: Y(A,TKE) = Y(A)xP(TKE|A) where P(TKE) is a Gaussian
 * function whose parameters are mass-dependent.
 * I. Stetcu: this systematics use data from thermal, scalled to the correct
 * incident energy
 ******************************************************************************/
void FissionFragments::buildYATKEfromSystematics (Yields *ffy) {
  
  TKEmin=999; TKEmax=0;
  maxYieldATKE=0.0;
  double s=0.;
  double yieldTKE[NUMTKE];
  for (int i=Amin; i<=Amax; i++) {
    if (YA[i]>0.0) {
      std::fill_n(yieldTKE,NUMTKE,0.); // default is constant
      for (int j=0; j<NUMTKE; j++) {
        yieldTKE[j]= ffy->yieldTKE(i, j );
        if( yieldTKE[j]>0.0){
          if(TKEmin>j) TKEmin=j ;
          if(TKEmax<j) TKEmax=j ;
        }
      }
      double ss=0.;
      for(int j=0; j<NUMTKE;j++)
        ss+=yieldTKE[j];
      for(int j=0; j<NUMTKE;j++){
        YATKE[i][j]=YA[i]*yieldTKE[j]/ss;
        s+=YATKE[i][j];
      }
    }
  }
  
  for(int i=Amin; i<=Amax; i++){
    for(int j=TKEmin; j<=TKEmax;j++){
      YATKE[i][j]/=s;
      if( maxYieldATKE<YATKE[i][j] ) maxYieldATKE = YATKE[i][j];
    }
  }
  
  cout << "maxYieldATKE=" << maxYieldATKE << endl;
  
}



/*******************************************************************************
 * buildWahlYA
 *------------------------------------------------------------------------------
 * Builds the mass yields Y(A) from Wahl's systematics.
 *------------------------------------------------------------------------------
 * Reference:
 * A.C.Wahl, "Systematics of Fission-Product Yields", Los Alamos Technical
 * Report, LA-13928 (2002).
 *------------------------------------------------------------------------------
 * From Zeke Blain, blaine2@rpi.edu, F95 original version (03/28/2012).
 ******************************************************************************/
void FissionFragments::buildWahlYA (void) {
  
  // Initialize Variables
  // double del7[4], Y1, Y3, Y5, D[7], del2, del4; // Ask Zeke
  double T, PE, NTtot, A4tot, Abar;
  double Y[7], sigma[7], delta[7];
  double sig1P[3], del5P[3], NTP[3];
  
  double Yield2[NUMA];
  
  int stop;
  int D1, F, G, H, J, K, M;
  double ECOR, NT1, Ip, Cp, Sp, Id, Cd, Sc, Dh, Sh, SYM, TH, Fz, Fn, nus1, Rf, Rj, sum, sum2;
  
  double Nus[NUMA], Nuh[NUMA], NuL[NUMA], R[NUMA], Nut[NUMA];
  
  // allocate (yields%Ya(Amin:Amax))
  // yields%Ya = 0.0_rk
  
  // Input parameters for calculating FN,FZ,sigz,delz
  double Sn = B1n[Zc][Ac];
  PE = incidentEnergy + Sn;
  
  double sig1[6] = {2.808, 8.685, -0.0454, 0.372, -0.620, -0.0122};
  double sig6[4] = {3.17, 0.0, 0.303, 0.0};
  double del5[6] = {25.34, 18.55, -0.0402, -1.220, -1.732, 0.0};
  double A4[4]   = {136.66, -0.177, 0.060, -0.038};
  double Y2[4]   = {43.00, -1.91, -3.41, 0.0};
  double NT[6]   = {1.563, 16.66, -0.00804, 0.0918, 0.0, 0.0};
  
  // Calculate Values for nu
  double En = incidentEnergy;
  Ip=20.0;
  if (Zt==98 && Ac==252) Ip=0;
  Cp=72.17;
  Id=-5.0;
  Cd=18.17;
  //     D1=78
  //     F=104
  //     G=107
  //     H=117
  //     J=127
  //     K=130
  //     M=156
  Fz=pow(-1.0,Zt);
  Fn=pow(-1.0, Ac-Zt);
  TH=11.47-0.166*Zt*Zt/Ac+0.093*(2-Fz-Fn)-Sn;
  if (Zt==98 && Ac==252) {
    NT1=2.286+0.147*(Zt-92)+0.054*(Ac-236)+0.040*(2-Fz-Fn)+(0.145-0.0043*(Ac-236));
  } else {
    NT1=2.286+0.147*(Zt-92)+0.054*(Ac-236)+0.040*(2-Fz-Fn)+(0.145-0.0043*(Ac-236))*(En-TH);
  }
  ECOR=exp(-0.05*PE);
  if ((NT1+1)>4) {
    nus1=NT1+1;
  } else {
    nus1=4.0;
  }
  Sp=(Ac-nus1)/2;
  Dh=Id/sqrt(3.14159*Cd);
  SYM=(Ac-nus1)/2.0;
  Sh=Id*(exp(-1*pow(130.0-SYM,2)/Cd))/sqrt(3.14159*Cd);
  Sc=Ac-130.0-(NT1-Dh-Sh);
  Rf=0.9-0.0075*PE;
  Rj=0.1+0.0075*PE;
  for (int i=Amin; i<=Amax; i++) {
    Nus[i] = Ip*exp(-1*pow(i-Sp,2)/Cp)/sqrt(3.14159*Cp);
    Nuh[i] = Id*exp(-1*pow(i-130.,2)/Cd)/sqrt(3.14159*Cd);
    NuL[i] = Id*exp(-1*pow(i-Sc,2)/Cd)/sqrt(3.14159*Cd);
    Nut[i] = (NT1+ECOR*(Nus[i]+Nuh[i]+NuL[i]));
  }
  D1 = (int) (Ac-156-Nut[156]);
  F = (int) (Ac-135-Nut[135]);
  G = (int) (Ac-130-Nut[130]);
  H = (int) ((128+(Ac-128-Nut[128]))/2);
  J = 130;
  K = 135;
  M = 156;
  for (int i=Amin; i<=Amax; i++) {
    if (i<D1) {
      R[i]=0.2;
    } else if (i<=F && i>=D1) {
      R[i]=0.2+(i-D1)*(Rf-0.2)/(F-D1);
    } else if (i<G && i>F) {
      R[i]=Rf;
    } else if (i<=H && i>=G) {
      R[i]=0.5+(H-i)*(Rf-0.5)/(H-G);
    } else if (i<J && i>H) {
      R[i]=0.5-(i-H)*(0.5-Rj)/(J-H);
    } else if (i<=K && i>=J) {
      R[i]=Rj;
    } else if (i<M && i>K) {
      R[i]=Rj+(i-K)*(0.8-Rj)/(M-K);
    } else {
      R[i]=0.8;
    }
    Nut[i] = R[i]*(NT1+ECOR*(Nus[i]+Nuh[i]+NuL[i]));
  }
  
  // Calculate input parameter for any given fissioning system
  for (int i=0; i<3; i++) {
    sig1P[i] = sig1[i]+sig1[i+3]*(Zt-92);
    del5P[i] = del5[i]+del5[i+3]*(Zt-92);
    NTP[i]   = NT[i]+NT[i+3]*(Zt-92);
  }
  sigma[0] = sig1P[0]+(sig1P[1]-sig1P[0])*(1-exp(sig1P[2]*PE));
  sigma[1] = 2.45;
  sigma[2] = 8.6;
  if (sigma[0]>sigma[2]) sigma[2]=sigma[0];
  sigma[3] = sigma[1];
  sigma[4] = sigma[0];
  sigma[5] = sig6[0]+sig6[2]*(Zt-92);
  sigma[6] = sigma[5];
  NTtot = NTP[0]+(NTP[1]-NTP[0])*(1-exp(NTP[2]*PE));
  A4tot = (A4[0]+A4[2]*(Zt-92)+(A4[1]+A4[3]*(Zt-92))*PE);
  A4tot = A4tot-Nut[int(A4tot)];
  delta[4] = del5P[0]+(del5P[1]-del5P[0])*(1-exp(del5P[2]*PE));
  delta[0] = -delta[4];
  delta[3] = A4tot-(Ac-NT1)/2;
  delta[1] = -delta[3];
  delta[2] = 0;
  delta[6] = 30.31;
  delta[5] = -delta[6];
  Abar = (Ac-NT1)/2;
  if (PE>8.0 && PE<=20.0) {
    Y[5] = 6.8-(6.8/12.0)*(PE-8.0);
  } else if (PE<8.0) {
    Y[5] = 6.8;
  } else if (Zt==93) {
    Y[5] = Y[5]/2.0;
  } else if (Zt<93 || PE>20.0) {
    Y[5] = 0;
  }
  Y[6] = Y[5];
  Y[1] = Y2[0]+Y2[2]*(Zt-92)+(Y2[1]+Y2[3]*(Zt-92))*PE;
  if (Y[1]<0) Y[1]=0;
  Y[3] = Y[1];
  if (PE<11.96) {
    Y[2] = 4.060*exp(0.470*(PE-11.96));
  } else {
    T = -0.030+0.0050*(Ac-236);
    Y[2] = 4.060+86.02*(1.0-exp(T*(PE-11.96)));
  }
  Y[0]=0.0;
  Y[4]=Y[0];
  stop=0;
  sum2=0;
  while (stop==0) {
    sum=0;
    for (int i=Amin; i<=Amax; i++) {
      Yield2[i]=0;
      for (int j=0; j<7; j++) {
        Yield2[i] = Yield2[i]+Y[j]*exp(-(pow(i-Nut[i]-Abar+delta[j],2))/(2*pow(sigma[j],2)))/(sigma[j]*sqrt(2*3.14159));
      }
    }
    for (int i=Amin; i<Amax; i++) {
      sum = sum+min(Yield2[i],Yield2[i+1])+abs(Yield2[i]-Yield2[i+1])/2.0;
    }
    if (sum>=199.5 && sum<=200.5) {
      stop=1;
    } else {
      Y[0]=Y[0]+0.1;
      Y[4]=Y[0];
    }
  }
  for (int i=Amin; i<=Amax; i++) {
    Yield2[i]=0;
    for (int j=0; j<7; j++) {
      Yield2[i] = Yield2[i]+Y[j]*exp(-(pow(i-Nut[i]-Abar+delta[j],2))/(2*pow(sigma[j],2)))/(sigma[j]*sqrt(2*3.14159));
    }
  }
  
  for (int i=Amin; i<=Amax; i++) {
    YA[i] = Yield2[i];
    sum2 += Yield2[i];
  }
  for (int i=Amin; i<=Amax; i++) {
    YA[i] /= sum2;
  }
  
}

/*------------------------------------------------------------------------------
 *                   SUBROUTINE buildSystematicsTKEvsMass
 *------------------------------------------------------------------------------
 * based on work by Ezekiel Blaine, RPI, blaine2@rpi.edu, 4/9/2012.
 * Using Hambsch data for U-235, Pu-239 and Cf-252.
 *----------------------------------------------------------------------------*/
void FissionFragments::buildSystematicsTKEvsMass (void) {
  
  int i;
  double A1, A2, A3, A4, A5, P1;
  
  P1=1.00+0.0112*(Zc-94)+0.0009*(Ac-240)-0.00135*pow(float(Zc-94),2);
  A3=15.20048736;
  A4=-0.0798617477;
  A5=0.0001558700726;
  A1=39639.39298*(1+0.000045*(Zc-94)+0.000022*(Ac-240)-0.0000018*pow(float(Ac-240),2));
  A2=-1272.328785*(1+0.0000013*(Zc-94)+0.0000006*(Ac-240)-0.000000055*pow(float(Ac-240),2));
  
  meanTKE[Asym]=A1+A2*Asym+A3*pow(float(Asym),2)+A4*pow(float(Asym),3)+A5*pow(float(Asym),4);
  
  for (i=1; i<=15; i++) {
    meanTKE[Asym+i] = A1+A2*(Asym+i)+A3*pow(float(Asym+i),2)+A4*pow(float(Asym+i),3)+A5*pow(float(Asym+i),4);
    meanTKE[Asym-i] = meanTKE[Asym+i];
  }
  
  for (i=21; i<=Amax-Asym; i++) {
    meanTKE[Asym+i] = 324.5*P1-1.028*(Asym+i);
    meanTKE[Asym-i] = meanTKE[Asym+i];
    if (Asym-i<0) { cerr << "[buildSystematicsTKEvsMass] error; negative array index!\n"; exit(-1); }
  }
  
  for (i=16; i<=20; i++) {
    meanTKE[Asym+i] = (meanTKE[Asym+21]-meanTKE[Asym+15])/(21-15)*(i-15)+meanTKE[Asym+15];
    meanTKE[Asym-i] = meanTKE[Asym+i];
  }
  
}

/*******************************************************************************
 * checkDistributions
 *------------------------------------------------------------------------------
 * Check the distributions of fission fragments produced by the method
 * 'generateInitialFissionFragmentHistories'.
 ******************************************************************************/
void FissionFragments::checkDistributions (string inputFilename, string outputFilename) {
  
  ifstream inputFile;
  ofstream outputFile;
  
  //  const int NUME = 501; // number of points on energy grid
  const int numberEnergies = 401; // number of points on energy grid
  
  int i;
  
  int Z, A, iTKE, Al=0, Ah, Zl=0, Zh;
  double KE, U,spin;
  int parity;
  
  double Ul=0.0, Uh=0.0;
  double KEl=0.0, KEh=0.0, TKE;
  double Q;
  
  double massYields[NUMA];
  double chargeYields[NUMZ];
  double tkeYields[NUMTKE];
  
  double excitationEnergies[NUMA];
  double spins[NUMA];
  
  double TKEvsMass[NUMA];
  double TXEvsMass[NUMA];
  
  double initialEnergiesLF[numberEnergies];
  double initialEnergiesHF[numberEnergies];
  double initialTXE[numberEnergies];
  
  double **YieldsATKE;
  YieldsATKE = new double *[NUMA];
  for (int i=0; i<NUMA; i++) {
    YieldsATKE[i] = new double [NUMTKE];
    for (int j=0; j<NUMTKE; j++) {
      YieldsATKE[i][j] = 0.0;
    }
  }
  
  int numberEventsPerMass[NUMA];
  int numberEventsPerCharge[NUMZ];
  int numberEventsPerTKE[NUMTKE];
  
  int totalNumberEvents;
  bool isSymmetric;
  
  isSymmetric=false;
  
  for (int i=0; i<NUMA; i++) {
    excitationEnergies[i]=0.0;
    spins[i]=0.0;
    TKEvsMass[i]=0.0;
    TXEvsMass[i]=0.0;
    numberEventsPerMass[i]=0;
    massYields[i]=0.0;
  }
  
  std::fill_n(chargeYields,NUMZ,0);
  std::fill_n(tkeYields,NUMTKE,0);
  
  std::fill_n(numberEventsPerCharge,NUMZ,0);
  std::fill_n(numberEventsPerTKE,NUMTKE,0);
  
  std::fill_n(initialEnergiesLF,numberEnergies,0.0);
  std::fill_n(initialEnergiesHF,numberEnergies,0.0);
  std::fill_n(initialTXE,numberEnergies,0.0);
  
  double **ZAYields;
  ZAYields = new double * [NUMZ];
  for (int i=0; i<NUMZ; i++) {
    ZAYields[i] = new double [NUMA];
    std::fill_n (ZAYields[i], NUMA, 0.0);
  }
  
  // define energy grid
  double energyGrid[numberEnergies];
  double deltaE = 0.25;
  for (i=0; i<numberEnergies; i++) {
    energyGrid[i] = i*deltaE;
  }
  
  totalNumberEvents=0;
  inputFilename = WORKDIR + inputFilename;
  inputFile.open(&inputFilename[0]);
  if (!inputFile) { cerr << "[checkDistributions] ERROR: cannot find history file: " << inputFilename << endl; exit(-1); }
  
  double mass;
  while (!inputFile.eof()) {
    
    inputFile >> Z >> A >> KE >> U >> spin >> parity; // >> mass;
    
    if (A<Asym) {
      Ul=U;
    } else {
      Uh=U;
    }
    
    excitationEnergies[A]+=U;
    spins[A]+=spin;
    
    i = int(U/deltaE)+1;
    if (i>numberEnergies) {
      cout << "U = " << U << "\n";
      cerr << "[checkDistributions] ERROR: should increase numberEnergies!" << endl;
      exit(-1);
    }
    
    if (A<Asym) {
      initialEnergiesLF[i]++;
    } else {
      initialEnergiesHF[i]++;
    }
    
    if (A>=Asym) {
      i = int((Ul+Uh)/deltaE)+1;
      initialTXE[i]++;
      TXEvsMass[A] += (Ul+Uh);
    }
    
    numberEventsPerMass[A]++;
    numberEventsPerCharge[Z]++;
    
    massYields[A]++;
    chargeYields[Z]++;
    totalNumberEvents++;
    
    ZAYields[Z][A]++;
    
    
    if (A<=Asym && !isSymmetric) {
      KEl=KE;
      Al=A;
      Zl=Z;
      if (A==Asym) isSymmetric=true;
    } else {
      KEh=KE;
      Ah=A;
      Zh=Z;
      TKE=KEl+KEh;
      iTKE= (int) floor(TKE+0.5);
      numberEventsPerTKE[iTKE]++;
      
      tkeYields[iTKE]++;
      TKEvsMass[Al]+= TKE;
      TKEvsMass[Ah]+= TKE;
      
      YieldsATKE[Al][iTKE]++;
      YieldsATKE[Ah][iTKE]++;
      
      // Q-values
      Q = (mass_excess(Zc,Ac)-mass_excess(Zl,Al)-mass_excess(Zh,Ah)+incidentEnergy+SnCompound);
      QfA[Al] += Q;
      QfA[Ah] += Q;
      QfTKE[iTKE] += Q;
      
      if (A==Asym) isSymmetric=false; // reset
      
    }
    
    if (A>=Asym) { Ul=0.0; } else { Uh=0.0; } // reset
    
  }
  inputFile.close();
  
  // 2 fragments (2 lines) <=> 1 event
  totalNumberEvents /= 2;
  
  //---------------------------------------------------------------------------------------
  outputFilename = WORKDIR + outputFilename;
  outputFile.open(&outputFilename[0]);
  
  int gnuplotIndex=-1;
  
  // Mass Yields Y(A)
  outputFile << "# [gnuplot " << ++gnuplotIndex << "] Mass Yields Y(A)" << endl << endl;
  for (i=0; i<NUMA; i++) {
    if (numberEventsPerMass[i]!=0) {
      massYields[i] /= totalNumberEvents;
      outputFile << i << " " << massYields[i] << endl;
    }
  }
  outputFile << endl;
  
  // Charge Yields Y(Z)
  outputFile << "# [gnuplot " << ++gnuplotIndex << "] Charge Yields Y(Z)" << endl << endl;
  for (i=0; i<NUMZ; i++) {
    if (numberEventsPerCharge[i]!=0) {
      chargeYields[i] /= totalNumberEvents;
      outputFile << i << " " << chargeYields[i] << endl;
    }
  }
  outputFile << endl;
  
  // TKE Yields Y(TKE)
  outputFile << "# [gnuplot " << ++gnuplotIndex << "] TKE   Y(TKE)  <Qf>(TKE)" << endl << endl;
  for (i=0; i<NUMTKE; i++) {
    if (numberEventsPerTKE[i]!=0) {
      tkeYields[i] /= totalNumberEvents;
      QfTKE[i] /= numberEventsPerTKE[i];
      outputFile << i << " " << tkeYields[i] << " " << QfTKE[i] << endl;
    }
  }
  outputFile << endl;
  
  // P(U)
  outputFile << "# [gnuplot " << ++gnuplotIndex << "] Initial Excitation Energy Probability P(Ui)" << endl << endl;
  double sum=0, sumLF=0, sumHF=0;
  for (i=0; i<numberEnergies; i++) {
    sum   += initialTXE[i];
    sumLF += initialEnergiesLF[i];
    sumHF += initialEnergiesHF[i];
  }
  sum   *= deltaE;
  sumLF *= deltaE;
  sumHF *= deltaE;
  for (i=0; i<numberEnergies; i++) {
    initialTXE[i]        /= sum;
    initialEnergiesLF[i] /= sumLF;
    initialEnergiesHF[i] /= sumHF;
    outputFile << energyGrid[i] << " " << initialEnergiesLF[i] << " " <<
    initialEnergiesHF[i] << " " << initialTXE[i] << endl;
  }
  outputFile << endl;
  
  // <U>=f(A)
  outputFile << "# [gnuplot " << ++gnuplotIndex << "] <U>=f(A) " << endl << endl;
  for (i=0; i<NUMA; i++) {
    if (numberEventsPerMass[i]!=0) { excitationEnergies[i] /= numberEventsPerMass[i]; }
    outputFile << i << " " << excitationEnergies[i] << endl;
  }
  outputFile << endl;
  
  // <spin>=f(A)
  outputFile << "# [gnuplot " << ++gnuplotIndex << "] <spin>=f(A) " << endl << endl;
  for (i=0; i<NUMA; i++) {
    if (numberEventsPerMass[i]!=0) {
      spins[i] /= numberEventsPerMass[i];
      outputFile << i << " " << spins[i] << endl;
    }
  }
  outputFile << endl;
  
  // P(TXE)
  outputFile << "# [gnuplot " << ++gnuplotIndex << "] Total Excitation Energy Probability P(TXE)" << endl << endl;
  sum=0;
  for (i=0; i<=numberEnergies; i++) {
    sum += initialTXE[i];
  }
  for (i=0; i<=numberEnergies; i++) {
    initialTXE[i] /= sum;
    outputFile << energyGrid[i] << " " << initialTXE[i] << endl;
  }
  outputFile << endl;
  
  // <TKE>(A)
  outputFile << "# [gnuplot " << ++gnuplotIndex << "] A   <TKE>(A) <Q>(A)" << endl << endl;
  for (i=0; i<NUMA; i++) {
    if (numberEventsPerMass[i]!=0) {
      TKEvsMass[i] /= numberEventsPerMass[i];
      QfA[i] /= numberEventsPerMass[i];
      outputFile << i << " " << TKEvsMass[i] << " " << QfA[i] << " " << numberEventsPerMass[i] << endl;
    }
  }
  outputFile << endl;
  
  // <TXE>(A)
  outputFile << "# [gnuplot " << ++gnuplotIndex << "] <TXE>=f(A) " << endl << endl;
  for (i=0; i<NUMA; i++) {
    if (numberEventsPerMass[i]!=0) {
      TXEvsMass[i] /= numberEventsPerMass[i];
      outputFile << i << " " << TXEvsMass[i] << endl;
    }
  }
  outputFile << endl;
  
  outputFile.close();
  
}

/*******************************************************************************
 * computeNeutronSeparationEnergies
 *------------------------------------------------------------------------------
 * Compute one and two-neutron separation energies from Audi table of mass
 * excess.
 ******************************************************************************/
void FissionFragments::computeNeutronSeparationEnergies (void) {
  
  double Mn = mass_excess(0,1); // neutron mass excess (MeV)
  
  for (int i=Amin; i<=Amax; i++) {
    for (int j=Zmin; j<=Zmax; j++) {
      if (YZA2[j][i]!=0.0) {
        S1n[j][i] = -mass_excess(j,i)+mass_excess(j,i-1)+Mn;
        S2n[j][i] = -mass_excess(j,i)+mass_excess(j,i-2)+2*Mn;
      }
    }
  }
  
}


/*------------------------------------------------------------------------------
 *                       SUBROUTINE computePfactors
 *------------------------------------------------------------------------------
 * Computes Pfactors for every fission fragment yield, and returns a mass-
 * dependent array Pf(Amin:Ac/2) that contains the Y(Z|A)-weighted Pfactors.
 *----------------------------------------------------------------------------*/
void FissionFragments::computePfactors (void) {
  int Zpl, Zph;
  for (int i=Amin; i<=Asym; i++) {
    Zpl = (int) (Z0[i]+0.5);
    Zph = Zc-Zpl;
    for (int j=-dZ; j<=dZ; j++) {
      Pfactors[i][Zpl+j]    = Pfactor(Zpl+j, i-Zpl-j);
      Pfactors[Ac-i][Zph-j] = Pfactor(Zph-j, Ac-i-Zph+j);
    }
  }
  
  ofstream out;
  out.open("Pfactors.out");
  double pfl, pfh, sum;
  double y;
  for (int i=Amin; i<=Asym; i++) {
    Zpl = (int) (Z0[i]+0.5);
    Zph = Zc-Zpl;
    pfl = 0.0;
    pfh = 0.0;
    sum = 0.0;
    for (int j=-dZ; j<=dZ; j++) {
      y = YZA2[Zpl+j][i];
      pfl = pfl + y * Pfactors[i][Zpl+j];
      pfh = pfh + y * Pfactors[Ac-i][Zph-j];
      sum += y;
    }
    if (sum!=0.0) { pfl /= sum; pfh /= sum; }
    
    //    pfl = Pfactors[i][Zpl];
    //    pfh = Pfactors[Ac-i][Zc-Zpl];
    
  }
  out.close();
  
}

/*******************************************************************************
 * computeTemperatures
 *------------------------------------------------------------------------------
 * Computes temperatures for a given nucleus (Z,A) on a given energy grid
 * (with NUME number of energy points).
 * The temperature is calculated as: 1/T = ( d/dU ) log( LevelDensity( U )
 ******************************************************************************/
void FissionFragments::computeTemperatures (int Z, int A) {
  
  //	cout << "# Computing Temperatures for " << Z << " " << A << "\n\n";
  
  Nucleus nucleus;
  
  nucleus.za.setZA(Z,A);
  nucleus.ndisc = riplReadDiscreteLevels(&nucleus.za,nucleus.lev,reassign);
  statFixDiscreteLevels(&nucleus);
  statSetupEnergyBin(&nucleus);
  statSetupLevelDensityParameter(&nucleus,&nucleus.ldp);
  
  double *ld;
  ld = new double [NUME];
  
  for (int k=0; k<NUME; k++) {
    ld[k] = ldLevelDensity(ldEnergyGrid[k], A, &nucleus.ldp);
  }
  
  double de;
  double *temp;
  temp = new double [NUME];
  
  double *ldp  = new double [NUME];
  
  nucleus.za.setZA(Z,A);
  nucleus.ndisc = riplReadDiscreteLevels(&nucleus.za,nucleus.lev,reassign);
  statFixDiscreteLevels(&nucleus);
  statSetupEnergyBin(&nucleus);
  statSetupLevelDensityParameter(&nucleus,&nucleus.ldp);
  
  for (int k=0; k<NUME; k++) {
    ldp[k] = ldDensityParameter(ldEnergyGrid[k], A, &nucleus.ldp);
  }
  
  for (int k=2; k<NUME-2; k++) {
    de = ldEnergyGrid[k] - ldEnergyGrid[k-1];
    temp[k] = ( -log(ld[k+2]) + 8.0*log(ld[k+1]) - 8.0*log(ld[k-1]) + log(ld[k-2]) ) / (12.*de);
  }
  
  de = ldEnergyGrid[1]-ldEnergyGrid[0];
  temp[0] = ( -25.0*log(ld[0]) + 48.0*log(ld[1]) - 36.0*log(ld[2]) + 16.0*log(ld[3]) - 3.0*log(ld[4]) ) / (12.0*de);
  
  de = ldEnergyGrid[2]-ldEnergyGrid[1];
  temp[1] = ( -25.0*log(ld[1]) + 48.0*log(ld[2]) - 36.0*log(ld[3]) + 16.0*log(ld[4]) - 3.0*log(ld[5]) ) / (12.0*de);
  
  de = ldEnergyGrid[NUME-2]-ldEnergyGrid[NUME-3];
  temp[NUME-2] = ( 25.0*log(ld[NUME-2]) - 48.0*log(ld[NUME-3]) + 36.0*log(ld[NUME-4]) - 16.0*log(ld[NUME-5]) + 3.0*log(ld[NUME-6]) ) / (12.0*de);
  
  de = ldEnergyGrid[NUME-1]-ldEnergyGrid[NUME-2];
  temp[NUME-1] = ( 25.0*log(ld[NUME-1]) - 48.0*log(ld[NUME-2]) + 36.0*log(ld[NUME-3]) - 16.0*log(ld[NUME-4]) + 3.0*log(ld[NUME-5]) ) / (12.0*de);
  
  for (int k=0; k<NUME; k++) temp[k] = 1.0/temp[k];
  
  /*  for (int k=0; k<NUME; k++) {
   cout << k << " " << ldEnergyGrid[k] << " " << ld[k] << " " << temp[k] << "\n";
   }
   */
  
  cout << 1000*Z+A << " ";
  
  /*	for (int k=0; k<NUME; k++) {
   cout << fixed << setprecision(2) << temp[k] << " ";
   }
   cout << endl;
   */
  
  for (int k=0; k<41; k+=4) { // 0-10 MeV per step of 1.0 MeV
    cout << fixed << setprecision(2) << ldp[k] << " ";
  }
  for (int k=48; k<81; k+=8) { // 11-20 MeV per step of 2.0 MeV
    cout << fixed << setprecision(2) << ldp[k] << " ";
  }
  for (int k=100; k<NUME; k+=20) { // 25-100 MeV per step of 5.0 MeV
    cout << fixed << setprecision(2) << ldp[k] << " ";
  }
  
  cout << endl;
  
  
  delete ld;
  delete temp;
  
}


/*******************************************************************************
 * computeTemperatureTables
 *------------------------------------------------------------------------------
 * Computes temperature tables for (A,Z) on a given energy grid (with NUME
 * number of energy points), which are later used to compute the spin
 * parameter. 1/T = ( d/dU ) log( LevelDensity( U )
 ******************************************************************************/
void FissionFragments::computeTemperatureTables () {
  
  cout << "[computeTemperatureTables] Running... ";
  
  double de;
  
  temperatures = new double **[NUMA];
  for (int i=Amin; i<=Amax; i++) {
    temperatures[i] = new double *[NUMdZ];
    for (int j=0; j<=2*dZ; j++) {
      temperatures[i][j] = new double [NUME];
      //if (YZA[j][i]!=0.0) { << 1.0.6 >>
      
      for (int k=2; k<NUME-2; k++) {
        de = ldEnergyGrid[k] - ldEnergyGrid[k-1];
        temperatures[i][j][k] =
        ( -log(levelDensities[i][j][k+2]) + 8.0*log(levelDensities[i][j][k+1]) -
         8.0*log(levelDensities[i][j][k-1]) + log(levelDensities[i][j][k-2]) ) / (12.*de);
        //         cout << " ... " << k << "  " << temperatures[i][j][k] << "\n";
      }
      
      de = ldEnergyGrid[1]-ldEnergyGrid[0];
      temperatures[i][j][0] =
      ( -25.0*log(levelDensities[i][j][0]) + 48.0*log(levelDensities[i][j][1]) -
       36.0*log(levelDensities[i][j][2]) + 16.0*log(levelDensities[i][j][3]) -
       3.0*log(levelDensities[i][j][4]) ) / (12.0*de);
      
      de = ldEnergyGrid[2]-ldEnergyGrid[1];
      temperatures[i][j][1] =
      ( -25.0*log(levelDensities[i][j][1]) + 48.0*log(levelDensities[i][j][2]) -
       36.0*log(levelDensities[i][j][3]) + 16.0*log(levelDensities[i][j][4]) -
       3.0*log(levelDensities[i][j][5]) ) / (12.0*de);
      
      de = ldEnergyGrid[NUME-2]-ldEnergyGrid[NUME-3];
      temperatures[i][j][NUME-2] =
      ( 25.0*log(levelDensities[i][j][NUME-2]) - 48.0*log(levelDensities[i][j][NUME-3]) +
       36.0*log(levelDensities[i][j][NUME-4]) - 16.0*log(levelDensities[i][j][NUME-5]) +
       3.0*log(levelDensities[i][j][NUME-6]) ) / (12.0*de);
      
      de = ldEnergyGrid[NUME-1]-ldEnergyGrid[NUME-2];
      temperatures[i][j][NUME-1] =
      ( 25.0*log(levelDensities[i][j][NUME-1]) - 48.0*log(levelDensities[i][j][NUME-2]) +
       36.0*log(levelDensities[i][j][NUME-3]) - 16.0*log(levelDensities[i][j][NUME-4]) +
       3.0*log(levelDensities[i][j][NUME-5]) ) / (12.0*de);
      
      for (int k=0; k<NUME; k++) {
        temperatures[i][j][k] = 1.0/temperatures[i][j][k];
      }
      
      //}
    }
  }
  
  cout << "OK\n";
  
}

/*******************************************************************************
 * computeWahlParametersForChargeDistribution
 *------------------------------------------------------------------------------
 * Calculate the Wahl Parameters {Fn(A), Fz(A), delz(A), sigz(A)} required for
 * computing the charge distribution Y(Z|A) for any given fissioning isotope.
 *------------------------------------------------------------------------------
 * Reference:
 * A.C.Wahl, "Systematics of Fission-Product Yields", Los Alamos Technical
 * Report, LA-13928 (2002).
 *------------------------------------------------------------------------------
 * From Zeke Blain, blaine2@rpi.edu, F95 original version (March 2012).
 ******************************************************************************/
void FissionFragments::computeWahlParametersForChargeDistribution (void) {
  
  double sigz140t, delz140t, Fz140t, Fn140t, sigzSLt, delzSLt, FzSLt, SL50t,
  sigz50t, delzmaxt, sigzSLWt, delzSLWt, FzSLWt, FnSLWt;
  
  double delz[NUMA], sigz[NUMA], Fz[NUMA], Fn[NUMA];
  
  double PE = incidentEnergy + B1n[Zc][Ac];
  
  double sigz140[5] = {0.566, 0.0, 0.0064, 0.0109, 0.0};
  double delz140[5] = {-0.487, 0.0, 0.0180, 0.0, -0.00203};
  //	double Fz140[5]   = {1.207, 0.0, -0.0420, 0.0, 0.0022};
  double Fz140[5]   = {1.242, 0.0, -0.0183, -0.0152, 0.0};
  double Fn140[5]   = {1.076, 0.0, 0.0, 0.0, 0.0};
  double sigzSL[5]  = {-0.0038, 0.0, 0.0, 0.0, 0.0};
  double delzSL[5]  = {-0.0080, 0.0, 0.0, 0.0, 0.0};
  double FzSL[5]    = {0.0030, 0.0, 0.0, 0.0, 0.0};
  double SL50[5]    = {0.191, 0.0, -0.0076, 0.0, 0.0};
  double sigz50[5]  = {0.356, 0.060, 0.0, 0.0, 0.0};
  double delzmax[5] = {0.699, 0.0, 0.0, 0.0, 0.0};
  double sigzSLW[5] = {-.045, 0.0094, 0.0, 0.0, 0.0};
  double delzSLW[5] = {0.0, -0.0045, 0.0, 0.0, 0.0};
  double FzSLW[5]   = {0.159, -0.028, 0.0, 0.0, 0.0};
  double FnSLW[5]   = {0.039, 0.0, 0.0, 0.0, 0.0};
  
  std::fill_n(Z0,NUMA,0.0);
  std::fill_n(sZ0,NUMA,0.0);
  std::fill_n(FZZ0,NUMA,0.0);
  std::fill_n(FNZ0,NUMA,0.0);
  
  //============================================================================
  // TEMPORARY... FOR ANYTHING BUT n+U-235 FISSION
  //============================================================================
  /*  if (Ac!=236) {
   
   for (int i=Amin; i<=Amax; i++) {
   delz[i]  = 0.5; // dZ
   if(i>Asym) delz[i] = -0.5; // dZ
   sigz[i] = sigmaZ; // width of distribution
   Fz[i] = 1.0; // FZ - odd-even factor
   Fn[i] = 1.0; // FN - odd-even factor
   }
   
   //-- save Wahl parameters
   double x = float(Zt)/Ac;
   for (int i=Amin; i<=Amax; i++) {
   Z0[i]   = (x*i+delz[i]);
   sZ0[i]  = sigz[i];
   FZZ0[i] = Fz[i];
   FNZ0[i] = Fn[i];
   //      cout << i << " " << Z0[i] << " " << sZ0[i] << "\n";
   }
   
   //-- Find true (Zmin, Zmax)
   Zmin = (int) floor(Z0[Amin]+0.5)-dZ;
   Zmax = (int) floor(Z0[Amax]+0.5)+dZ;
   
   return;
   
   }
   */
  
  
  // Calculate input parameter for any given fissioning system
  // >> SHOULD 6.551 BE REPLACED BY neutronSeparationEnergy? I THINK SO...
  // WHAT ABOUT (92,236)... IS IT JUST VALID FOR U235?
  if (PE <= 8.0) {
    
    sigz140t = sigz140[0]+sigz140[1]*(Zt-92)+sigz140[2]*(Ac-236)+sigz140[3]*(PE-6.551)+sigz140[4]*(Ac-236)*(Ac-236);
    delz140t = delz140[0]+delz140[1]*(Zt-92)+delz140[2]*(Ac-236)+delz140[3]*(PE-6.551)+delz140[4]*(Ac-236)*(Ac-236);
    Fz140t   = Fz140[0]+Fz140[1]*(Zt-92)+Fz140[2]*(Ac-236)+Fz140[3]*(PE-6.551)+Fz140[4]*(Ac-236)*(Ac-236);
    Fn140t   = Fn140[0]+Fn140[1]*(Zt-92)+Fn140[2]*(Ac-236)+Fn140[3]*(PE-6.551)+Fn140[4]*(Ac-236)*(Ac-236);
    sigzSLt  = sigzSL[0]+sigzSL[1]*(Zt-92)+sigzSL[2]*(Ac-236)+sigzSL[3]*(PE-6.551)+sigzSL[4]*(Ac-236)*(Ac-236);
    delzSLt  = delzSL[0]+delzSL[1]*(Zt-92)+delzSL[2]*(Ac-236)+delzSL[3]*(PE-6.551)+delzSL[4]*(Ac-236)*(Ac-236);
    FzSLt    = FzSL[0]+FzSL[1]*(Zt-92)+FzSL[2]*(Ac-236)+FzSL[3]*(PE-6.551)+FzSL[4]*(Ac-236)*(Ac-236);
    SL50t    = SL50[0]+SL50[1]*(Zt-92)+SL50[2]*(Ac-236)+SL50[3]*(PE-6.551)+SL50[4]*(Ac-236)*(Ac-236);
    sigz50t  = sigz50[0]+sigz50[1]*(Zt-92)+sigz50[2]*(Ac-236)+sigz50[3]*(PE-6.551)+sigz50[4]*(Ac-236)*(Ac-236);
    delzmaxt = delzmax[0]+delzmax[1]*(Zt-92)+delzmax[2]*(Ac-236)+delzmax[3]*(PE-6.551)+delzmax[4]*(Ac-236)*(Ac-236);
    sigzSLWt = sigzSLW[0]+sigzSLW[1]*(Zt-92)+sigzSLW[2]*(Ac-236)+sigzSLW[3]*(PE-6.551)+sigzSLW[4]*(Ac-236)*(Ac-236);
    delzSLWt = delzSLW[0]+delzSLW[1]*(Zt-92)+delzSLW[2]*(Ac-236)+delzSLW[3]*(PE-6.551)+delzSLW[4]*(Ac-236)*(Ac-236);
    FzSLWt   = FzSLW[0]+FzSLW[1]*(Zt-92)+FzSLW[2]*(Ac-236)+FzSLW[3]*(PE-6.551)+FzSLW[4]*(Ac-236)*(Ac-236);
    FnSLWt   = FnSLW[0]+FnSLW[1]*(Zt-92)+FnSLW[2]*(Ac-236)+FnSLW[3]*(PE-6.551)+FnSLW[4]*(Ac-236)*(Ac-236);
    
  } else if (PE<=20.0) {
    
    sigz140[0] = 0.542; sigz140[1] = 1.310; sigz140[2] = 0.033; sigz140[3] = 0.0; sigz140[4] = -0.005;
    delz140[0] = -0.428; delz140[1] = 0.0; delz140[2] = 0.0; delz140[3] = 0.164; delz140[4] = -0.0116;
    SL50[0] = 0.191; SL50[1] = 0.0; SL50[2] = -0.0076; SL50[3] = 0.0; SL50[4] = 0.0;
    sigz50[0] = 0.542; sigz50[1] = 1.310; sigz50[2] = 0.033; sigz50[3] = 0.0; sigz50[4] = -0.005;
    
    sigz140t = (sigz140[0]+sigz140[2]*(Zt-92))+((sigz140[1]+sigz140[3]*(Zt-92))-(sigz140[0]+sigz140[2]*(Zt-92)))*(1.0-exp(-sigz140[4]*PE));
    delz140t = (delz140[0]+delz140[2]*(Zt-92))+((delz140[1]+delz140[3]*(Zt-92))-(delz140[0]+delz140[2]*(Zt-92)))*(1.0-exp(-delz140[4]*PE));
    Fz140t   = 1.0;
    Fn140t   = 1.0;
    sigzSLt  = 0.0;
    delzSLt  = 0.0;
    FzSLt    = 0.0;
    SL50t    = (SL50[0]+SL50[2]*(Zt-92))+((SL50[1]+SL50[3]*(Zt-92))-(SL50[0]+SL50[2]*(Zt-92)))*(1.0-exp(-SL50[4]*PE));
    sigz50t  = (sigz50[0]+sigz50[2]*(Zt-92))+((sigz50[1]+sigz50[3]*(Zt-92))-(sigz50[0]+sigz50[2]*(Zt-92)))*(1.0-exp(-sigz50[4]*PE));
    delzmaxt = 0.0;
    sigzSLWt = 0.0;
    delzSLWt = 0.0;
    FzSLWt   = 0.0;
    FnSLWt   = 0.0;
    
    
  } else {
    
    cerr << "[computeWahlParameters] Wahl parameters only defined up to Einc+Sn=20 MeV!" << endl;
    exit(-1);
    
  }
  
  // Calculate the region values
  
  double F1=floor((250.-Ac)/14.+0.5);  // F1=anint((250.-Ac)/14.); in F95
  double F2=1.-F1;
  double AK1=50.0*float(Ac)/Zt-delzmaxt/SL50t;
  double AK2=(50.0-delzmaxt)*float(Ac)/Zt;
  double Apmax=F1*AK1+F2*AK2;
  
  int B1=70;
  int B2= (int) floor(77+0.036*(Ac-236)+0.5);      // nint() in F95
  int B4= (int) floor((delzmaxt-delz140t+Apmax*SL50t+140*delzSLt)/(SL50t+delzSLt)+0.5); // nint() in F95
  int B3=(Ac-B4);
  int B5=(Ac-B2);
  int B6=(Ac-B1);
  int Bb=int(Apmax+0.5); // nint() in F95
  int Ba=int(Ac-Apmax+0.5); // nint() in F95
  
  // Calculate values of sigz, delz, Fn, Fz for give Ap values
  for (int i=Amin; i<=Amax; i++) {
    // Peak Regions
    if ((i>=B2 && i<=B3) || (i>=B4 && i<=B5)) {
      if (i>Ac/2.0) {
        delz[i]=delz140t+delzSLt*(i-140.0);
        sigz[i]=sigz140t+sigzSLt*(i-140.0);
        Fz[i]=Fz140t+FzSLt*(i-140.0);
        Fn[i]=Fn140t+FzSLt*(i-140.0);
      } else if (i<Ac/2.0) {
        delz[i]=-1*delz140t+delzSLt*(i-(Ac-140));
        sigz[i]=sigz140t-sigzSLt*(i-(Ac-140));
        Fz[i]=Fz140t-FzSLt*(i-(Ac-140));
        Fn[i]=Fn140t-FzSLt*(i-(Ac-140));
      }
      // Fn(i)=Fn140t
    }
  }
  
  for (int i=Amin; i<=Amax; i++) {
    // Near Symmetry Region
    if (i>B3 && i<B4) {
      Fn[i]=1.;
      Fz[i]=1.;
      if (i>B3 && i<=Ba) {
        delz[i]=delz[B3]-SL50t*(i-B3);
        sigz[i]=sigz50t;
      } else if (i>Ba && i<Bb) {
        delz[i]=delz[Ba]+(i-Ba)*(2.0*delz[Ba])/(Ba-Bb);
        sigz[i]=sigz140t-sigzSLt*(140-Bb);
      } else if (i>=Bb && i<B4) {
        delz[i]=delz[B4]+SL50t*(B4-i);
        sigz[i]=sigz50t;
      }
    }
  }
  
  for (int i=Amin; i<=Amax; i++) {
    // Wing Regions
    if ((i>B1 && i<B2) || (i>B5 && i<B6)) {
      if (i>Ac/2.0) {
        delz[i]=delz[B5]-delzSLWt*(i-B5);
        sigz[i]=sigz[B5]+sigzSLWt*(i-B5);
        Fz[i]=Fz140t+FzSLWt*(i-B5);
        Fn[i]=Fn140t+FnSLWt*(i-B5);
      } else if (i<Ac/2.0) {
        delz[i]=delz[B2]+delzSLWt*(B2-i);
        sigz[i]=sigz[B5]+sigzSLWt*(B2-i);
        Fz[i]=Fz140t+FzSLWt*(B2-i);
        Fn[i]=Fn140t+FnSLWt*(B2-i);
      }
    }
    
    // Far Wing Regions
    if (i<=B1 || i>=B6) {
      if (i>Ac/2.0) {
        delz[i]=delz[B5];
        sigz[i]=sigz[B5];
        Fz[i]=Fz140t;
        Fn[i]=Fn140t;
      } else if (i<Ac/2.0) {
        delz[i]=delz[B2];
        sigz[i]=sigz[B5];
        Fz[i]=Fz140t;
        Fn[i]=Fn140t;
      }
    }
  }
  
  //-- save Wahl parameters
  double x = float(Zt)/Ac;
  for (int i=Amin; i<=Amax; i++) {
    Z0[i]   = (x*i+delz[i]);
    sZ0[i]  = sigz[i];
    FZZ0[i] = Fz[i];
    FNZ0[i] = Fn[i];
    //    FZZ0[i] = 1.0;
    //    FNZ0[i] = 1.0;
    //    cout << i << " " << Z0[i] << " " << sZ0[i] << "\n";
  }
  
  //-- Find true (Zmin, Zmax)
  Zmin = (int) floor(Z0[Amin]+0.5)-dZ;
  Zmax = (int) floor(Z0[Amax]+0.5)+dZ;
  
}

/*------------------------------------------------------------------------------
 *                             FUNCTION Pfactor
 *------------------------------------------------------------------------------
 * Computes the Pfactor for a nucleus of charge Z and neutron number N.
 *----------------------------------------------------------------------------*/
double FissionFragments::Pfactor (int Z, int N) {
  int Nn, Np; //-- valence neutrons and protons
  int i;
  double p;
  
  const int magicNumbers[3] = { 26, 50, 82 }; //!-- spherical closed shells
  
  Nn=99; Np=99;
  for (i=0; i<3; i++) {
    Nn = min(Nn,abs(N-magicNumbers[i]));
    Np = min(Np,abs(Z-magicNumbers[i]));
  }
  
  if (Np+Nn==0) {
    p=0.0;
  } else {
    p = Np*Nn/float(Np+Nn);
  }
  
  return (p);
}

/*******************************************************************************
 * readDeformations
 *------------------------------------------------------------------------------
 * Read ground-state deformations from FRDM95 calculations stored in the file
 * 'mass-frdm95.dat'.
 ******************************************************************************/
void FissionFragments::readDeformations (void) {
  
  string filename = DATADIR ;
  filename.append("/mass-frdm95.dat");
  
  ifstream inputFile (&filename[0]);
  if (!inputFile) {
    cerr << "[readDeformations] file not found, exiting program" << endl ;
    exit (-1) ;
  }
  
  int a, z;
  for (z=0 ; z<NUMZ; z++) {
    for (a=0; a<NUMA; a++) {
      beta2[z][a]=-100.;
    }
  }
  
  char space=' ';
  string line, str;
  do {
    getline(inputFile,line);
  } while (line.find("#")!=string::npos);
  
  do {
    str=line.substr(0,4); z=atoi(&str[0]); if (z>99) break;
    str=line.substr(4,4); a=atoi(&str[0]); if (a>199) continue;
    str=line.substr(44,7);
    if (str[6]!=space) beta2[z][a]=atof(&str[0]);
  } while (getline(inputFile,line));
  
  inputFile.close() ;
  
}

/*******************************************************************************
 * readNeutronBindingEnergies
 *------------------------------------------------------------------------------
 * read one-neutron and two-neutrons separation energies from Audi-Wapstra
 * tables. Ref. G.Audi et al., NPA 729 (2003) 337-676.
 * 'AUDI2003DATAFILE' is set in config.h.
 ******************************************************************************/
void FissionFragments::readNeutronBindingEnergies (void) {
  
  ifstream dataFile;
  string str;
  string string1, string2;
  
  string dataFilename = "~/Development/CGMF/Data/neutronSeparationEnergies.Audi2003.dat";
  
  int A, Z;
  char Symbol[2];
  double B1, B2;
  
  // to replace with BOOST arrays
  for (int i=0; i<NUMA; i++) {
    for (int j=0; j<NUMZ; j++) {
      B1n[j][i] = -99.0;
      B2n[j][i] = -99.0;
    }
  }
  
  dataFile.open(&dataFilename[0]);
  if (!dataFile) { cerr << "ERROR: Cannot find neutron separation energies from Audi 2003 file!" << endl; exit(-1); }
  
  while (!dataFile.eof()) {
    
    getline(dataFile, str);
    if (str.substr(0,1)!="#" && !str.empty()) {
      
      istringstream(str) >> A >> Symbol >> Z;
      
      string1 = str.substr(12,9);
      string2 = str.substr(32,9);
      
      if (string1.find("*")==string::npos) {
        istringstream(string1) >> B1;
        B1n[Z][A] = B1 * 1e-3; // keV to MeV
      }
      if (string2.find("*")==string::npos) {
        istringstream(string2) >> B2;
        B2n[Z][A] = B2 * 1e-3; // keV to MeV
      }
      
    }
    
    //    cout << Z << " " << A << " " << B1n[Z][A] << " " << B2n[Z][A] << endl;
    
  }
  
  
  dataFile.close();
  
}

/*******************************************************************************
 *                         subroutine readMasses
 *------------------------------------------------------------------------------
 * Read excess masses from Audi-Wapstra 2003 atomic mass evaluation.
 * Ref. G.Audi et al., NPA 729 (2003) 3-128.
 *
 * >> April 18, 2012: now reads directly from mass_excess(Z,A) from CGM.
 *
 ******************************************************************************/
void FissionFragments::readMasses (void) {
  
  ifstream dataFile;
  string str;
  int i,j,k;
  
  // to replace with BOOST arrays
  for (i=0; i<NUMA; i++) {
    for (j=Zmin; j<=Zmax; j++) {
      ExpMassMeV[j][i] = -99.0;
      ExpMassAMU[j][i] = -99.0;
    }
  }
  
  // replaced with CGM-retrieved values
  for (i=Amin; i<=Amax; i++) {
    for (k=-dZ; k<=dZ; k++) {
      j=int(Z0[i]+0.5)+k;
      ExpMassAMU[j][i] = float(i)+mass_excess(j,i)*1e6/amu;
      ExpMassMeV[j][i] = mass_excess(j,i);
    }
  }
  
  //-- special masses as stored in physics.h class
  /*  ExpMassAMU[0][1] = neutronMass;
   ExpMassAMU[1][1] = protonMass;
   ExpMassAMU[1][2] = deuteronMass;
   ExpMassAMU[1][3] = tritonMass;
   ExpMassAMU[2][3] = he3Mass;
   ExpMassAMU[2][4] = alphaMass;
   */
  
  massExcessCompound = mass_excess(Zc, Ac);
  
  SnCompound = -mass_excess(Zc,Ac)+mass_excess(Zc,Ac-1)+mass_excess(0,1);
  if (incidentEnergy==0.0) SnCompound=0.0; // SPONTANEOUS FISSION
  
  // no need for 'Duflo-Zuker' anymore... if not found in Audi table, use FRLDM95
  
}

/*******************************************************************************
 * readYieldsATKE
 *------------------------------------------------------------------------------
 * Reads fission fragment yields in the form of Y(A,TKE).
 ******************************************************************************/
void FissionFragments::readYieldsATKE (string yieldsFilename) {
  
  ifstream yieldsFile;
  int A, TKE;
  double yields;
  
  yieldsFile.open(&yieldsFilename[0]);
  if (!yieldsFile) { cerr << "ERROR: Cannot find fission fragment yields file: " <<
    yieldsFilename << endl; exit(-1); }
  
  double norm=0.0;
  TKEmin=999; TKEmax=0; Amin=300; Amax=0;
  while (!yieldsFile.eof()) {
    yieldsFile >> A >> TKE >> yields;
    //    TKE = TKE + 1; // TODO: CHECK FOR ALL [CF-252SF only for now]
    YATKE[A][TKE] = yields;
    norm += yields;
    if (yields!=0.0) {
      if (TKE<TKEmin) TKEmin=TKE;
      if (TKE>TKEmax) TKEmax=TKE;
      if (A<Amin) Amin=A;
      if (A>Amax) Amax=A;
    }
  }
  
  yieldsFile.close();
  
  Amin = Amin-5;
  Amax = Amax+5;
  
  // symmetrize yields, if needed
  for (int i=Asym; i<Amax; i++) {
    for (int j=0; j<NUMTKE; j++) {
      YATKE[i][j] = YATKE[Ac-i][j];
    }
  }
  
  //-- renormalize yields Y[A][TKE] to 1.0
  for (int i=0; i<NUMA; i++) {
    for (int j=0; j<NUMTKE; j++) {
      YATKE[i][j] /= norm;
    }
  }
  
}

/*******************************************************************************
 * readYieldsAZTKE
 *------------------------------------------------------------------------------
 * Reads fission fragment yields in the form of [A, Z, TKE, Yield].
 ******************************************************************************/
void FissionFragments::readYieldsAZTKE (string yieldsFilename) {
  
  ifstream yieldsFile;
  int A, Z, TKE;
  double yields;
  
  yieldsFile.open(&yieldsFilename[0]);
  if (!yieldsFile) { cerr << "ERROR: Cannot find fission fragment yields [A,Z,TKE,Y] file: " <<
    yieldsFilename << endl; exit(-1); }
		
  double norm=0.0;
  TKEmin=999; TKEmax=0; Amin=300; Amax=0; Zmin=100; Zmax=0;
  while (!yieldsFile.eof()) {
    yieldsFile >> A >> Z >> TKE >> yields;
    YATKE[A][TKE] += yields;
    YZATKE[Z][A][TKE] += yields;
    YZA2[Z][A] += yields;
    norm += yields;
    if (yields!=0.0) {
      if (TKE<TKEmin) TKEmin=TKE;
      if (TKE>TKEmax) TKEmax=TKE;
      if (A<Amin) Amin=A;
      if (A>Amax) Amax=A;
      if (Z<Zmin) Zmin=Z;
      if (Z>Zmax) Zmax=Z;
    }
  }
  
  yieldsFile.close();
  
  // symmetrize yields, if needed
  for (int i=Asym; i<Amax; i++) {
    for (int j=0; j<NUMTKE; j++) {
      YATKE[i][j] = YATKE[Ac-i][j];
    }
  }
  
  //-- renormalize yields Y[A][TKE] to 1.0
  for (int i=0; i<NUMA; i++) {
    for (int j=0; j<NUMTKE; j++) {
      YATKE[i][j] /= norm;
    }
  }
  
  double Ymax;
  int jmax;
  // determine YZA from YZA2
  for (int i=Amin; i<=Amax; i++) {
    Ymax=-100;
    jmax=0;
    for (int j=1; j<NUMZ; j++) {
      if (YZA2[j][i]>Ymax) { Z0[i]=j; jmax=j; Ymax=YZA2[j][i]; }
    }
    Z0[i]=jmax;
  }
  
  for (int i=Amin; i<Amax; i++) {
    for (int j=0; j<=2*dZ; j++) {
      YZA[j][i] = YZA2[int(Z0[i])+j-dZ][i];
    }
  }
  
}


/*******************************************************************************
 * spinParameter
 *------------------------------------------------------------------------------
 * Returns the spin parameter B entering in the spin population formula
 * (2J+1)exp[-(J(J+1)/2B^2].
 * 'B' depends on A, Z, and excitation energy U.
 *
 * Moment of Inertia of rigid-body ellipsoid:
 *
 * I = 2/5 * M * R^2 * (1 + 0.31*beta + 0.44*beta^2 + ...)
 *
 * with M = A * m_nucleon
 *      R = r * A^(1/3)
 *
 ******************************************************************************/
double FissionFragments::spinParameter (double U, int A, int Z)
{
  
  if (beta2[Z][A]<-99.) { cerr << "Warning: spin parameter not defined!\n"; return (-1.); }
  
  // rigid-body moment of inertia
  double momInertia = 0.4*pow((double) A, 5./3.)*1.2*1.2*.5*(neutronMass+protonMass)*amuMeV;
  momInertia *= (1.0+0.31*beta2[Z][A]+0.44*pow(beta2[Z][A],2.));
  //	momInertia *= (1.0-0.8*exp(-0.693*U/5.0)); // energy-dependence
  
  // ad-hoc correction factor
  momInertia *= alphaI;
  
  int iU; // = findEnergyIndex(U, ldEnergyGrid);
  
  int ZAID = 1000*Z+A;
  int iZA=0; while (ZAID!=ZAIDlist[iZA]) { iZA++; }
  
  iU = findEnergyIndex(U,ldEnergyGrid2);
  
  double betaTemp = (temp[iZA][iU+1]-temp[iZA][iU])/(ldEnergyGrid2[iU+1]-ldEnergyGrid2[iU])*(U-ldEnergyGrid2[iU])+temp[iZA][iU];
  //  double betaTemp = temperatures[A][Z-int(Z0[A]+0.5)+dZ][iU] ;
  
  double Bspin = sqrt(momInertia*betaTemp/hbarc/hbarc);
  
  return (Bspin);
  
}

/*******************************************************************************
 * Finds the index corresponding to the value 'x0' in the array 'xarray'.
 ******************************************************************************/
template <int N> int FissionFragments::findEnergyIndex (double x0, double (&xarray)[N]) {
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
 *  Sample the probability of multi-chance fission and constructs the YA, YZ, YAZ
 ********************************************************************************/
void FissionFragments::constructYields (Yields *ffy){
  
  // reset first nucleus to original compound system
  incidentEnergy=incidentEnergy0 ;
  Ac=Ac0;
  Zc=Zc0;
  SnCompound=SnCompound0;
  
  double excitEnergy=ncl[0].max_energy;
  int nemit_now=-1;
  do{
    nemit_now=getPrefissionNeutrons (&excitEnergy);
  }while (nemit_now<0);
  
  if( nemit_now > 0 ){
    pfnInfo[index_save]=new double[nemit_now];
    for(int i=0;i<nemit_now;i++)
      pfnInfo[index_save][i]=energy_pfn[i];
  }
  
  Ac=Ac0-nemit_now;
  massExcessCompound = mass_excess(Zc, Ac);
  SnCompound = -mass_excess(Zc,Ac)+mass_excess(Zc,Ac-1)+mass_excess(0,1);
  incidentEnergy = excitEnergy-SnCompound;

	setSpinStrength(incidentEnergy);
	
  ffy->setFissioningSystem(ZAIDc-nemit_now, incidentEnergy);
  Amin=300; Amax=-1;
  maxYieldA=0.;
  for (int a=0; a<=NUMA; a++) {
    YA[a]=ffy->yieldA(a);
    if(YA[a]>0. && a<Amin)
      Amin=a;
    if(YA[a]>0. && a>Amax)
      Amax=a;
  }
  //renormalize the mass distribution
  double s=0.;
  for(int a=Amin; a<=Amax; a++){
    s+= YA[a];
  }
  for(int a=Amin; a<=Amax; a++){
    YA[a] /= s;
    if(YA[a]>maxYieldA)
      maxYieldA=YA[a];
  }
  ffy->setRescaleTKE(YA,Amin,Amax);
  computeWahlParametersForChargeDistribution();
  buildYZA();
  buildYZ ();
  
  for (int i=0; i<NUMA; i++) {
    maxYieldZ[i]=0.0;
    for (int j=0; j<NUMZ; j++) {
      if (YZA2[j][i]>maxYieldZ[i]) { maxYieldZ[i] = YZA2[j][i]; }
    }
  }
  
  return;
  
}


/*******************************************************************************
 *  Read multi-chance fission probability tables, calculated with CoH.
 *
 * Output:
 *
 * - barrier[]:   fission barrier heights for all fissioning nuclei
 * - emissprob[]: multi-chance fission probabilities
 *
 ********************************************************************************/
void FissionFragments::readMultichanceFissionData (void) {
  
#include <stdio.h>
#include <string.h>
  
  string f = DATADIR;
  f += "multichancefission.dat";
  ifstream fissdata (&f[0],ios::in);
  
  if (!fissdata) {
    cerr << "[readMultichanceFissionData] multi-chance fission probabilities data file could not be found" << endl;
    exit(-1);
  }
  
  int i;
  string line;
  int id, nn=0, nenergies=0;
  bool found = false;
  
  while (!found) {
    getline(fissdata,line);
    if (atoi(line.substr(0,5).c_str())==ZAIDc-1) {
      found = true;
      i = int(line.find(" "));
      id = atoi(line.substr(0,i).c_str());
      line = line.substr(i+1);
      i = int(line.find(" "));
      nenergies = atoi(line.substr(0,i).c_str());
      line = line.substr(i+1);
      i = int(line.find(" "));
      nn = atoi(line.c_str());
    }
    if (fissdata.eof()) {
      cerr << "[readMultichanceFissionData] could not find multi-chance fission data for this nucleus/energy" << endl;
      exit(-1);
    }
  }
  
  if (found) {
    
    barrier=new double[nn+1];
    for (i=0; i<=nn; i++) fissdata >> barrier[i];
    
    double Pf[nenergies][nn+1]; // fission probabilities
    double energies[nenergies];
    
    for (i=0; i<nenergies; i++) {
      fissdata >> energies[i];
      for (int j=0; j<=nn; j++) fissdata >> Pf[i][j];
    }
    
    emissprob = new double[nn+1];
    
    if (incidentEnergy<energies[0] || incidentEnergy>energies[nenergies-1]) {
      cerr << "[readMultichanceFissionData] could not find multi-chance fission data for this incident energy" << endl;
      exit(-1);
    }
    
    for (i=0; i<nenergies; i++) {
      if (incidentEnergy<=energies[i]) {
        double slope;
        for (int j=0; j<=nn; j++) {
          slope=(Pf[i][j]-Pf[i-1][j])/(energies[i]-energies[i-1]);
          if (slope!=0.0) {
            emissprob[j] = (incidentEnergy-energies[i-1])*slope+Pf[i-1][j];
          } else {
            emissprob[j] = Pf[i][j];
          }
        }
        break;
      }
    }
    
    // Find max. number of pre-fission neutrons
    int k=nn;
    for (int j=0; j<=nemit; j++) {
      if (emissprob[j]<=0.0) {
        k=j-1;
        break;
      }
    }
    nemit=k;
    
    // Find max. of emissprob (needed for later sampling)
    emissprob_max = 0.0;
    double s=0.0;
    for (int j=0; j<=nemit; j++) {
      s+=emissprob[j];
    }
    for (int j=0; j<=nemit; j++) {
      emissprob[j] /= s;
      if (emissprob[j]>emissprob_max) emissprob_max = emissprob[j];
    }
    
    fissdata.close();
    
    // pfn == "pre-fission neutrons"
    if (nemit>0){ // else?
      energy_pfn  = new double [nemit]; // energies of the pre-fission neutrons
      sp_pfn_mult = new int * [nemit];  // multiplicity-dependent pre-fission neutron spectra
      for (int i=0; i<nemit; ++i) {
        sp_pfn_mult[i] = new int [MAX_ENERGY_BIN];
      }
    }
    
  }
  
  return;
  
}



/*******************************************************************************
 * Returns the number of pre-fission neutrons. Also calculates their energies
 * (energy_pfn), spectrum (spn_pfn) and multiplicity-dependent spectrum
 * (spn_pfn_mult).
 *******************************************************************************/
int FissionFragments::getPrefissionNeutrons (double *excitEnergy) {
  
  if (nemit==0) return 0; // no pre-fission neutrons
  
  // Sample the multi-chance fission probabilities
  int n;
  do {
    n=(int) floor(genrand_real2()*(nemit+1));
  } while (emissprob[n]<emissprob_max*genrand_real1());
  if (n==0) return 0;
  
  // if the max. residual excitation energy is lower than the fission barrier
  // height in the final residual nucleus, fission cannot occur.
  double eps_max=ncl[n].excitation[0]-barrier[n];
  if (eps_max<0.) return -1;
  
  int kfmax=(int) ceil(eps_max/ENERGY_BIN);
  
  // sample the first energy from the spectrum sp1 (sum of evaporation and
  // pre-equilibrium), which was stored at the beginning of the calculation
  
  double max_prob=0.;
  
  for(int i=0;i<kfmax;i++){
    if(sp1[i]>max_prob)max_prob=sp1[i];
  }
  
  // Sample the first neutron-out spectrum 'sp1' (contains evaporation and
  // pre-equilibrium contributions)
  int kf;
  do{
    kf=(int) floor(genrand_real2()*kfmax);
  }while(sp1[kf]<max_prob*genrand_real3());
  
  energy_pfn[0]=energy_grid1[kf]; // continuum-to-continuum energy transfer
  
  for (int k=1; k<MAX_ENERGY_BIN; k++) {
    if (energy_pfn[0]<(k-0.5)*ENERGY_BIN) {
      sp_pfn[k-1]++;
      sp_pfn_mult[n-1][k-1]++;
      break;
    }
  }
  
  // Sample higher-order multi-chance fission probabilities (evaporation spectra
  // only)
  int kstart=kf;
  for (int i=2; i<=n; ++i) {
    
    evapInterface (pdt+neutron,tc,i-1,kstart,sp);
    eps_max=ncl[n].excitation[kstart]-barrier[n];
    kfmax=kstart+(int) floor(eps_max/ENERGY_BIN)+1;
    
    // << 1.0.6 >> Can this happen? ...
    if(eps_max<0.){
      cout << i << " nemit_now=" << n << " *eps_max=" << eps_max << " kstart=" << kstart << endl;
    }
    
    max_prob=0.;
    for(int k=kstart;k<kfmax;k++){
      if(sp[k]>max_prob)max_prob=sp[k];
    }
    
    do{
      kf=(int) (kstart+(kfmax-kstart)*genrand_real2());
    }while(sp[kf]<max_prob*genrand_real3());
    if(kf<0)cout << i << " nemit_now=" << n << " kstart=" << kstart << " kf=" << kf << endl;
    
    energy_pfn[i-1]= tc[kf].ecms;
    
    kstart=kf;
    for (int k=1; k<MAX_ENERGY_BIN; k++) {
      if (energy_pfn[i-1]<(k-0.5)*ENERGY_BIN) {
        sp_pfn[k-1]++;
        sp_pfn_mult[n-1][k-1]++;
        break;
      }
    }
  }
  *excitEnergy=ncl[n].excitation[kf];
  
  
  return n;
  
}

void FissionFragments::generateYieldsForSensitivityStudies (string YATKEfilename) { // for work with Ramona (LLNL) and Jorgen (LBNL)
  
  string basedir = "/Users/talou/calculus/FFdecay/98-Cf-252sf/data/GLS/YATKEsamples/";
  string fname = basedir + YATKEfilename;
  
  int i;
  int A[200];
  double TKE[200];
  
  // read mass values
  string AgridFile = basedir+"Agrid.dat";
  ifstream f;
  f.open(&AgridFile[0],ios::in);
  i=0;
  while (!f.eof()) {
    f >> A[i++];
  }
  f.close();
  
  int nA = i-1;
  
  // read TKE values
  string TKEgridFile = basedir+"TKEgrid.dat";
  f.open(&TKEgridFile[0],ios::in);
  i=0;
  while (!f.eof()) {
    f >> TKE[i++];
  }
  f.close();
  
  int nTKE = i-1;
  
  // read Y(A,TKE) values
  double **yatke;
  yatke = new double* [nA];
  for (int i=0; i<nA; i++) {
    yatke[i] = new double [nTKE];
  }
  
  ifstream dataFile;
  
  fname = basedir+ "test.dat";
  cout << fname << endl;
  
  string line;
  dataFile.open(&fname[0],ios::in);
  
  for (int i=0; i<nA; i++) {
    for (int j=0; j<nTKE; j++) {
      dataFile >> yatke[i][j];
    }
  }
  
  dataFile.close();
  
  // Build Y(Z|A)
  Amin = A[0];
  Amax = A[nA-1];
  
  ofstream file;
  file.open(&(basedir+"B")[0],ios::out);
  file << setprecision(3) << setiosflags(ios::scientific);
  int counter=0;
  for (int i=0; i<nA; i++) {
    file << "# [" << counter << "] " << i+72 << endl << endl;;
    for (int j=0; j<nTKE; j++) {
      file << j+140.0 << " " << yatke[i][j] << endl;
    }
    counter++;
    file << endl;
  }
  
  
  dZ=3;
  
  buildYZA();
  
  int ZAID, mass;
  double y1;
  double y2 [nTKE];
  double TKEmin[1000], TKEmax[1000]; // ZAID-dependence
  
  double eps=1.0e-7;
  bool found;
  
  ofstream out;
  
  
  out.open(&(basedir+"test.out")[0],ios::out);
  out << setprecision(3) << setiosflags(ios::scientific);
  
  double sum=0.0;
  counter=-1;
  for (int i=0; i<nA; i++) {
    mass = A[i];
    if (mass<126.0) continue;
    //    cout << i << " " << A[i] << " " << mass << " " << Z0[mass] << endl;
    for (int j=-dZ; j<=dZ; j++) {
      ZAID = 1000*int(Z0[mass]+j+0.5)+mass;
      counter++;
      y1 = YZA[j+dZ][A[i]];
      TKEmin[counter]=0.0;
      TKEmax[counter]=0.0;
      found = false;
      for (int k=0; k<nTKE; k++) {
        y2[k] = yatke[i][k]*y1;
        //				cout << i << " " << j << " " << k+140 << " " << y2[k] << " " << found << endl;
        if (found and y2[k]<eps) {
          TKEmax[counter]=TKE[k-1];
          break;
        }
        if (!found and y2[k]>=eps) {
          TKEmin[counter]=TKE[k];
          found=true;
        }
      }
      if (found and TKEmax[counter]==0.0) TKEmax[counter]=TKE[nTKE-1];
      
      if (TKEmin[counter]!=0.0) {
        for (int k=0; k<nTKE; k++) {
          if (y2[k]>=eps) {
            sum+=y2[k];
          }
        }
      }
    }
  }
  
  if (sum==0.0) { cerr << "sum=0!\n"; exit(-1); }
  sum=sum; // removed since dealing with heavy fragments only... 2.0; // to normalize yields to 2.0
  
  double s=0.0;
  counter=-1;
  for (int i=0; i<nA; i++) {
    mass = A[i];
    if (mass<126.0) continue;
    //    cout << i << " " << A[i] << " " << mass << " " << Z0[mass] << endl;
    for (int j=-dZ; j<=dZ; j++) {
      ZAID = 1000*int(Z0[mass]+j+0.5)+mass;
      y1 = YZA[j+dZ][A[i]];
      counter++;
      
      //			cout << counter << " " << ZAID << " " << TKEmin[counter] << " " << TKEmax[counter] << endl;
      
      if (TKEmin[counter]!=0.0) {
        out << ZAID << " " << int(TKEmin[counter]) << " " << int(TKEmax[counter]) << endl;
        for (int k=TKEmin[counter]-140; k<=TKEmax[counter]-140; k++) {
          y2[k] = yatke[i][k]*y1;
          s+=y2[k];
          out << y2[k] << " ";
        }
        out << endl;
      }
    }
  }
  
  cout << s/sum << endl;
  
  out.close();
  
  // Sep. 14, 2015 -------------------------------------------------------------
  // Check normalization
  // Solution: renormalized full yields to 2.0
  //
  // Issue with yields truncated because of (A,TKE)-grid, and not because of eps=1e-7.
  // Solution: change TKEgrid to [140:230] MeV
  //
  // Ordering is with A, not ZAID... Is this a problem for Jorgen?
  //
  // Could I limit the file to heavy fragments only? If symmetric, should save
  // half the space...
  // Solution: OK, done. Keep normalization at 1.0.
  
  exit(0);
  
  
}


// Ionel: new systematics to tune the alpha parameter below the 2nd-chance
// fission threshold, only for Pu-240 at this point.
void FissionFragments::setSpinStrength (double einc) {
	
	switch(ZAIDc){
  case 94240:
			
			alphaI=2.16-.633*exp(-.233*einc);
			//    cout << "alphaI="<< alphaI;
			
			break;
			
  default:
			break;
	}
	
}



