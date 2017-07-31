
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>

using namespace std;

#include "config-ff.h"
#include "ffd.h"

// --------------------------------------------------------------------------------

const int MAX_NUMBER_EVENTS = 10000;
const int MAX_MULT = 50;

// --------------------------------------------------------------------------------

int countsA [NUMA];

double spectrumEnergyGrid[NUMBER_SPECTRUM_ENERGY_GRID];
double spectrumEnergySteps[NUMBER_SPECTRUM_ENERGY_GRID];

struct emittedParticleType
{
  
  char particleType; // 'n' (neutron) or 'g' (gamma)
  
  // to keep track of total number of emitted particles, also for light and heavy fragments
  int counts, lfCounts, hfCounts;
  
  // multiplicity and multiplicity distributions
  double nubar, lfNubar, hfNubar;
  double Pnu [MAX_MULT], lfPnu [MAX_MULT], hfPnu [MAX_MULT];   // total, light fragment and heavy fragment P(nu)
  
  // average energies in center-of-mass and lab frames
  double Ecm, lfEcm, hfEcm;
  double Elab, lfElab, hfElab;
  
	double Thetacm, Thetalab;
	
  // Angle between neutron and light fragment
  double neutronLightFragmentAngularDistribution [NUMANGLES];
  double neutronNeutronAngularDistribution [NUMANGLES];
  
  // results as a function of fragment mass
  double EcmA [NUMA], ElabA [NUMA], nuA [NUMA];
  double EcmAcut [NUMA], nuAcut [NUMA];
  
  // results as a function of TKE
  double EcmTKE [NUMTKE], ElabTKE [NUMTKE], nuTKE [NUMTKE];
  
  // spectra in the center-of-mass and lab frames
  double cmSpectrum [NUMBER_SPECTRUM_ENERGY_GRID];
  double labSpectrum [NUMBER_SPECTRUM_ENERGY_GRID];
	
  // average spectra for light and heavy fragments
  double *cmLightFragmentSpectrum, *cmHeavyFragmentSpectrum;
  double *labLightFragmentSpectrum, *labHeavyFragmentSpectrum;
  
  // c.m. spectrum as a function of mass of primary fragment
  double **cmSpectrumA;
	
  // average c.m. energy as a function of mass and TKE
  double **EcmATKE;
  int **countsATKE;
  double **nubarATKE;
  
  // exclusive spectra for given multiplicities
  double **cmExclusiveSpectra; // function of multiplicity and Eout
	
  // two-particle energy correlations
  double **energyCorrelations;
  
  // two-particle angular correlations
  double **angularCorrelations;
  
  // Y(E,theta) in the lab frame
  double **labEnergyAngleDistribution;
  
};


struct emissionType {
	int multiplicity;
	double cmEnergies   [MAX_NUMBER_PARTICLES];
	double labEnergies  [MAX_NUMBER_PARTICLES];
	double cmAngles     [MAX_NUMBER_PARTICLES];
	double labAngles    [MAX_NUMBER_PARTICLES];
	int transitionTypes [MAX_NUMBER_PARTICLES];
};

struct fragmentEventType {
	int A;
	int Z;
	double KE;
	double Ui;
	float Ji;
	int Pi;
	emissionType emissions[3]; // neutrons [0]. gammas [1], and internal conversion [2]
};

double YZ [NUMZ];
double YA [NUMA];

double YUl   [NUMBER_EXCITATION_ENERGY_GRID];
double YUh   [NUMBER_EXCITATION_ENERGY_GRID];
double YUtot [NUMBER_EXCITATION_ENERGY_GRID];

double YKEl [NUMBER_EXCITATION_ENERGY_GRID];
double YKEh [NUMBER_EXCITATION_ENERGY_GRID];
double YTKE [NUMBER_EXCITATION_ENERGY_GRID];

double incidentNeutronEnergy; // [MeV]

// --------------------------------------------------------------------------------

static void		readHistoryFile		(string);
static void		initArrays			(void);
//static void     computeResults      (void);
static void		computeResults		(emittedParticleType *);
//static void		updateHistograms	(emittedParticleType *);
static void		displaySummary		(string);
static void     saveResultsInFile   (emittedParticleType *, string);
static void     initEmittedParticles (emittedParticleType *);

// --------------------------------------------------------------------------------

fragmentEventType lf[MAX_NUMBER_EVENTS], hf[MAX_NUMBER_EVENTS];

emittedParticleType neutrons, gammas;

int numberEvents; // Total number of fission events read from file


/**********************************************************/
/*                         MAIN                           */
/**********************************************************/
int main(int argc, char *argv[])
{
	
	
  initArrays ();
  initEmittedParticles(&neutrons);
  initEmittedParticles(&gammas);
 
  readHistoryFile ("histories.CGMF");

//  computeResults ();
  computeResults (&neutrons);
  computeResults (&gammas);
  
  saveResultsInFile (&neutrons, "neutrons.cgmf.test");
  saveResultsInFile(&gammas, "gammas.cgmf.test");
  
  displaySummary ("");
	
/*	for (int i=0; i<NUMA; i++) {
		cout << i << " " << neutrons.nuA[i] << " " << gammas.nuA[i] << " " << countsA[i] << endl;
	}*/
	
	
  cout << "\n*** WELL DONE! ***\n\n";
	
  return 0;
	
}

/*==============================================================================================
 Save results for a type of particles (n or g) in file (user input).
 ===============================================================================================*/
void saveResultsInFile (emittedParticleType *particles, string outputFilename)
{
  
  ofstream fout;
  fout.open (&outputFilename[0]);
  
  fout << setprecision(4) << setiosflags(ios::scientific);
  fout << "#\n# CGMF Results\n#\n";
  
  int gnuplotIndex=-1;
  
  //-- Laboratory PFNS
  fout << "# [gnuplot #" << ++gnuplotIndex << "] lab spec, cm spec, LF lab spec, LF cm spec, HF lab spec, HF cm spec\n\n";
  for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) {
    cout << i << " " << spectrumEnergyGrid[i] << "\n";
    fout << spectrumEnergyGrid[i] << " " << particles->labSpectrum[i] << " " << particles->cmSpectrum[i] << " " <<
    particles->labLightFragmentSpectrum[i] << " " << particles->cmLightFragmentSpectrum[i] << " " <<
    particles->labHeavyFragmentSpectrum[i] << " " << particles->cmHeavyFragmentSpectrum[i] << "\n";
  }
  
  //-- cm PFNS and multiplicity-dependent exclusive spectra
  fout << "\n# [gnuplot #" << ++gnuplotIndex << "] c.m. spectrum\n";
  fout << "# Energy    Total       nu=1      nu=2      nu=3       nu=4      nu=5\n";
  fout << "# (MeV)     (1/MeV)\n\n";
  for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++)
    fout << spectrumEnergyGrid[i] << " " << particles->cmSpectrum[i] << " " <<
    particles->cmExclusiveSpectra[1][i] << " " << particles->cmExclusiveSpectra[2][i] << " " <<
    particles->cmExclusiveSpectra[3][i] << " " << particles->cmExclusiveSpectra[4][i] << " " <<
    particles->cmExclusiveSpectra[5][i] << "\n";
  
  //-- Multiplicity distributions P(nu), and same for light and heavy fragments resp.
  fout << "\n# [gnuplot #" << ++gnuplotIndex << "] P(nu) P_LF(nu) P_HF(nu)\n\n";
  for (int i=0; i<20; i++)
    fout << i << " " << particles->Pnu[i] << " " << particles->lfPnu[i] << " " << particles->hfPnu[i] << "\n";
  
  //-- Fragment mass-dependent results: mass Yields A Y(A), <nu>(A), <Ecm>(A), <Ecm>(A) w/ cut, etc.
  fout << "\n# [gnuplot #" << ++gnuplotIndex << "] A  Y(A)   <nu>(A)   <Ecm>(A)  <Ecm>(A)|cut\n\n";
  for (int i=0; i<NUMA; i++)
    fout << i << " " << YA[i] << " " << particles->nuA[i] << " "
    << particles->EcmA[i] << " " << particles->EcmAcut[i] << "\n";
  
  //-- Z-dependent results: Y(Z), ...
  /*out << "\n#  [gnuplot #" << ++gnuplotIndex << "] Z   Y(Z) \n\n";
   for (int i=0; i<NUMZ; i++) {
   out << i << " " << YZ[i] << "\n";
   }*/
  
  
  //-- Angular distributions
  fout << "\n# [gnuplot #" << ++gnuplotIndex << "] Theta   n-LF    n-n\n";
  fout << "#       (deg)   (arb. u.)    (arb. u.)\n\n";
  if (particles->particleType=='n') {
    for (int i=0; i<NUMANGLES; i++) {
      fout  << i*dTheta << " " << particles->neutronLightFragmentAngularDistribution[i] << " " << particles->neutronNeutronAngularDistribution[i] << " " << sin(i*dTheta) << "\n";
    }
  } else {
    fout << "0 0 0\n";
  }
  
  //-- Y(TKE), <nu>(TKE), <Ecm>(TKE), <Elab>(TKE)
  fout << "\n# [gnuplot #" << ++gnuplotIndex << "] TKE   Y(TKE)   <nu>  <Ecm>  <Elab>\n";
  fout << "#             (MeV)   (n/f)  (MeV)   (MeV)   (MeV) \n\n";
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    fout << excitationEnergyGrid[i] << " " << YTKE[i] << " " << particles->nuTKE[i] <<
    " " << particles->EcmTKE[i] << " " << particles->ElabTKE[i] << "\n";
  }
  
  //-- Einc, <nu>, <nu_l>, <nu_h>, <Ecm>, <Ecm_l>, <Ecm_h>, <Elab>, <Elab_l>, <Elab_h>
  fout << "\n# [gnuplot #" << ++gnuplotIndex << "] Einc, <nu>, <nu_l>, <nu_h>, <Ecm>, <Ecm_l>, <Ecm_h>, <Elab>, <Elab_l>, <Elab_h>\n\n";
  fout << incidentNeutronEnergy << " " << particles->nubar << " " << particles->lfNubar
  << " " << particles->hfNubar << " " << particles->Ecm << " " << particles->lfEcm
  << " " << particles->hfEcm << " " << particles->Elab << " " << particles->lfElab
  << " " << particles->hfElab << "\n";
  
  fout.close();
  
}

/*==============================================================================================
 Display a summary of the main results.
 ===============================================================================================*/
void displaySummary (string outputFilename)
{
	
	std::ostream* fp = &cout;
	std::ofstream fout;
	
	if (outputFilename!="") {
		fout.open(&outputFilename[0]);
		fp = &fout;
	}
	
	*fp << setprecision(4) << setiosflags(ios::scientific);
	
	*fp << "*******************************************************" << endl;
	*fp << "             C G M F  Summary of Results               " << endl;
	*fp << "*******************************************************" << endl;
/*	if (incidentNeutronEnergy!=0.0) {
		*fp << "                n ( " << incidentNeutronEnergy << " MeV ) + " << Zc*1000+Ac-1 << endl;
	} else {
		*fp << "      Spontaneous fission of " << Zc*1000+Ac << endl;
	}*/
	
	*fp << "*******************************************************" << endl << endl;
	
	*fp << "-------------------------------------------------------" << endl;
	*fp << " Neutrons \n";
	*fp << "-------------------------------------------------------" << endl;
	
	*fp << setiosflags(ios::fixed) << setprecision(4);
	
	*fp << "<nu>      = " << neutrons.nubar   << endl;
	*fp << "<nu> (lf) = " << neutrons.lfNubar << endl;
	*fp << "<nu> (hf) = " << neutrons.hfNubar << endl << endl;
	
	*fp << "<Ecm>      = " << neutrons.Ecm   << " MeV " << endl;
	*fp << "<Ecm> (lf) = " << neutrons.lfEcm << " MeV " << endl;
	*fp << "<Ecm> (hf) = " << neutrons.hfEcm << " MeV " << endl << endl;
	
	*fp << "<Elab>      = " << neutrons.Elab   << " MeV " << endl;
	*fp << "<Elab> (lf) = " << neutrons.lfElab << " MeV " << endl;
	*fp << "<Elab> (hf) = " << neutrons.hfElab << " MeV " << endl << endl;
	
	*fp << endl;
	*fp << "-------------------------------------------------------" << endl;
	*fp << " Gamma Rays \n";
	*fp << "-------------------------------------------------------" << endl;
	
	*fp << "<Ng>      = " << gammas.nubar   << endl;
	*fp << "<Ng> (lf) = " << gammas.lfNubar << endl;
	*fp << "<Ng> (hf) = " << gammas.hfNubar << endl << endl;
	
	*fp << "<Ecm>      = " << gammas.Ecm   << " MeV " << endl;
	*fp << "<Ecm> (lf) = " << gammas.lfEcm << " MeV " << endl;
	*fp << "<Ecm> (hf) = " << gammas.hfEcm << " MeV " << endl << endl;
	
	if (fout.is_open()) { fout.close(); }
	
}


/*==============================================================================================
 Compute final results for neutrons and gammas.
 ===============================================================================================*/
void computeResults (emittedParticleType *particles)
{
	int i;
	double sum;
	int counts, lfCounts, hfCounts;

	
	counts   = particles->nubar;
	lfCounts = particles->lfNubar;
	hfCounts = particles->hfNubar;
	
	// Average Quantities ------------------------------------------------------------------------

	particles->nubar   /= numberEvents;
	particles->lfNubar /= numberEvents;
	particles->hfNubar /= numberEvents;
	
	particles->Ecm     /= counts;
	particles->lfEcm   /= lfCounts;
	particles->hfEcm   /= hfCounts;
		
	// Multiplicity Distributions ----------------------------------------------------------------

	sum = 0.0;
	for (i=0; i<MAX_MULT; i++) { sum += particles->Pnu[i]; }
	for (i=0; i<MAX_MULT; i++) { particles->Pnu[i] /= sum; }

	sum = 0.0;
	for (i=0; i<MAX_MULT; i++) { sum += particles->lfPnu[i]; }
	for (i=0; i<MAX_MULT; i++) { particles->lfPnu[i] /= sum; }

	sum = 0.0;
	for (i=0; i<MAX_MULT; i++) { sum += particles->hfPnu[i]; }
	for (i=0; i<MAX_MULT; i++) { particles->hfPnu[i] /= sum; }

	// As a function of fission fragment mass -----------------------------------------------------

	for (i=0; i<NUMA; i++) { if (countsA[i]!=0.0) { particles->nuA[i] /= countsA[i]; } else { particles->nuA[i]=0.0; } }
	
	
}



/*==============================================================================================
 Initialize most arrays to zero.
 ================================================================================================*/
void initArrays ()
{
  
  for (int i=0; i<NUMA; i++) { countsA[i]=0; }

	for (int j=0; j<3; j++) {
		for (int i=0; i<MAX_NUMBER_EVENTS; i++) {
			lf[i].emissions[j].multiplicity=0;
			for (int k=0; k<MAX_MULT; k++) {
				lf[i].emissions[j].cmEnergies[k]=0;
				lf[i].emissions[j].labEnergies[k]=0;
				lf[i].emissions[j].cmAngles[k]=0;
				lf[i].emissions[j].labAngles[k]=0;
				hf[i].emissions[j].cmEnergies[k]=0;
				hf[i].emissions[j].labEnergies[k]=0;
				hf[i].emissions[j].cmAngles[k]=0;
				hf[i].emissions[j].labAngles[k]=0;
      }
		}
	}
	
  // define spectrum energy grid -----------------------------------------------
  spectrumEnergyGrid[0]=1.0e-5;
  for (int i=1; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) {
    spectrumEnergyGrid[i] = spectrumEnergyGrid[i-1]*1.025;
    spectrumEnergySteps[i] = spectrumEnergyGrid[i]-spectrumEnergyGrid[i-1];
    //  cout << i << " " << spectrumEnergyGrid[i] << " " << spectrumEnergySteps[i] << "\n";
  }
  spectrumEnergySteps[0] = spectrumEnergySteps[1];
  
	return;
}



/*==============================================================================================
 Read all Monte Carlo histories from CGMF output file 
 Lines go by pairs: one for the light fragment, and the following line for the heavy fragment
================================================================================================*/

void readHistoryFile (string filename)
{
	
	fragmentEventType *event;
	int A,Z;
	int Ji, Pi;
	double KE,Ui;
	int nMult, gMult, icMult;
	double Ecm, Elab;
	double Thetacm, Thetalab;

	
	ifstream historyFile;
	historyFile.open (&filename[0]);
	
	int c=-1;
	while (!historyFile.eof()) {

		c++; // increment event counter
		
		for (int f=0; f<=1; f++) { // loop over light and heavy fragments -------------------
			
			if (f==0) { event = &lf[c]; } else { event = &hf[c]; }
			
//			historyFile >> event->A >> event->Z >> event->Ui >> event->Ji >> event->Pi >> event->KE >>
//				event->emissions[0].multiplicity >> event->emissions[1].multiplicity >> event->emissions[2].multiplicity;

			historyFile >> A >> Z >> Ui >> Ji >> Pi >> KE >> nMult >> gMult >> icMult;
			
			neutrons.nubar += nMult;
			if (f==0) { neutrons.lfNubar += nMult; }
			
			
			for (int i=0; i<nMult; i++) { historyFile >> event->emissions[0].transitionTypes[i]; }
			for (int i=0; i<gMult; i++) { historyFile >> event->emissions[1].transitionTypes[i]; }
			
			for (int i=0; i<nMult; i++) {
				historyFile >> Ecm >> Thetacm >> Elab >> Thetalab;
				neutrons.Ecm += Ecm;
				neutrons.Thetacm += Thetacm;
				neutrons.Elab += Elab;
				neutrons.Thetalab += Thetalab;
			}

			for (int i=0; i<gMult; i++) {
				historyFile >> Ecm >> Thetacm >> Elab >> Thetalab;
				gammas.Ecm += Ecm;
				gammas.Thetacm += Thetacm;
				gammas.Elab += Elab;
				gammas.Thetalab += Thetalab;
			}

			cout << nMult << " " << gMult << "\n";
			event = NULL;
			
		} // end loop over light and heavy fragments ------------------------------------------

		//if (c>0) { exit(0); }
		
		if (c>MAX_NUMBER_EVENTS) { cerr << "ERROR: number of fission events greater than MAX_NUMBER_EVENTS!\n"; cout << c << endl; exit(-1); }
		
	} // end reading file
	historyFile.close();
	numberEvents = c;
	
	cout << "Number fission events read in: " << numberEvents << endl;

/*
	fragmentEventType *lightFragment, *heavyFragment;
	
	for (int c=0; c<numberEvents; c++) {
		
		lightFragment = &lf[c];
		neutrons.lfNubar += lightFragment->emissions[0].multiplicity;

		
		
		
		gammas.lfNubar += lightFragment->emissions[1].multiplicity;
		
		
	}
	
	*/
	
}

/*******************************************************************************
 Initialize the neutron and gamma-ray structures (emittedParticleType).
 ******************************************************************************/

void initEmittedParticles (emittedParticleType *particles) {
  
  particles->counts   = 0;
  particles->lfCounts = 0;
  particles->hfCounts = 0;
  
  particles->Ecm   = 0.0;
  particles->lfEcm = 0.0;
  particles->hfEcm = 0.0;
  
  particles->Elab   = 0.0;
  particles->lfElab = 0.0;
  particles->hfElab = 0.0;
    
  for (int i=0; i<NUMA; i++) {
    particles->EcmA[i] = 0.0;
    particles->EcmAcut[i] = 0.0;
    particles->ElabA[i] = 0.0;
    particles->nuA[i] = 0.0;
    particles->nuAcut[i] = 0.0;
  }
    
  for (int i=0; i<NUMBER_EXCITATION_ENERGY_GRID; i++) {
    particles->EcmTKE[i] = 0.0;
    particles->ElabTKE[i] = 0.0;
    particles->nuTKE[i] = 0.0;
  }
  
  for (int i=0; i<NUMMULT; i++) {
    particles->Pnu[i] = 0.0;
    particles->lfPnu[i] = 0.0;
    particles->hfPnu[i] = 0.0;
  }
  
  for (int i=0; i<NUMBER_SPECTRUM_ENERGY_GRID; i++) {
    particles->cmSpectrum[i] = 0.0;
    particles->labSpectrum[i] = 0.0;
  }
  
  particles->cmLightFragmentSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  particles->cmHeavyFragmentSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  
  particles->labLightFragmentSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  particles->labHeavyFragmentSpectrum = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  
  particles->cmSpectrumA = new double * [NUMA];
  for (int i=0; i<NUMA; i++) {
    particles->cmSpectrumA[i] = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  }
  
  particles->cmExclusiveSpectra = new double * [10];
  for (int i=0; i<10; i++) {
    particles->cmExclusiveSpectra[i] = new double [NUMBER_SPECTRUM_ENERGY_GRID];
  }
  
/*  particles->energyCorrelations = new double * [NUMBER_OF_PRIVATE_GRID];
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) {
    particles->energyCorrelations[i] = new double [NUMBER_OF_PRIVATE_GRID];
  }
  */
  
  particles->angularCorrelations = new double * [NUMANGLES];
  for (int i=0; i<NUMANGLES; i++) {
    particles->angularCorrelations[i] = new double [NUMANGLES];
  }
  
  /*
  particles->neutronLightFragmentAngularDistribution = new double [NUMANGLES];
  std::fill_n (particles->neutronLightFragmentAngularDistribution, NUMANGLES, 0.0);
  
  particles->neutronNeutronAngularDistribution = new double [NUMANGLES];
  std::fill_n (particles->neutronNeutronAngularDistribution, NUMANGLES, 0.0);
  
  particles->labEnergyAngleDistribution = new double * [NUMBER_OF_PRIVATE_GRID];
  for (int i=0; i<NUMBER_OF_PRIVATE_GRID; i++) {
    particles->labEnergyAngleDistribution[i] = new double [NUMANGLES];
  }
   */
  
  
  // To compare with Vorobyev's data <Ecm>(A,TKE)
  particles->EcmATKE    = new double * [NUMA];
  particles->nubarATKE  = new double * [NUMA];
  particles->countsATKE = new int    * [NUMA];
  for (int i=0; i<NUMA; i++) {
    particles->EcmATKE[i]    = new double [NUMTKE];
    particles->nubarATKE[i]  = new double [NUMTKE];
    particles->countsATKE[i] = new int    [NUMTKE];
  }
  
}


