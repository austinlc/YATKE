/*
 *  CGMF
 *  Version: 1.0.6
 *
 *  [ FissionFragments.h ]
 * 
 *  Generates primary fission fragment yields to be sampled in CGMF. It either
 *  reads data from input file, or produce them from systematics.
 *
 */

#ifndef __BROSAYIELDS_H__
#define __BROSAYIELDS_H__
#endif

// screw arrays, i'm using std::vectors
#include <vector>
#include <string.h>

using namespace std;

class BrosaYields
{
	/* PRIVATE */
private:
	// The data that is important for BrosaYields
	// the parameters that are necessary for the yield constructions
	// these will be "push_back(ed)" in the constructor
	

	/* PUBLIC */
public:




	vector<double> w;
	vector<double> Abar;
	vector<double> sigA;
	vector<double> dmin;
	vector<double> dmax;
	vector<double> ddec;

	int ZAIDc, Ac, Zc; // Original (before any pre-fission neutron emission) compound fissioning nucleus

	

	

	//   7/26/17 values - Cf 252


	static const int NUMMODE = 3;
	static const int NUMA= 300;
	static const  int NUMTKE = 300;

  	static const  int    NUMdZ  =   21; // [-dZ:+dZ] if dZ=10 for charge distribution around most probable Zp[A] 
  	static const  int    NUMZ=  100; // number of charges Z


	// vectors of A and TKE
	// initialized in constructor
	int A[NUMA];
	int TKE[NUMTKE];

	int Amin;
	int Amax;
	int TKEmin;
	int TKEmax;
	int Zmin, Zmax;
	int Zt;
	double PE;

	// Wahl's parameters

	int dZ;
  	double Z0[NUMA];
  	double sZ0[NUMA];
  	double FZZ0[NUMA];
  	double FNZ0[NUMA];

  	double sigmaZ;

	double YA[NUMA];
  	double YZA[NUMdZ][NUMA];
  	double YZA2[NUMZ][NUMA];


	// array of TKE|A for each mode 
	// first index is mode number
	// second index is A (row)
	// third index is probability of TKE|A
	double YTKEA[NUMMODE][NUMA][NUMTKE];

	static double e2(){
		return 1.4399643929;// MeV fm
	}

// Constructor: when implemented, a file is read with the necessary parameters for the Brosa modes
	// necessary data in input file: h (equivalently the weight w), Abar, sigA, dmax, dmin, ddec
	// for each mode
	BrosaYields (string inputFilename, int ZAID, int A_min);

	~BrosaYields (void);

	// function to read file for input parameters
	void readfile(string inputFilename);

	// function to build YTKEA
	void buildYTKEA();



	void buildYA();

	// A method that would be necessary would be finding out what mode we are dealing with
	int samplemode(); // returns a mode number according to the order of that in the input file
	// I intend for the order to be S2, S1, SL

	// let's have a public method that can sample an A
	// ALWAYS SAMPLES THE HEAVY FRAGMENT
	int sampleA(int i_mode);

	// And another method that samples Y(TKE|A)
	int sampleTKEA(int A, int mode);

	void computeWahlParameters();

	void buildYZA();

	int sampleZA(int A);

};

