/*
 *  CGMF
 *  Version: 1.0.7.2
 *  Austin Carter
 *
 *  [ BrosaYields.cpp ]
 *
 */

// some libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <string.h>


#include "config.h"
//for the random numbers
#include "mt19937ar.h"
// has constants: pi, etc.
#include "physics.h"

// our namespace
using namespace std;

#include "BrosaYields.h"

// CONSTRUCTOR
BrosaYields::BrosaYields(string filein, int ZAID, int A_min){
 
  // initialize stuff
  Amin = A_min;
  Amax = NUMA - 1;
  TKEmin = 100;
  TKEmax = NUMTKE - 1;
  Zt = Zc;
  PE = 0.0;
  dZ = 3;
  PE = 0.0; 


  
	ZAIDc = ZAID;
	Zc=int(ZAIDc/1000.0);
  	Ac=ZAIDc-1000*Zc;
  	Zt = Zc;
  	PE = 0.0;
	
	// populate A and TKE vectors
	for (int i = 0; i < NUMA; i++){
		A[i] = i;
	}

	for (int i = 0; i < NUMTKE; i++){
		TKE[i] = i;
	}

	// put parameter data into vectors from file
	readfile(filein);
  //cout << "readfile good" << endl;

	// build the YA array
	buildYA();
  //cout << "buildYA good" << endl;
	// build the Y(TKE|A)_m array
	buildYTKEA();
  //cout << "buildYTKEA good" << endl;

	// build the Y(Z|A) array
	computeWahlParameters();
  //cout << "computeWahlParameters good" << endl;
	buildYZA();
  //cout << "buildYZA good" << endl;




}

// Destructor
BrosaYields::~BrosaYields () {}

void BrosaYields::readfile(string inputFilename){

  	string f = DATADIR;
	f += inputFilename;
	ifstream brosadata (&f[0],ios::in);

	if (!brosadata) {
    	cerr << "[BrosaYields::readfile(inputFilename)] Brosa Yields parameters data file could not be found" << f <<  endl;
    	exit(-1);
  	}

  	string line;
  	string strtmp;
  	double dubtmp;

  	while(getline(brosadata,line)){
  		stringstream ss(line);
  		//found the w
  		if(int(line.find("w")) == 0){
  			ss >> strtmp;
  			for (int i = 0; i < NUMMODE; i++){ss >> dubtmp; w.push_back(dubtmp);}
  		}
  		if(int(line.find("dmin")) == 0){
  			ss >> strtmp;
  			for (int i = 0; i < NUMMODE; i++){ss >> dubtmp; dmin.push_back(dubtmp);}
  		}
  		if(int(line.find("dmax")) == 0){
  			ss >> strtmp;
  			for (int i = 0; i < NUMMODE; i++){ss >> dubtmp; dmax.push_back(dubtmp);}
  		}
  		if(int(line.find("ddec")) == 0){
  			ss >> strtmp;
  			for (int i = 0; i < NUMMODE; i++){ss >> dubtmp; ddec.push_back(dubtmp);}
  		}
  		if(int(line.find("Abar")) == 0){
  			ss >> strtmp;
  			for (int i = 0; i < NUMMODE; i++){ss >> dubtmp; Abar.push_back(dubtmp);}
  		}
  		if(int(line.find("sigA")) == 0){
  			ss >> strtmp;
  			for (int i = 0; i < NUMMODE; i++){ss >> dubtmp; sigA.push_back(dubtmp);}
  		}
  	}

}

void BrosaYields::buildYA(){
  std::fill_n (YA, NUMA, 0.0);
	for(int i = Amin; i < NUMA; i++){
		YA[i] = 0.0;
		for (int j = 0; j < 3; j++){
			YA[i] += 1/sqrt(8* pi * sigA[j]*sigA[j])*
			(exp(-(A[i] - Abar[j])*(A[i] - Abar[j])/2./sigA[j]/sigA[j]) +
				exp(-(A[i] - Ac +Abar[j]) * (A[i] - Ac +Abar[j]) / 2. / sigA[j] / sigA[j]));
		}
	}
}

void BrosaYields::buildYTKEA(){
	// okay, the way this large array is structured is as follows:
	// 3 modes that act as pages; the first index element
	// then for the 2d array (contained within each mode):
	//     TKE ---------->
	// A
	// |
	// |
	// |
	// ^
	// where the distributions are normalized to 1.0 along the row which are sampled in sampleTKEA

	double T;
	double TKEAtotal;

	for (int i = 0; i < NUMMODE; i++){
		for (int j = 0; j < NUMA; j++){
      // 0 the array
      for (int k = 0; k < NUMTKE; k++){
        YTKEA[i][j][k]= 0.0;
      }
			TKEAtotal = 0.;
			for (int k = TKEmin; k < NUMTKE; k++){
				T = ((double)Zc/(double)Ac)*((double)Zc/(double)Ac) * ((double)Ac - (double)A[j])*(double)A[j] * e2() / (double)TKE[k] - dmin[i];
				if(T > 0.){
					YTKEA[i][j][k] = (200./(double)TKE[k])*(200./(double)TKE[k])* exp(2. * (dmax[i] - dmin[i])/ ddec[i] - 
						T/ddec[i] - (dmax[i] - dmin[i])*(dmax[i] - dmin[i]) /T /ddec[i])
						* (double)YA[j];
				} else {
					YTKEA[i][j][k] = 0.;
				}
				TKEAtotal += YTKEA[i][j][k];

			}

			// avoid the nans
			if (TKEAtotal == 0.0){
				TKEAtotal = 0.01;
			}
			// normalize to 1.0 for each row
			for (int k = TKEmin; k < NUMTKE; k++){
				YTKEA[i][j][k]= YTKEA[i][j][k] / TKEAtotal;
			}
		}
	}

}

int BrosaYields::samplemode(){
	// get a random number
	double r1 = genrand_real3();

	// sample probability mass dist for mode
	double wtot = 0.;
	int i_mode;
	for(int i = 0; i < NUMMODE; i++){
		wtot += w[i];
		if (r1 <= wtot){
			i_mode = i;
			break;
		}
	}
	return i_mode;

}


int BrosaYields::sampleA(int i_mode){

	// sample the gaussian to find A for the HEAVY Fragment
	double r2 = genrand_real3();
	double r3 = genrand_real3();

	double A_s = Abar[i_mode] + sigA[i_mode] * sqrt(-2. * log(r2)) * cos(2. * pi * r3);


	// ALWAYS SAMPLES THE HEAVY FRAGMENT

	// return the heavy fragment as rounded integer
	return int(A_s+0.5);
}

int BrosaYields::sampleTKEA(int A, int i_mode){

	// we use YTKEA as a probability mass distribution
	// find the index of A
	int i_A = (int)(A);

	double r4 = genrand_real3();

	double KEtot = 0.0;

	for (int i_TKE = 0; i_TKE < NUMTKE; i_TKE++){
		KEtot += YTKEA[i_mode][i_A][i_TKE];
		if (r4 <= KEtot){
      //cout << "sampled TKE|A" << endl;
			return (int)(i_TKE);
		}
	}

	// should not get here
	cerr << "ERROR: Could not sample KE distribution in BrosaModes class." << endl;
    exit(-1);
}

int BrosaYields::sampleZA(int A){
	int i_A = int(A);
	double r5 = genrand_real3();

	// normalize the column vector in YZA (because we sample Z based on A)
	double YZA_tot = 0.0;
	for (int i = 0; i < NUMZ; i++){
		YZA_tot += YZA2[i][i_A];
	}
	//create a new vector to be sampled from
	double YZA_temp[NUMZ];
	for (int i = 0; i<NUMZ; i++){
		YZA_temp[i] = YZA2[i][i_A] / YZA_tot;
	}

	// sample the probability mass distribution
	double Ztot = 0.0;
	for (int i_Z = 0; i_Z < NUMZ; i_Z++){
		Ztot += YZA_temp[i_Z];
		if(r5 <= Ztot){
			return i_Z;
		}
	}
	// should never get here
	cerr << "ERROR: Could not sample Z distribution in BrosaModes class." << endl;
    exit(-1);
}

void BrosaYields::buildYZA(){
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

void BrosaYields::computeWahlParameters(){
	double sigz140t, delz140t, Fz140t, Fn140t, sigzSLt, delzSLt, FzSLt, SL50t,
  sigz50t, delzmaxt, sigzSLWt, delzSLWt, FzSLWt, FnSLWt;
  
  double delz[NUMA], sigz[NUMA], Fz[NUMA], Fn[NUMA];
  
  double sigz140[5] = {0.566, 0.0, 0.0064, 0.0109, 0.0};
  double delz140[5] = {-0.487, 0.0, 0.0180, 0.0, -0.00203};
  //  double Fz140[5]   = {1.207, 0.0, -0.0420, 0.0, 0.0022};
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

