//
//  terminate.cpp
//  CGMF
//
//  Created by Patrick Talou on 9/2/15.
//
//

#include <stdio.h>
#include <string>
#include <iostream>
#include <stdlib.h>

using namespace std;

#include "terminate.h"
#include "cgm.h"


/**********************************************************/
/*     Emergency Stop                                     */
/**********************************************************/
int cgmTerminateCode(std::string msg)
{
	/*** Release global storage */
	cgmDeleteAllocated();
	
	/*** Exit code */
	cerr << "ERROR     :" << msg << endl;
	exit(-1);
}

int cgmTerminateCode(std::string msg, int n)
{
	/*** Release global storage */
	cgmDeleteAllocated();
	
	/*** Exit code */
	cerr << "ERROR     :" << msg << " : " << n << endl;
	exit(-1);
}

int cgmTerminateCode(std::string msg, double x)
{
	/*** Release global storage */
	cgmDeleteAllocated();
	
	/*** Exit code */
	cerr << "ERROR     :" << msg << " : " << x << endl;
	exit(-1);
}

/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void cgmDeleteAllocated()
{
	delete [] ncl;
	for(int i=0 ; i<SPECTRA_OUTPUT ; i++) delete [] spc[i];
}


/**********************************************************/
/*     Memory Allocation                                  */
/**********************************************************/
void cgmAllocateMemory()
{
	try{
		/*** compound nucleus */
		ncl = new Nucleus [MAX_COMPOUND];
		
		/*** calculated results */
		for(int i=0 ; i<4 ; i++){
			spc[i] = new double[MAX_ENERGY_BIN];
		}
	}
	catch(bad_alloc){
		cerr << "ERROR     :memory allocation error";
		exit(-1);
	}
}

