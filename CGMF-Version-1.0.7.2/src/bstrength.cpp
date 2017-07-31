// $Id: bstrength.cpp 5 2011-08-09 02:08:58Z kawano $
/******************************************************************************/
/*  bstrength.cpp                                                             */
/*        beta-decay strength from FRDM calculations                          */
/******************************************************************************/

#include <string>
#include <sstream>
#include <ostream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "structur.h"
#include "bstrength.h"
#include "config.h"

static string datadir = DATADIR;

/***********************************************************/
/*      Read Beta Decay Strength Distribution from a File  */
/***********************************************************/
int betaDataRead(ZAnumber *za, Beta *beta)
{
  ZAnumber      za1;
  ostringstream os;
  ifstream      fp;
  string        str,file,d1,d2;

  os << setw(3) << setfill('0') << za->getZ() << setw(3) << setfill('0') << za->getA();
  file = os.str();
  str = "strint" + file + ".dat";

  /*** try current directry first */
  fp.open(&str[0]);

  /*** then system data area */
  if(!fp){
    str = datadir + "/" + BETADIRECTORY + str;
    fp.open(&str[0]);
  }

  if(!fp){
    beta->nstate=0;
    return(beta->nstate);
  }

  getline(fp,str);
  getline(fp,str);
  getline(fp,str);
//d1 = str.substr(42,10); double qbeta = atof(&d1[0]);

  int k = 0;
  while(getline(fp,str)){
    d1 = str.substr( 0,10);
    d2 = str.substr(10,10);

    beta->br[k].setVal(atof(&d1[0]),atof(&d2[0]));
    k++;
    if(k >= MAX_BETA_BRANCH){
      k--;
      break;
    }
  }
  fp.close();

  beta->nstate = k;


  /*** renormalize the probability */
  double sum = 0.0;
  for(int j=0 ; j<beta->nstate ; j++) sum += beta->br[j].getR();
  if(sum==0.0) beta->nstate = 0;

  for(int j=0 ; j<beta->nstate ; j++) beta->br[j].scaleR(1.0/sum);
/*
  cout << "# GT Sum " << sum << "  "<< beta->nstate << endl;
  for(int i=0 ; i<beta->nstate ; i++){
    cout << "#GT "
         << setprecision(4) << setiosflags(ios::scientific) << setw(11)
         << beta->br[i].getE() << setw(11) << beta->br[i].getR() << endl;
  }
*/
  return(beta->nstate);
}
