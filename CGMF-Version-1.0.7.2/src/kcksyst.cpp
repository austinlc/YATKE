// $Id: kcksyst.cpp 5 2011-08-09 02:08:58Z kawano $
/******************************************************************************/
/*  kcksyst.cpp                                                               */
/*        level density parameters based on KTUY mass model                   */
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
#include "kcksyst.h"
#include "terminate.h"
#include "config.h"

static string datadir = DATADIR;

/***********************************************************/
/*      Read Level Density Paraters from a File            */
/***********************************************************/
int kckDataRead(ZAnumber *za, LevelDensity *ldp)
{
  ifstream  fp;
  string    d;
  int       z,a;
  
  string str = KCKDATAFILE;
  
  /*** try current directry */
  fp.open(&str[0]);
  
  /*** then system data area */
  if(!fp){
    str = datadir + str;
    fp.open(&str[0]);
  }
  
  if(!fp) cgmTerminateCode("level density parameter file not found");
  
  bool found=false;
  while(getline(fp,str)){
    if(str[0] == '#') continue;
    
    d=str.substr( 0, 5);    z = atoi(&d[0]);
    d=str.substr( 5, 6);    a = atoi(&d[0]);
    ZAnumber za1(z,a);
    
    if((za1.getZ()==za->getZ()) && (za1.getA()==za->getA())){
      
      d=str.substr(11,13);  ldp->pairing_energy = atof(&d[0]);
      d=str.substr(24,13);  ldp->shell_correct  = atof(&d[0]);
      d=str.substr(37,10);  ldp->match_energy   = atof(&d[0]);
      d=str.substr(47,10);  ldp->a              = atof(&d[0]);
      
      if(ldp->match_energy > 0.0){
        d=str.substr(57,10);  ldp->temperature  = atof(&d[0]);
        d=str.substr(67,10);  ldp->E0           = atof(&d[0]);
      }else{
        d=str.substr(77,10);  ldp->temperature  = atof(&d[0]);
        d=str.substr(87,10);  ldp->E0           = atof(&d[0]);
      }
      found = true;
      break;
    }
  }
  fp.close();
  
  
  if(!found){
    ostringstream os;
    os << "level density parameter for Z " << za->getZ() << " - A " << za->getA() << " not found";
    cgmTerminateCode(os.str());
  }
  
  return(0);
}


/***********************************************************/
/*      Systematics for Asymptotic Level Density Parameter */
/***********************************************************/
double kckAsymptoticLevelDensity(double a)
{
  double astar = 0.126181*a + 7.52191e-05*a*a;
  return(astar);
}


/***********************************************************/
/*      Systematics for Spin-Cutoff Parameter              */
/***********************************************************/
double kckSpinCutoff(double a)
{
  double sigma = sqrt(0.0465 *a + 4.296);
  return(sigma);
}


/***********************************************************/
/*      Systematics for T Parameter                        */
/***********************************************************/
double kckTemperature(double a, double dw)
{
  double gamma = 0.1;
  double tsyst = 47.0582 * pow(a,-0.893714) * sqrt(1.0 - gamma*dw);
  return(tsyst);
}


/***********************************************************/
/*      Systematics for Spin-Cutoff Parameter              */
/***********************************************************/
double kckE0(double a, double p, double dw)
{
  double gamma = 0.16;
  double tsyst = kckTemperature(a,dw);
  double esyst = p - gamma*dw + tsyst*(-0.01380*a - 1.379);
  return(esyst);
}
