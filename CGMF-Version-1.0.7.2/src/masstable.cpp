// $Id: masstable.cpp 5 2011-08-09 02:08:58Z kawano $
/******************************************************************************/
/*  masstable.cpp                                                             */
/*        find nuclear mass from Mass Table                                   */
/******************************************************************************/

#include <string>
#include <sstream>

using namespace std;

#include "masstable.h"
#include "terminate.h"

//#include "masstable_aw95.h"      // Audi Wapstra 1995 table
//#include "masstable_ripl2.h"     // AW95 + FRDM95 from RIPL2
//#include "masstable_ripl3.h"     // AW03 + FRDM95 from RIPL3
//#include "masstable_audi2011.h"    // AW11 + FRDM95 from RIPL3
#include "masstable_audi2012_frdm2012.h" // Audi-Wapstra 2012 + FRDM-2012 

double mass_excess(int z, int a)
{
  double    mx  = 0.0;
  unsigned int za = z*1000+a;

  bool found = false;
  for(int i=0 ; i<nMassTable ; i++){
    if(MassTable[i].za == za){
      found = true;
      mx = MassTable[i].mass;
      break;
    }
  }

  if(!found){
    ostringstream os;
    os << "mass data for Z " << z << " - A " << a << " not found";
    cgmTerminateCode(os.str());
  }

  return(mx);
}

