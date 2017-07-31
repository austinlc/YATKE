// $Id: ripl2levels.cpp 5 2011-08-09 02:08:58Z kawano $
/******************************************************************************/
/*  ripl2levels.cpp                                                           */
/*        read in discrete level data from file                               */
/******************************************************************************/

#include <string>
#include <sstream>
#include <ostream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "mt19937ar.h"


using namespace std;

#include "structur.h"
#include "ripl2levels.h"
#include "terminate.h"
#include "config.h"

static void riplSelectCandidate  (string, double *, int *);
static int  riplFixTransition    (int *, double *, double *, int k , Level * );
static int  riplNormalizeBranch  (int, double *, double * );
static int  riplSetNmax          (MaxLevelCtl, int, int, int);
static int  riplFixNmax          (int, Level *);
static int  riplGSonly           (unsigned int, unsigned int, Level *);

static string datadir = DATADIR;


/**********************************************************/
/*      Read Discrete Level Data from RIPL Database       */
/**********************************************************/
int riplReadDiscreteLevels(ZAnumber *za, Level *lev, MaxLevelCtl ml)
{
  ostringstream os;
  ifstream      fp;
  string        str,file,d;
  int           nlevel = 0;
  int           nf[MAX_GAMMA_BRANCH];
  double        pe[MAX_GAMMA_BRANCH],ic[MAX_GAMMA_BRANCH];
  
  os << setw(4) << setfill('0') << za->getZ();

  file = os.str();
  file[0] = 'z';
  str = file + ".dat";
  
  /*** try current directry first */
  fp.open(&str[0]);
  
  /*** then system data area */
  if(!fp){
    str = datadir + LEVELDIRECTORY + file + ".dat";
    fp.open(&str[0]);
  }
  
  if(!fp){
    cgmTerminateCode("discrete level data file not found");
    nlevel = riplGSonly(za->getZ(), za->getA(), lev);
    return(nlevel);
  }
  
  for(int i=0 ; i<MAX_LEVELS ; i++) lev[i].energy = lev[i].spin = 0.0;
  
  bool found=false;
  int  a=0, z=0, nol=0, nc=0, nmax=0;
  while(getline(fp,str)){
    
    /*** search for Z and A entry in the file */
    d=str.substr( 5, 5);  a    = atoi(&d[0]);
    d=str.substr(10, 5);  z    = atoi(&d[0]);
    d=str.substr(15, 5);  nol  = atoi(&d[0]);
    d=str.substr(20, 5);//nog  = atoi(&d[0]);
    d=str.substr(25, 5);  nmax = atoi(&d[0]);
    d=str.substr(30, 5);  nc   = atoi(&d[0]);
    ZAnumber za1(z,a);
    
    if((za1.getZ()==za->getZ()) && (za1.getA()==za->getA())) found = true;
    
    if(nmax >= MAX_LEVELS){
      os.str("");
      os << "discrete levels for Z " << za1.getZ() << " - A " << za1.getA() << " too many (" << nol << ")";
      fp.close();
      cgmTerminateCode(os.str());
    }
    
    /*** set Nmax by the value of MaxLevelCtl */
    nlevel = riplSetNmax(ml, nol, nmax, nc);
        
    /*** for all discrete levels */
    double elev=0.0, s=0.0, thlf=0.0;
    int    p=0, ng=0;
    for(int i=0 ; i<nol ; i++){
      getline(fp,str);
      d=str.substr( 4,10);  elev = atof(&d[0]);
      d=str.substr(15, 5);  s    = atof(&d[0]);
      d=str.substr(20, 3);  p    = atoi(&d[0]);
      d=str.substr(25, 9);  thlf = atof(&d[0]);
      d=str.substr(34, 3);  ng   = atoi(&d[0]);
            
      /*** if no gamma is assigned */
      if(ng==0){
        nf[0] = 0;
        pe[0] = 1.0;
        ic[0] = 0.0;
      }
      /*** for gamma-ray branches */
      else{
        for(int j=0 ; j<ng ; j++){
          getline(fp,str);
          d=str.substr(39, 4);  nf[j] = atoi(&d[0]) -1;
          d=str.substr(66,10);  pe[j] = atof(&d[0]);
          d=str.substr(77,10);  ic[j] = atof(&d[0]);
        }
      }
      
      /*** if found, copy data into lev array */
      if(found && i<nlevel){
        
        if( (s < 0.0) && (i == 0) ) riplSelectCandidate(str, &s, &p);
        if(ng==0){
          ng=riplFixTransition(nf,pe,ic,i,lev);
        }
        ng = riplNormalizeBranch(ng,pe,ic);
        
        lev[i].energy   = elev;
        lev[i].spin     = s;
        lev[i].halflife = thlf;
        lev[i].parity   = p;
        lev[i].ngamma   = ng;
        
        for(int j=0 ; j<ng ; j++){
          lev[i].fstate[j] = nf[j];
          lev[i].branch[j] = pe[j];
          lev[i].gratio[j] = 1.0/(1.0 + ic[j]); // Ig/(Ig + Ie) = 1/(1+ICC)
        }
      }
    }
    if(found) break;
  }
  fp.close();
    
  if(!found){
    if(ml != reassign){
      os.str("");
      os << "discrete level data for Z " << za->getZ() << " - A " << za->getA() << " not found";
      cgmTerminateCode(os.str());
    }
    else{
      nlevel = riplGSonly(za->getZ(), za->getA(), lev);
    }
  }
  else{
    if( (lev[0].spin < 0.0) || lev[0].parity == 0 ){
      if(ml == reassign) nlevel = riplGSonly(za->getZ(), za->getA(), lev);
      else{
        os.str("");
        os << "discrete level data for Z " << za->getZ() << " - A " << za->getA() << " found but spin/parity not assigned";
        cgmTerminateCode(os.str());
      }
    }
  }
  
  if(ml == extended){
    nlevel = riplFixNmax(nlevel,lev);
    if(nlevel == 0){
      os.str("");
      os << "discrete level data for Z " << za->getZ() << " - A " << za->getA() << " not complete";
      cgmTerminateCode(os.str());
    }
  }
    
  return(nlevel);
}

/**********************************************************/
/*      Read Discrete Level Data from RIPL Database       */
/* IS: in this version, I read all the nuclear levels from A to A-nemit */
/**********************************************************/
void riplReadDiscreteLevels(Nucleus *ncl, MaxLevelCtl ml , int nemit )
{
  ostringstream os;
  ifstream      fp;
  string        str,file,d;
  int           nlevel = 0;
  int           nf[MAX_GAMMA_BRANCH];
  double        pe[MAX_GAMMA_BRANCH],ic[MAX_GAMMA_BRANCH];
  int           indx[nemit+1] ;
	
  if(nemit==0){
    ncl[0].ndisc=riplReadDiscreteLevels(&ncl[0].za, ncl[0].lev, ml);
    return ;
  }
	
  os << setw(4) << setfill('0') << ncl[0].za.getZ();
  //  cout << " Z=" << za->getZ() << " A=" << za->getA() << endl ;
  file = os.str();
  file[0] = 'z';
  str = file + ".dat";
	
  /*** try current directry first */
  fp.open(&str[0]);
	
  /*** then system data area */
  if(!fp){
    str = datadir + LEVELDIRECTORY + file + ".dat";
    fp.open(&str[0]);
  }
	
  if(!fp){
    cgmTerminateCode("discrete level data file not found");
  }
	
  int amin = ncl[nemit].za.getA() ;
  int amax = ncl[0].za.getA() ;
  for(int i=0; i<=nemit; ++i)
    indx[ncl[i].za.getA() - amin] = i ;
	
  for(int n=0; n<=nemit; ++n)
    for(int i=0 ; i<MAX_LEVELS ; i++) ncl[n].lev[i].energy = ncl[n].lev[i].spin = 0.0;
	
  //  cout << "searching for za=" << za_final.getZ() << " " << za_final.getA() << endl ;
  bool reading_data = false ;
  int  a=0, z=0, nol=0, nc=0, nmax=0;
  int indx_read ;
  while(getline(fp,str)){
		
    /*** search for Z and A entry in the file */
    d=str.substr( 5, 5);  a    = atoi(&d[0]);
    d=str.substr(10, 5);  z    = atoi(&d[0]);
    d=str.substr(15, 5);  nol  = atoi(&d[0]);
    d=str.substr(20, 5);//nog  = atoi(&d[0]);
    d=str.substr(25, 5);  nmax = atoi(&d[0]);
    d=str.substr(30, 5);  nc   = atoi(&d[0]);
		
    if(a>amax)break ;
    if( amin <= a ){
      indx_read = indx[a-amin] ;
      if(indx_read<0)
				reading_data=false ;
      else
				reading_data = true ;
    }else{
      reading_data = false ;
    }
		
    if(nmax >= MAX_LEVELS){
      ZAnumber za1(z,a);
      os.str("");
      os << "discrete levels for Z " << za1.getZ() << " - A " << za1.getA() << " too many (" << nol << ")";
      fp.close();
      cgmTerminateCode(os.str());
    }
		
    /*** set Nmax by the value of MaxLevelCtl */
    nlevel = riplSetNmax(ml, nol, nmax, nc);
		
    //    cout << "Reading Z=" << z << " A=" << a << endl ;
		
    /*** for all discrete levels */
    double elev=0.0, s=0.0, thlf=0.0;
    int    p=0, ng=0;
    for(int i=0 ; i<nol ; i++){
      getline(fp,str);
      d=str.substr( 4,10);  elev = atof(&d[0]);
      d=str.substr(15, 5);  s    = atof(&d[0]);
      d=str.substr(20, 3);  p    = atoi(&d[0]);
      d=str.substr(25, 9);  thlf = atof(&d[0]);
      d=str.substr(34, 3);  ng   = atoi(&d[0]);
			
      /*** if no gamma is assigned */
      if(ng==0){
        nf[0] = 0;
        pe[0] = 1.0;
        ic[0] = 0.0;
      }
      /*** for gamma-ray branches */
      else{
        for(int j=0 ; j<ng ; j++){
          getline(fp,str);
          d=str.substr(39, 4);  nf[j] = atoi(&d[0]) -1;
          d=str.substr(66,10);  pe[j] = atof(&d[0]);
          d=str.substr(77,10);  ic[j] = atof(&d[0]);
        }
      }
			
      /*** if found, copy data into lev array */
      if(reading_data && i<nlevel){
        if( (s < 0.0) && (i == 0) ) riplSelectCandidate(str, &s, &p);
				if(ng==0){
					ng=riplFixTransition(nf,pe,ic,i,ncl[indx_read].lev);
        }
        ng = riplNormalizeBranch(ng,pe,ic);
				
        ncl[indx_read].lev[i].energy   = elev;
        ncl[indx_read].lev[i].spin     = s;
        ncl[indx_read].lev[i].halflife = thlf;
        ncl[indx_read].lev[i].parity   = p;
        ncl[indx_read].lev[i].ngamma   = ng;
				
        for(int j=0 ; j<ng ; j++){
          ncl[indx_read].lev[i].fstate[j] = nf[j];
          ncl[indx_read].lev[i].branch[j] = pe[j];
          ncl[indx_read].lev[i].gratio[j] = 1.0/(1.0 + ic[j]); // Ig/(Ig + Ie) = 1/(1+ICC)
        }
      }
    }
    if(reading_data){
      ncl[indx_read].ndisc=nlevel ;
      indx[a-amin]=-1 ;
    }
    if(a==amax)break ;
  }
  fp.close();
	
  for(int n=0;n<=nemit;++n){
    if(indx[n] > 0 )
      ncl[n].ndisc = riplGSonly(ncl[n].za.getZ(), ncl[n].za.getA(), ncl[n].lev);
  }
	
  return ;
}


/**********************************************************/
/*      Get Spin and Parity from Candidates               */
/**********************************************************/
void riplSelectCandidate(string str, double *s, int *p)
{
  string csrc = str.substr(46,18);
  char   cdst[] = "                  ";
  
  int i0 = csrc.find("(");
  int i1 = csrc.find_first_of(",)");
  
  if(i0 == (signed int)string::npos) return;
  
  for(int i=i0+1 ; i<i1 ; i++) cdst[i-i0-1] = csrc[i];
  cdst[i1-i0-1] = '\0';
  
  *s = (strchr(cdst,'/') != NULL) ?  0.5 : 1.0;
  *p = (strchr(cdst,'-') != NULL) ? -1   : 1;
  
  i1 = strlen(cdst);
  for(int i=0 ; i<i1 ; i++){
    if( (cdst[i] == '+') || (cdst[i] == '-') || (cdst[i] == '/') ){
      cdst[i] = '\0';
      break;
    }
  }
  *s *= atof(cdst);
}

/**********************************************************/
/*      Fix Transition If No Decay                        */
/**********************************************************/
int riplFixTransition(int nf[], double pe[], double ic[], int k, Level *lev)
{
  int n = 0;
  
  /*** if this is the first excited state */
  if(k == 1){
    nf[0] = 0;
    pe[0] = 1.0;
    ic[0] = 0.0;
    n = 1;
    return(n);
  }
  
  /*** first, look for states to which decay is possible by E1 */
  int    p = lev[k].parity;
  double s = lev[k].spin;
  
  int dm = 100;
  for(int i=0 ; i<k ; i++){
    if( (s == 0.0) && (lev[i].spin == 0.0) ) continue;
    int dj = (int)fabs(lev[i].spin - s );
    if(dj < dm) dm = dj;
    /*** E1, different parity and |J0 - J1| <= 1.0 */
    if( (lev[i].parity != p) && (dj <= 1) ){
      nf[n] = i;
      pe[n] = 1.0;
      ic[n] = 0.0;
      n++;
      if(n == MAX_GAMMA_BRANCH-1) break;
    }
  }
  /*** if no transition still,
   find the closest spin state with the lowest energy */
  if(n == 0){
    for(int i=0 ; i<k ; i++){
      if( (s == 0.0) && (lev[i].spin == 0.0) ) continue;
      int dj = (int)fabs(lev[i].spin - s);
      if(dj == dm){
        nf[0] = i;
        pe[0] = 1.0;
        ic[0] = 0.0;
        n = 1;
        break;
      }
    }
  }
  
  return (n) ;
}

/**********************************************************/
/*      Renormalize Transition Probability                */
/**********************************************************/
int riplNormalizeBranch(int ng, double pe[], double ic[])
{
  double s=0.0;
  
  for(int i=0 ; i<ng ; i++) s += pe[i];
  
  if(s == 0.0){
    if(ng > 1){
      s = 1.0/(double)ng;
      for(int i=0 ; i<ng ; i++){
        pe[i] = s;
        ic[i] = 0.0;
      }
    }
    else{
      pe[0] = 1.0;  // decay 100% to the first line in the branches
      ic[0] = 0.0;
      ng = 1;
    }
  }else{
    s = 1.0/s;
    for(int i=0 ; i<ng ; i++) pe[i] *= s;
  }
  
  return(ng);
}


/**********************************************************/
/*      Determine Nmax from Input                         */
/**********************************************************/
int riplSetNmax(MaxLevelCtl ml, int nol, int nmax, int nc)
{
  int nlevel = 0;
  
  switch(ml){
    case normal:   // Nc is used for the highest level
      nlevel = nc;
      break;
    case extended: // Nmax is used, but upto unknown spin/parity state
    case reassign: // Nmax is used, and unknown level spin/parity will be re-assigned later
      nlevel = nmax;
      break;
    case all:      // include all levels given
      nlevel = nol;
      break;
    default:
      nlevel = nc;
      break;
  }
  
  return(nlevel);
}


/**********************************************************/
/*      Renormalize Transition Probability                */
/**********************************************************/
int riplFixNmax(int nl, Level *lev)
{
  for(int i=0 ; i<nl ; i++){
    /*** spin or parity not assigned */
    if( (lev[i].spin < 0.0) || (lev[i].parity == 0) ){
      nl = i;
      break;
    }
  }
  
  return(nl);
}


/**********************************************************/
/*      Assume Ground State Spin                          */
/**********************************************************/
int riplGSonly(unsigned int z, unsigned int a, Level *lev)
{
  unsigned int n = a -z;
  
  /*** even-even */
  if(((z%2) == 0) && ((n%2) == 0) ){
    lev[0].spin   = 0.0;
    lev[0].parity = 1;
  }
  /*** odd-odd, 1+ assumed */
  else if(((z%2) != 0) && ((n%2) != 0) ){
    lev[0].spin   = 1.0;
    lev[0].parity = 1;
  }
  /*** even-odd, 1/2+ assumed */
  else{
    lev[0].spin   = 0.5;
    lev[0].parity = 1;
  }
  
  return(1);
}
