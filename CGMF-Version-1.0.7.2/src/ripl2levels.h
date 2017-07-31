// $Id: ripl2levels.h 5 2011-08-09 02:08:58Z kawano $

//#define LEVELDIRECTORY "levels/"
#define LEVELDIRECTORY "levels-2012/"

/**************************************/
/*      Selection for Nmax            */
/**************************************/
typedef enum {normal=0, extended=1, reassign=2, all=3} MaxLevelCtl;


/**************************************/
/*      ripl2levels.cpp               */
/**************************************/
int     riplReadDiscreteLevels (ZAnumber *, Level *, MaxLevelCtl);
void    riplReadDiscreteLevels (Nucleus *, MaxLevelCtl , int );
