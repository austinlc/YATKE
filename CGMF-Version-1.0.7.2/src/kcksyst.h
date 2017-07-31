// $Id: kcksyst.h 5 2011-08-09 02:08:58Z kawano $
/*
   kcksys.h : 
        define level density parameter file location
        prototype of function to read level density parameters
 */

#define KCKDATAFILE "kcksyst.dat"

/**************************************/
/*      kcksyst.cpp                   */
/**************************************/
int     kckDataRead                     (ZAnumber *, LevelDensity *);
double  kckAsymptoticLevelDensity       (double);
double  kckSpinCutoff                   (double);
double  kckTemperature                  (double, double);
double  kckE0                           (double, double, double);
