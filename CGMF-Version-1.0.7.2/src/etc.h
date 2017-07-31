// $Id: etc.h 5 2011-08-09 02:08:58Z kawano $
/*
   etc.h :
        prototype of miscellaneous functions
        complex variable manipulation carried over from old C code
*/


#ifndef __COMPLEX_H__
#define __COMPLEX_H__
#include "complex.h"
#endif

/**************************************/
/*      etc.cpp                       */
/**************************************/
double  jvol_to_dep           (double, double, double, double);
double  jsrf_to_dep           (double, double, double, double);
double  dep_to_jvol           (double, double, double, double);
double  dep_to_jsrf           (double, double, double, double);

#ifndef HAVE_MINMAX
int     min                   (int   ,int   );
int     max                   (int   ,int   );
#endif

double  cfmin                 (double,double);
double  cfmax                 (double,double);

double  gaussian_weight       (double, double, double);
double  laguerre              (int,    double, double);
double  gam                   (double);
double  loggamma              (double);
double  legendre              (int, double);
double  legendre1             (int, double);
double  assocLegendrePol      (int,int,double);
double  bessi2                (int, double);
double  bessk2                (int, double);

