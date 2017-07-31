// $Id: masstable.h 5 2011-08-09 02:08:58Z kawano $
/*
   masstable.h : 
        define mass excess data file location
        prototype of function to read mass table
 */

class MassExcess{
 public:
  unsigned int za;    // Z*1000 + A
  float        mass;  // mass excess
};

/**************************************/
/*      masstable.cpp                 */
/**************************************/
int     read_masstable        (int, int, double *);
double  mass_excess           (int, int);

