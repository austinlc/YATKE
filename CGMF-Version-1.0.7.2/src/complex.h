// $Id: complex.h 5 2011-08-09 02:08:58Z kawano $
/*
   complex.h : 
        definition of complex number manipulation
        basically this is not needed in C++
 */

/****************************/
/*      complex value       */
/****************************/
class Complex{ 
 public:
    double real;
    double imag;

    Complex(){
      real = 0.0;
      imag = 0.0;
    }
    Complex(double r, double i){
      real = r;
      imag = i;
    }
};


/**************************************/
/*      etc.cpp                       */
/**************************************/
Complex rational              (double, double, double, double);
Complex compsqrt              (Complex *);
Complex comppow               (Complex *, double);
double  absolute              (Complex *);

