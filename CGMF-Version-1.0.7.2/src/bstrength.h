// $Id: bstrength.h 5 2011-08-09 02:08:58Z kawano $
/*
   bstrength.h : 
        define beta strength file location
        prototype of function to read beta decay parameters
 */

#define BETADIRECTORY   "bstrength/"
#define ENSDFDIRECTORY  "ENSDF/"

static const int MAX_BETA_BRANCH = 1000; // maximum beta branches

/****************************/
/*   Beta Decay Branching   */
/****************************/
class Branch{
 private:
  double energy;
  double ratio;
 public:
  Branch(){
    energy = 0.0;
    ratio  = 0.0;
  }
  double getE(void){
    return(energy);
  }
  double getR(void){
    return(ratio);
  }
  void setE(double e){
    energy = e;
  }
  void setR(double b){
    ratio  = b;
  }
  void scaleR(double x){
    ratio  *= x;
  }
  void setVal(double e, double b){
    energy = e;
    ratio  = b;
  }
};


/****************************/
/*   Beta Decay Parameter   */
/****************************/
class Beta{
 public:
    int     nstate           ;     /* nubmer of final states            */
    Branch  *br              ;     /* branching ratio                   */

    Beta(){
      nstate = 0;
      br     = NULL;
    }
};



/**************************************/
/*      bstrength.cpp                 */
/**************************************/
int     betaDataRead                    (ZAnumber *, Beta *);

/**************************************/
/*      ensdf.cpp                     */
/**************************************/
int     betaENSDFRead                   (ZAnumber *, Beta *);
