/****************************/
/*    System Data           */
/****************************/
class System{
 public:
    ZAnumber   target        ;     /* target Z and A numbers           */
    ZAnumber   compound      ;     /* compound nucleus Z and A         */
    Pdata      incident      ;     /* incident particle                */
    double     ex_total      ;     /* total excitation energy          */
    double     cms_energy    ;     /* incident energy                  */
    double     lab_energy    ;     /* incident Lab energy              */

    System(){
      init();
    }

    void init(){
      target.setZA(0,0);
      compound.setZA(0,0);
      incident.particle.setZA(0,0);
      incident.particleID = unknown;
      cms_energy     = 0.0;
      lab_energy     = 0.0;
      ex_total       = 0.0;
    }
};



/**************************************/
/*      excinterface.cpp              */
/**************************************/
void    excitonInterface (const int, const int, const double, double **);


/**************************************/
/*      exciton.cpp                   */
/**************************************/
void    preqExcitonModel(Nucleus *, System *, Pdata *, Transmission **, double **);
