#define OUT_SPEC       0x0001    //  1
#define OUT_INIPOP     0x0002    //  2
#define OUT_HIST       0x0004    //  4
#define OUT_BROADENED  0x0008    //  8
#define OUT_INDIVSPEC  0x0010    // 16
//#define              0x0020    // 32
//#define             0x0040
//#define             0x0080
#define CALC_MC     0x0100
#define CALC_BETA   0x0200
#define CALC_TRAN   0x0400

typedef enum {SINGLE, LEVDEN, BETA, TRANSMISSION} InitialPop;

#define DEFAULT_CAL (OUT_SPEC)


class Calculation{ 
 public:
    bool print_history    ;
    bool print_spectrum   ;
    bool print_init_pop   ;
    bool print_broadened  ;
    bool print_indivspec  ;
    bool calc_montecarlo  ;
    bool calc_betadecay   ;
    bool calc_entrance    ;
    InitialPop init_pop   ;
};

extern Calculation ctl;

extern string fileExt;
extern double timeGate;
