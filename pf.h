#ifndef _pf_H
#define _pf_H 1

#include "TGraph.h"
#include "TRandom3.h"

Double_t func(Double_t x, Double_t bkg, Double_t Speak, Double_t pos, Double_t res, Double_t &gaus);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
Int_t generate_sample(TRandom3 *xr, Double_t time, Double_t invT12);
Int_t initfit(Int_t in_debug,const char *fn, const char *fn2, Int_t sys_bits,const char *fncuthist);
Int_t scanfit(Double_t *best_xinv, Double_t *limit90, Double_t *like0,
	      Double_t **xi, Double_t **fcns, Double_t **fcns_tlimit, Int_t* Nstep);
Double_t profile1(Double_t xinv_true, Double_t &best_xinv);
Double_t profile(Double_t xinv_true, Double_t &best_xinv);
Int_t   read_data(const char *fn);
int fixsystematic();
int relsystematic();
int nsignal();
int get_n1sig();
int get_n2sig();
double get_distQbb();
double get_binwidth();
double *fcndebug();


#define QBB 2039.061
#define DeltaQBB 0.2
#define  NA (6.022141*0.693147/7.56)
#define GRAD    0  // calculate FCN gradient to speed up
#define CL90CUT 2.7055
#define MAXEVT 500

struct detpartition{
  char name[10];      // detector name
  double resolution;  // energy resolution sigma !!
  double  resUnc;     // 
  double efficiency;
  double effUnc;
  double exposure;
  double shiftE;
  double shiftEUnc;
  int    fit_bindex;
  double bkgindex;
  char   bkgname[20];  // background set name
  char sysname[20];
  int sysnumber;
  double Eevt[MAXEVT];
};

#define DETMAX 143
struct partition{
  int t1,t2;
  int ndet;
  double exposureP;
  double avgeff;
  int nevt_bkg[DETMAX];
  int nevt_sig[DETMAX];
  struct detpartition dp[DETMAX];
};

struct bkginput{
  int num;       // increasing number of bkg set
  double value;  // background index
  double exposure;
  double nevt;
};

#endif

