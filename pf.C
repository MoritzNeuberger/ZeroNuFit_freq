// profile likelihood fit 
// 
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMinuit.h"
#include "TGraph.h"
#include <string>

#include "pf.h"
#include <map>
#include <string>

#define MAXSET 100
struct partition part[MAXSET];

// global parameters
Int_t     Nsets;     // number of partitions 
Int_t     Nbkgsets;  // number of background parameters

Double_t  TcutBinWidth;
Int_t     TcutNbins;
Double_t  *TcutValue;  // 90% CL cut threshold as fct of 1/T
TMinuit   *gMinuit;
Int_t     debug;
Int_t     sys_sigma,sys_pos,sys_eff;  // bits for systematic errors float yes/no

std::map<std::string, struct bkginput> bkgnames;
std::map<std::string, int> syspar;

Double_t  last_fcn[14];  // debug info only
int n1sig; // number of events within 1 sigma of Qbb, for debugging only
int n2sig; //  number of events within 2 sigma of Qbb, for debugging only

int ksys=0;  //number of systematic error groups

double    distQbb;  //  minimal distance event energy to Qbb, for debugging only

double  binwidth = 0.00075;  // global bin width  in 1/T in units of 1E-25/yr

int get_n1sig() {return n1sig;}
int get_n2sig() {return n2sig;}
double get_distQbb() {return distQbb;}
double get_binwidth(){return binwidth;}
//______________________________________________________________________________
Double_t func(Double_t x, Double_t bkg, Double_t Speak, Double_t pos, Double_t res, Double_t &gaus)
{
  // fit function: gauss + flat bkg
  gaus =   exp(-0.5*(x-pos)*(x-pos)/res/res)/res/2.5066283;
  Double_t value = bkg + Speak*gaus;

 return value;

}

//------------------------------------------------------------------------------

Double_t calculate_range(Double_t roi[][2], int size) {
    Double_t output = 0;
    for (int i = 0; i < size; i++) {
    	output += (roi[i][1] - roi[i][0]);
    }

    return output;
}

Double_t GERDA_ROI[3][2] = {
    {1930, 2099},
    {2109, 2114},
    {2124, 2190}
};

Double_t GERDA_ROI_length_minus_gammas = calculate_range(GERDA_ROI,3);
Double_t GERDA_ROI_length_plus_gammas = GERDA_ROI[2][1] - GERDA_ROI[0][0];

// Define the MJD_ROI array
Double_t MJD_ROI[3][2] = {
    {1950, 2099},
    {2109, 2114},
    {2124, 2350}
};

Double_t MJD_ROI_length_minus_gammas = calculate_range(MJD_ROI,3);
Double_t MJD_ROI_length_plus_gammas = MJD_ROI[2][1] - MJD_ROI[0][0];


//------------------------------------------------------------------------------

double* fcndebug(){
  return last_fcn;
}
//-------------------------------------------------------------------------------
int generate_sample(  TRandom3 *xr, Double_t tim, Double_t invT12 )
// throw random toy data samples and initialize fit parameters
//   event numbers are poisson distributed according to background index and T12
// tim = time for exposure (no longer used)
// invT12 = signal strength
{
  Double_t *E;
  Double_t arglist[2];
  Int_t nbkg,i,ierflg,iset,ipar;
  char snam[30]={"generate_sample"};
  Int_t nmax = 1000;
  std::string s;
  int  nsigtot,nbkgtot;  
  double tmp=0;
  double tmp2=0;
  double tmp3=0;
  n1sig=0;  // debug info only
  n2sig=0;  // debug info only
  
  if(debug) printf("%s: invT12 input = %f \n",snam,invT12);

  
  // initialize  invT12 for fit, parameter = 1
  arglist[0]=1;
  arglist[1]=0.1;
  gMinuit->mnexcm("SET PAR", arglist ,2,ierflg);
  if(ierflg) printf("%s:Error setting invT12\n",snam);

  typedef std::map<std::string,struct bkginput>::iterator iter;
  for(iter it = bkgnames.begin(); it != bkgnames.end(); ++it)
    it->second.nevt = 0;
  nsigtot=nbkgtot=0;
  distQbb = 1.e5; // distance of energy from Qbb in units of resolution
  

    for(iset=0; iset<Nsets; iset++)
      for(int idet=0; idet<DETMAX; idet++)
        for (Int_t ievt=0;ievt<MAXEVT; ievt++)
          part[iset].dp[idet].Eevt[ievt] = 0;

  // loop over all partitions and detectors
  for(iset=0; iset<Nsets; iset++) {

    for(int idet=0; idet<DETMAX; idet++) {

      if(part[iset].dp[idet].exposure < 1e-5) continue;
      
      // calc expected  background and signal counts, using background index and 1/T
      //Double_t bkgcnt = part[iset].dp[idet].bkgindex * part[iset].dp[idet].exposure *260.;
      Double_t bkgcnt = part[iset].dp[idet].bkgindex * part[iset].dp[idet].exposure * GERDA_ROI_length_plus_gammas;  // smaller window due to MJD/GERDA overlapp
      //printf("GERDA_ROI_length_plus_gammas: %f",GERDA_ROI_length_plus_gammas);
      if(strcmp(part[iset].dp[idet].name, "ppc") == 0){
        bkgcnt = part[iset].dp[idet].bkgindex * part[iset].dp[idet].exposure * MJD_ROI_length_plus_gammas;
	}
      Double_t Speak= part[iset].dp[idet].exposure*NA*invT12*part[iset].dp[idet].efficiency;

      int nbkg = xr->Poisson(bkgcnt);
      int nsignal = xr->Poisson(Speak);

      if(nbkg+nsignal >= MAXEVT) {printf("%s: buffer for Eevt too small, nbkg=%d nsignal=%d dataset %d\n",snam,nbkg,nsignal,iset); return -2;}
      part[iset].nevt_bkg[idet]=nbkg;    // number of bkg events in detector number idet
      part[iset].nevt_sig[idet]=nsignal; // number of sig events in detector number idet
      nsigtot += nsignal;
      nbkgtot += nbkg;

      s = std::string(part[iset].dp[idet].bkgname);
      iter it = bkgnames.find(s);
      if( it == bkgnames.end()) {
         printf("%s: did not find name of bkg_fit_param %s\n",snam,part[iset].dp[idet].bkgname); return -2;
       }
       it->second.nevt += nbkg;
      
      if(debug>1) printf("%s: data set %d  bkgcnt=%f Nevt %d Speak %f nsignal %d \n",
		       snam,iset,bkgcnt,nbkg,Speak,nsignal);

      for (Int_t ievt=0;ievt<nbkg; ievt++) {
	//	part[iset].dp[idet].Eevt[ievt]  = xr->Uniform(1930.,2190.); // GERDA background window

	part[iset].dp[idet].Eevt[ievt]  = xr->Uniform(GERDA_ROI[0][0],GERDA_ROI[2][1]); 
  if(strcmp(part[iset].dp[idet].name, "ppc") == 0){
	  part[iset].dp[idet].Eevt[ievt]  = xr->Uniform(MJD_ROI[0][0],MJD_ROI[2][1]);  // MJD background window 1950-2350 keV --> overlapp interval

  }
          //printf("E: %f - [%f,%f] \n",part[iset].dp[idet].Eevt[ievt],GERDA_ROI[0][0],GERDA_ROI[2][1]);

	double x= fabs(part[iset].dp[idet].Eevt[ievt] - QBB) / part[iset].dp[idet].resolution;
	if(x<1.) n1sig++;  // for debugging only
	if(x<2.) n2sig++;  // for debugging only
	if(x<distQbb) distQbb=x;  // for debugging only: event energy closest to Qbb
      }	

      for (int ievt=nbkg;ievt<nsignal+nbkg; ievt++) {
	double shift = part[iset].dp[idet].shiftE;  // shift of peak position from SEP
	part[iset].dp[idet].Eevt[ievt]= xr->Gaus(QBB-shift,part[iset].dp[idet].resolution);
	double x= fabs(part[iset].dp[idet].Eevt[ievt] - QBB) / part[iset].dp[idet].resolution;
	if(x<1.) n1sig++;
	if(x<2.) n2sig++;
	if(x<distQbb) distQbb=x;
      }

      if(debug>2) printf("%s: data set %d det %d: bkg %d  signal %d \n",snam,iset,idet,nbkg,
			 nsignal);
    } // idet
    
  } // iset


  if(debug)printf("%s: total signal %d  total bkg %d  E<1s of Qbb: %d ,2s: %d  closest %f\n",
		     snam,nsigtot,nbkgtot,n1sig,n2sig,distQbb);


  // initialize Miniut bkg fit parameter
  
  for(iter it = bkgnames.begin(); it != bkgnames.end(); it++) {
    arglist[0]=it->second.num+2;      // parameter number, starting from 1
    //arglist[1]=it->second.nevt/240./it->second.exposure;

    arglist[1]=it->second.nevt/GERDA_ROI_length_plus_gammas/it->second.exposure;  
    gMinuit->mnexcm("SET PAR", arglist ,2,ierflg);
    if(ierflg) {printf("%s: Error setting bkg for data set %d\n",snam,iset); return -1;}
    if(debug) printf("%s: bkg set %s  tot evt %f   exposure %f\n",snam,it->first.c_str(),it->second.nevt,
		     it->second.exposure);
  }

  return nsigtot+nbkgtot;

}
//-----------------------------------------------------------------------
int nsignal() {
  // return number of signal events  in gaussian
  int sum=0;
  for(int iset=0; iset<Nsets;iset++)
    for(int idet=0; idet<DETMAX; idet++) {
      if(part[iset].dp[idet].exposure<1e-4) continue;
      sum += part[iset].nevt_sig[idet];
    }
  return sum;
}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   Int_t i;

//calculate log likelihood = sum_bins (log likelih.)  with poisson statistics
// maximize product_bins  probability  or minimize    -2 * sum_bins log likel.
   Double_t chisq = 0;
   Double_t Speak, Bck, Nexpected;
   Double_t delta, mu, gaus;
   Double_t mean;
   Int_t ievt,iset,idet;
   Double_t *E;
   Int_t     Nevt_range;
   int evtbkgs[100];
   int sysres, sysshift, syseff;
   
   f=0;
//    if(GRAD) gin[0]=0;


   for(int i=0;i<100;i++) evtbkgs[i]=0;
   
   // loop over all  partitions / detectors
   for(iset=0;iset<Nsets;iset++) {
     for(idet=0; idet<DETMAX; idet++) {
     
       if(part[iset].dp[idet].exposure < 1.e-4) continue;

       int index = part[iset].dp[idet].fit_bindex+1; // index of minuit parameter (index=0 is 1/T)
       sysres   = 1 + Nbkgsets + 0 + part[iset].dp[idet].sysnumber*3;  //  fit parameter number fit for systematic error resolution
       sysshift = 1 + Nbkgsets + 1 + part[iset].dp[idet].sysnumber*3;  //  systematic error energy shift
       syseff   = 1 + Nbkgsets + 2 + part[iset].dp[idet].sysnumber*3;  //  systematic error efficiency

       
     // correct resolution, peak position and efficiency for systematic errors,
       // fit parameter is a scale factor around 1 that multiplies the  resolution uncertainty
       Double_t res = part[iset].dp[idet].resolution + part[iset].dp[idet].resUnc*par[sysres];
       // instead of shift event energy by +shiftE  I move peak position by -shiftE
       Double_t position = QBB - part[iset].dp[idet].shiftE - part[iset].dp[idet].shiftEUnc * par[sysshift];
       Double_t eff = part[iset].dp[idet].efficiency + part[iset].dp[idet].effUnc*par[syseff];
       
       Speak = part[iset].dp[idet].exposure *NA*par[0]* eff;  // expected count of peak events
       //Bck = 240.*par[index]*part[iset].dp[idet].exposure;    // count of bkg events
       Bck = GERDA_ROI_length_minus_gammas * par[index]*part[iset].dp[idet].exposure;    // count of bkg events, smaller fit window from 1950-2190 minus 2 gamma lines
      if(strcmp(part[iset].dp[idet].name, "ppc") == 0){
          Bck=MJD_ROI_length_minus_gammas * par[index]*part[iset].dp[idet].exposure;  // smaller energy window due to MJD overlapp
      }      

       Nexpected = Speak + Bck;

       if (!Nexpected) {
	 f = 1.e20; //std::numeric_limits<double>::max();
	 return;
       }


//      if(GRAD) gin[par_off]=0;
       mean = QBB;
       chisq=0;

     
     // loop over all events of this data set / partition
       Nevt_range = 0;
       for(ievt=0; ievt<part[iset].nevt_bkg[idet] + part[iset].nevt_sig[idet];  ievt++) {
       // filter data in a specific range
	 
      if(strcmp(part[iset].dp[idet].name, "ppc") == 0){
        if( part[iset].dp[idet].Eevt[ievt] <  MJD_ROI[0][0] || part[iset].dp[idet].Eevt[ievt] >  MJD_ROI[2][1]  ||    // reduced window, overlapp MJD
        ( part[iset].dp[idet].Eevt[ievt] >= MJD_ROI[0][1] && part[iset].dp[idet].Eevt[ievt] <= MJD_ROI[1][0]) ||
        ( part[iset].dp[idet].Eevt[ievt] >= MJD_ROI[1][1] && part[iset].dp[idet].Eevt[ievt] <= MJD_ROI[2][0]) )  {
          //printf(" removed set %d  idet %d   energy %f \n",iset,idet, part[iset].dp[idet].Eevt[ievt]);
          continue;
	      }
      } 
      else{
        if( part[iset].dp[idet].Eevt[ievt] <  GERDA_ROI[0][0] || part[iset].dp[idet].Eevt[ievt] >  GERDA_ROI[2][1]  ||    // reduced window, overlapp MJD
        ( part[iset].dp[idet].Eevt[ievt] >= GERDA_ROI[0][1] && part[iset].dp[idet].Eevt[ievt] <= GERDA_ROI[1][0]) ||
        ( part[iset].dp[idet].Eevt[ievt] >= GERDA_ROI[1][1] && part[iset].dp[idet].Eevt[ievt] <= GERDA_ROI[2][0]) )  {
          //printf(" removed set %d  idet %d   energy %f \n",iset,idet, part[iset].dp[idet].Eevt[ievt]);
          continue;
	      }
      }

	 Nevt_range++;  
	 evtbkgs[index]++;
	 
	 // normalized PDF = func
   Double_t pdf = 0;
      if(strcmp(part[iset].dp[idet].name, "ppc") == 0){
        pdf  = func(part[iset].dp[idet].Eevt[ievt], Bck/MJD_ROI_length_minus_gammas , Speak, position, res, gaus) / Nexpected;
      } 
      else{
        pdf  = func(part[iset].dp[idet].Eevt[ievt], Bck/GERDA_ROI_length_minus_gammas , Speak, position, res, gaus) / Nexpected;
      }
	  
       
	 chisq += log(pdf);  // Prod_bins PDF --> sum_bins ln(PDF)
       

	 if(part[iset].dp[idet].Eevt[ievt]>0 && part[iset].dp[idet].Eevt[ievt]<4044.&&debug>2) {
      if(strcmp(part[iset].dp[idet].name, "ppc") == 0){
	   printf("FCN: i=%d E=%f LogLik=%e sum=%f var %f %f %f\n",ievt,part[iset].dp[idet].Eevt[ievt],pdf,chisq,gaus,Speak,Bck/MJD_ROI_length_minus_gammas);
      } 
      else{
	   printf("FCN: i=%d E=%f LogLik=%e sum=%f var %f %f %f\n",ievt,part[iset].dp[idet].Eevt[ievt],pdf,chisq,gaus,Speak,Bck/GERDA_ROI_length_minus_gammas);
      }
       }
//        if(GRAD) {  // return gradient of log chi2
// 	 gin[par_off] -= 2./pdf;
// 	 gin[0] -= 2.*gaus*Exposure[iset]*NA*par[par_off*(iset+1)]/pdf;
//        }
       
       } // end event loop

       
       last_fcn[0]=-2.*chisq;  // debug only
     
       f += -2.*chisq;

       last_fcn[1] = f-last_fcn[0];  // debug only
   
       // constrained fit for total number of events
       // ln(poisson (nevt,mu) ) (forget nevt! term)
       if(debug>3) printf("FCN: evt loop data sets %d: %f ",iset,f);

       // poissonian for mu!=0 and mu==0
       if (Nevt_range)
	 chisq = Nevt_range*log(Nexpected ) - Nexpected;
       else 
	 chisq = (- Nexpected);

       last_fcn[2]=-2.*chisq;  // debug only
       last_fcn[3]=par[0];      // debug only
       last_fcn[4]=Nexpected;   // debug only
       
       f += -2.*chisq;
       //      if(GRAD) {
       //        gin[par_off] -= 2*Nevt_range*240./mu - 2.*240.;
       //        gin[0] -= 2*Nevt_range*Exposure[iset]*NA*par[par_off*(iset+1)]/mu 
       //   	       - 2.*Exposure[iset]*NA*par[par_off*(iset+1)]; 
       //      }
		 
       
       if(debug>2)
	 printf("FCN: data set %d Speak %f  bkg %f nevt %d  Out %f \n",
		iset,Speak,Bck,Nevt_range,f);
       // for the next experiment/data set
       
     } // loop over detectors
   } // loop over sets/partitions


   // add pull terms for systematic uncertainties, one for each "systematic index" set

   Double_t chisq1,chisq2,chisq3;
   chisq1=chisq2=chisq3=0;
   
   for(int i=0; i<ksys; i++){
     sysres   = 1 + Nbkgsets + 0 + i*3;  // fit index
     sysshift = 1 + Nbkgsets + 1 + i*3;
     syseff   = 1 + Nbkgsets + 2 + i*3;
  
     chisq1 += -0.5*par[sysres]*par[sysres];
		
     chisq2 += -0.5*par[sysshift]*par[sysshift];
       
     chisq3 += -0.5*par[syseff]*par[syseff];
   }
   
   if(debug>2) printf("FCN: Eres %g Qbb %g Eff %g chi2 %g %g %g\n",par[sysres],par[sysshift],par[syseff],chisq1,chisq2,chisq3);
       
   chisq=0;
   if(sys_sigma) chisq+=chisq1;
   if(sys_pos)   chisq+=chisq2;
   if(sys_eff)   chisq+=chisq3;
   f += -2.*chisq;   // final value returned to MINUIT

   if(debug>2)printf("evts %d %d %d %d %d %d \n",evtbkgs[0],evtbkgs[1],evtbkgs[2],evtbkgs[3],evtbkgs[4],evtbkgs[5]);
   
}
//----------------------------------------------------------------------------
 Int_t initfit(Int_t in_debug, const char *fn, const char *fn2, Int_t sys_bits, const char *fncuthist)
{  
// initialize Minuit, read in all information needed 
// fn:  file name with  information for data sets / partitions
// fn2:  file with information on  bkg index 
// sys_bits=bit pattern for systematic uncertainty
//     bit0= sigma of gauss, bit1=position gaus, bit2=efficiency
// fncuthist: root file with histogram for 90% CL cut threshold
  
  Int_t i;
  char line[200];
  
  char snam[50]="initfit";  
  typedef std::map<std::string, int>::iterator itersys;

  sys_sigma= (0x1 & sys_bits);
  sys_pos  = (0x2 & sys_bits);
  sys_eff  = (0x4 & sys_bits); 

  debug=in_debug;
   gMinuit = new TMinuit(34);  //initialize TMinuit 
   gMinuit->SetFCN(fcn);       // function name =  fcn

   Double_t arglist[10];
   Int_t ierflg = 0;
   Int_t ipar=2;
   char  nam[40];
   
   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
   if(ierflg) return ierflg;

   // ------------------------------------------------------------------------------
   // read in file with name/value of background parameters in the fit

   FILE *f = fopen(fn2,"r");
   if( f==0 ) {printf("error opening file background info"); return 2;}

   Nbkgsets=0;

   double value;
   struct bkginput k = {0};
   while( fscanf(f,"%s %lf",nam,&value)==2 ) {
     k.num = Nbkgsets;
     k.value = value; // bkg index
     bkgnames.insert(std::make_pair(nam,k));     
     Nbkgsets++;
   }

   fclose(f);

   
   // ------------------------------------------------------------------------------
   // read in information for different data sets / partitions
   f = fopen(fn,"r");
   if( f==0 ) {printf("error opening file data set info"); return 2;}

   Nsets=-1;
   while( fgets(line,200,f) ) {

     // new partition, starting with time stamp info for first, last
     if(debug>2) printf("%s \n",line);
     if(line[0] == 't') {
       // new partition, read in time stamp interval
       int t1,t2;
       char ch[10];
       Nsets++;
       int i= sscanf(line,"%s %d %d",ch,&t1,&t2);
       if(i!=3) {printf(" error reading time stamps for new partition %d \n",Nsets); return 2;}
       part[Nsets].t1 = t1;
       part[Nsets].t2 = t2;
       part[Nsets].exposureP=part[Nsets].ndet=0;
       for(int i=0; i<DETMAX; i++) part[Nsets].nevt_bkg[i]=part[Nsets].nevt_sig[i]=part[Nsets].dp[i].exposure=0;
     }
     else {
       // loop over all detectors in a partition
       int detnum;
       char detnam[20];
       double fwhm;
       double fwhm_err;
       double eff;
       double eff_err;
       double exposure;
       double shift;
       double shift_err;
       char bkgname[20];
       char nnn[20] = {"bkgp2"};
       char sysnam[20];
       
       int i=sscanf(line,"%d %s %lf %lf %lf %lf %lf %lf %lf %s %s",&detnum,detnam,&fwhm,&fwhm_err,
		    &eff,&eff_err,
		    &exposure,&shift,&shift_err,bkgname,sysnam);
       if(i!=11) {printf("error reading det info for partition\n"); return 2;}
       if(detnum<0 || detnum>=DETMAX)
	 {printf("error reading det info: illegal det num %d\n",detnum); return 2;}
       
       strcpy(part[Nsets].dp[detnum].name,detnam);
       part[Nsets].dp[detnum].resolution = fwhm/2.35;
       part[Nsets].dp[detnum].resUnc     = fwhm_err/2.35;
       part[Nsets].dp[detnum].efficiency = eff;
       part[Nsets].dp[detnum].effUnc     = eff_err;
       part[Nsets].dp[detnum].exposure   = exposure;
       part[Nsets].dp[detnum].shiftE     = shift;
       part[Nsets].dp[detnum].shiftEUnc  = shift_err;
       strcpy(part[Nsets].dp[detnum].sysname,sysnam);
       part[Nsets].exposureP            += exposure;

       // now look up fit parameter for background index belonging to name of bkgname
       typedef std::map<std::string,struct bkginput>::iterator iter;
       iter it = bkgnames.find(std::string(bkgname));
       if( it == bkgnames.end()) {
         printf(" did not find name of bkg_fit_param %s\n",bkgname); return 2;
       }

       it->second.exposure += exposure;
       part[Nsets].dp[detnum].fit_bindex = it->second.num;  // minuit fit parameter number
       part[Nsets].dp[detnum].bkgindex   = it->second.value;  // bkg index of data set
       strcpy(part[Nsets].dp[detnum].bkgname, bkgname);

       // systematic error
       itersys its = syspar.find(std::string(sysnam));
       if(its == syspar.end()){
	 syspar.insert(std::make_pair(sysnam,ksys));
	 part[Nsets].dp[detnum].sysnumber = ksys;
	 ksys++;
       }
       else {
	 part[Nsets].dp[detnum].sysnumber = its->second;
       }
       if(debug>1) printf("%s %s %d \n",part[Nsets].dp[detnum].name,part[Nsets].dp[detnum].sysname, part[Nsets].dp[detnum].sysnumber);
       
       //       
       part[Nsets].ndet++;    // number of detectors in partition set
     }
   }

   fclose(f);
   Nsets++;

   // ---------------------------------------------------------------------------
   // Set starting values and step sizes for parameters
   // 0=1/T12, 1=bkg data set 1, 2=bkg data set 2, ...

   gMinuit->mnparm(0, "invT12", 1.e-1, 0.02, -0.,50.,ierflg); // unit 1E-25 1/yr
   if(ierflg) return ierflg;

   typedef std::map<std::string,struct bkginput>::const_iterator iter;

   for(iter it = bkgnames.begin(); it != bkgnames.end(); ++it){
     ipar = 1+it->second.num;
     gMinuit->mnparm(ipar, it->first, it->second.value, 0.2*it->second.value,0,5*it->second.value,ierflg);
     if(ierflg) return ierflg;
     printf("%s: minuit number %d  bkg set %s  exposure %f  BI %f\n",snam,ipar,it->first.c_str(),it->second.exposure,it->second.value);
		 
   }

   // --------------------------------------------------------------
   // systematic errors: 3 scale factors for  systematic_name  / group 
   // multiply uncertainty by scale factor + add to resolution/peak position/efficiency
   for(itersys is = syspar.begin(); is != syspar.end(); ++is) {

     ipar = 1 + Nbkgsets + 3*is->second; // 3 system parameters for each group

     sprintf(line,"%s-res",is->first.c_str());
     gMinuit->mnparm(ipar++,line,0.,0.5,-3.,3,ierflg);
     if(ierflg) return ierflg;
     sprintf(line,"%s-pos",is->first.c_str());
     gMinuit->mnparm(ipar++,line,0,0.5, -3.,3.,ierflg);
     if(ierflg) return ierflg;
     sprintf(line,"%s-eff",is->first.c_str());
     gMinuit->mnparm(ipar++,line,0,0.5, -3.,3.,ierflg);
     if(ierflg) return ierflg;
   }
   
   if(debug<2) gMinuit->SetPrintLevel(-1); //  -1: no print out


   // now read in 90% CL threshold histogram for profile likelihood
   //
   TFile *ft= new TFile(fncuthist,"READ");
   if(!ft->IsOpen()) {
     // file does not exist ?
     printf(" could not open histogram file with 0.9 CL thresholds \n");
     TcutNbins=0;
     TcutBinWidth=0;
     return 0;
   }
   TH1F *tcut = (TH1F*) ft->Get("tcut");
   if(!tcut) {
     printf("failed to find histogram for likelihood cut, try other histo \n");
     tcut = (TH1F*) ft->Get("lim");
     if(!tcut) {printf("also failed to find other histo\n"); return -1;}
     else printf(" found other histogram\n");
   }
   TcutNbins=tcut->GetNbinsX();
   printf("reading in histo for cut on profile likelihood with %d bins\n",TcutNbins);

   TcutBinWidth=tcut->GetBinWidth(1);
   TcutValue = (Double_t*) malloc(sizeof(Double_t)*TcutNbins);
   for(i=0;i<TcutNbins;i++) TcutValue[i] = tcut->GetBinContent(i+1);
   if(debug) printf("    cut values for bins 0,1,2  %f %f %f %f bin width %f\n",TcutValue[0],
   	  TcutValue[1],TcutValue[2],TcutValue[3],TcutBinWidth);

   return 0;
}

//----------------------------------------------------------------------/
int read_data(const char *fn) {
   // read in data events from file name fn
  unsigned int  unixtime;
  Double_t  energy;
  char nam[40],dnam[40];
  char snam[30]={"read_data"};
  int i,j;

  FILE *f = fopen(fn,"r");
  if( f==0 ) {printf("%s: error opening file data event list",snam); return 2;}
   
  if(Nsets==0) {
    printf(" first read in definition of data sets \n"); return 2;}

  for(i=0;i<Nsets;i++)  for(j=0;j<DETMAX; j++) 	 part[i].nevt_bkg[j]=0;
  
  Int_t ok=0;
  // read in line by line and find the correct data set
  // according to the name (call initfit first)
  while(fscanf(f,"%s %u %lf  %s",dnam,&unixtime,&energy,nam) == 4) {
     Int_t nmatch=0;
     for(i=0;i<Nsets;i++) {
       if(unixtime < part[i].t1 || unixtime > part[i].t2) continue;
       for(j=0;j<DETMAX; j++) {
	 if(part[i].dp[j].exposure<1.e-4) continue;
	 if(strcmp(dnam,part[i].dp[j].name)==0) {  // match "detector" name

	   part[i].dp[j].Eevt[ part[i].nevt_bkg[j]++ ] = energy;
	   nmatch++;
	   if(debug>1) printf("%s: unixtime %u  energy %f data set %s number %d\n",
			      snam,unixtime,energy,nam,i);
	 }
       }
       if(nmatch!=1) 
	 printf("%s: error reading data events, set name does not match %s\n",
	      snam,dnam);
       else ok++;
       if(energy>2034 && energy<2044) printf(" evt in ROI: energy %f set %s\n",energy,nam);
     }
   }
   fclose(f);
   return 0;
}


// -------------------------------------------------------------------------
Int_t scanfit(Double_t *best_xinv, Double_t *limit90, Double_t *limit90low,
	      Double_t **xi, Double_t **fcns, Double_t **fcns_tlimit, Int_t *Nstep)
//
//   calculates test statistics for every 1/T=xinv until in steps of binwidth  until maxscan
//
// in       nothing
///         
// returns  best_xinv = best fit value invT12, 
//          limit90 = invT12 for which likelihood passes the 90% threshold 
//          limit90low = same 90% C.L. for low side 1/T in case best_xinv
//                     is larger than 0 and test statistics (1/T=0) is
//                     larger than limit at 0.
//          xi = array of xinv values 
//          fcns = array with value of test statistic for every xinv
//          fcns_limit = array of 90% limit for every xinv, from histogram
//          Nstep = number of value of xinv,fcns,fcns_limit
{

   Double_t arglist[10];
   Int_t ierflg = 0;

   Double_t  step;
   Double_t best_fcn=1.e8;
   Double_t last_fcn=1.e8;
   Double_t fcn0=-1.;
   Double_t maxscan=0.15; // max xinv in 1E-25 
   Double_t profilethreshold=CL90CUT;
   Int_t i;
   int    call_limit=3000;  // number of minuit calls
   Double_t best_amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
   Double_t  best_bkg,best_err,best_lo,best_hi;
   Int_t     best_i;
   TString   nam;
   char sn[20] = {"scanfit"};
   
   //   if(TcutBinWidth > step) step = TcutBinWidth;

   if(!TcutBinWidth) {
     // failed to read in histogram in initfit
     return -1;
   }
     
   step = binwidth;
   
   *limit90low = 0;
   // important: initialize *fcns = 0 in main program 
   // allocate memory for output arrays
   if(!(*fcns) ) {

     *Nstep = maxscan/step;
     *fcns = (Double_t *) malloc(sizeof(Double_t)* *(Nstep));
     //if(!fcns) return -99;
     *fcns_tlimit = (Double_t *) malloc(sizeof(Double_t)* *(Nstep));
     //if(!fcns) return -99;
     *xi = (Double_t *) malloc(sizeof(Double_t)* *(Nstep));
     //if(!xi) return -99;
   }


   // scan in 1/T, fix first parameter of fit = 1/T
   arglist[0]=1;
   gMinuit->mnexcm("FIX", arglist ,1,ierflg);
   if(ierflg) {printf("scanfit: command FIX failed %d\n",ierflg); return -1.*ierflg;}
   
// Now ready for minimization step

   
   for(i=0; i < *Nstep; i++){
     Double_t xinv = i*step+step/10.;

//   set 1/T used in MINUIT to value xinv
     arglist[0]=1;
     arglist[1]=xinv;
     gMinuit->mnexcm("SET PAR", arglist ,2,ierflg);
     if(ierflg) {printf("%s: command SET PAR failed %d \n",sn,ierflg); return -1.*ierflg;}

     
     // first: fix all 3 systematic errors, 1=1/T, 2=bkg index
     // 3=sigma of peak
     int irc = fixsystematic();
     if(irc!=0) return irc;

     // now pre-fit without systematics, only if any of the system. errors are studied
     arglist[0] = call_limit;
     arglist[1] = 1.e-3;   // convergence threshold
     ierflg=0;
     if(sys_sigma || sys_pos || sys_eff) gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

     if(ierflg) {  
       // fit was not successful, try again
       // repeat fit with lower precision for convergence
       arglist[0] = call_limit;
       arglist[1] = 1.e-2;
       gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
			
       if (ierflg)
	 {
	   // switch to SIMPLEX method
	   arglist[0] = call_limit;
	   arglist[1] = 1.e-3;
	   gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
	   if (ierflg)
	     // fit did not converge, stop and return with negative RC
	     printf("%s: pre-fit return code %d \n",sn,ierflg); return -1.*ierflg; 
	 }
     } // end pre-fit
     
     irc = relsystematic();
     if(irc!=0) return irc;

     // now the fit with the systematic error parameters floating
     arglist[0] = call_limit;
     arglist[1] = 1.e-6;   // precision of fit
     gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
     if(ierflg !=0 ) {
       // fit failed, try another fit
       // switch to SIMPLEX with less stringent precision
       printf("%s: MIGRAD return code %d, ",sn,ierflg);
       arglist[0] = call_limit;
       arglist[1] = 1.e-3;
       gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
       if(ierflg ==0 ) printf("  SIMPLEX succeeded with 10^-3 precision \n");
       if(ierflg !=0 ) { 
	 // switch to even poorer precision
	 arglist[0] = call_limit;
	 arglist[1] = 1.e-2;
	 gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
	 if(ierflg ==0 ) printf("  SIMPLEX succeeded with 10^-2 prec.\n");
	 if (ierflg !=0) /*printf("  SIMPLEX succeeded with less precision \n");*/
	   // fit failed, return with error code
	   printf("%s: last SIMPLEX return code %d \n",sn,ierflg); return -1.*ierflg; 
       }
     }

     // threshold for 90% limit (low/high side)
     // interpolate threshold between two bins from histogram
     Int_t nxx = xinv/TcutBinWidth;
     if(nxx>=TcutNbins-1)   profilethreshold = *(TcutValue+TcutNbins-1);
     else {
       Double_t w = (xinv-nxx*TcutBinWidth)/TcutBinWidth;
       profilethreshold = *(TcutValue+nxx) *(1-w) + *(TcutValue+nxx+1) *w; 
     }

     
     // Print results from fit
     gMinuit->mnstat(best_amin,edm,errdef,nvpar,nparx,icstat);

     if(xinv<step) fcn0 = best_amin;

     if(debug) {
       if(xinv<step) printf("%s: xinv %f: like=%f cut_t %f",sn,
		      xinv,best_amin,profilethreshold);
       else printf("%s: xinv %f: like=%f cut_t %f",sn,
		      xinv,best_amin-fcn0,profilethreshold);
     }

     *(*fcns+i)=best_amin;
     *(*fcns_tlimit+i)=profilethreshold;
     *(*xi+i)=xinv;
     
// 		 printf("fit value %f\n",best_amin);

     // search for minimum of FCN as function of 1/T=xinv
     if(best_amin < best_fcn) {
       best_fcn = best_amin;
       *best_xinv = xinv;
     }

     if(debug) {
       typedef std::map<std::string,struct bkginput>::iterator iter;
       printf(" best fit number bkg evt: ");
       for(iter it = bkgnames.begin(); it != bkgnames.end(); ++it) {
	 int p=it->second.num+1; 
	 gMinuit->mnpout(p,nam,best_bkg,best_err,best_lo,best_hi,best_i);

	 printf("bset %d b=%6.2e ",p-1,(float)best_bkg*220.*it->second.exposure);
       }
       printf("\n");
     }


   } // end loop over i<NStep, xinv

   if(debug) printf("  before final opt: best_xinv %f  like %f \n",
		    *best_xinv,best_fcn);

   // first find precisely the minimum by leaving xinv free

   // release  1/T and fix to value with lowest fcn
     arglist[0]=1;
     gMinuit->mnexcm("REL", arglist ,1,ierflg);
     if(ierflg) {printf("%s: command REL failed %d\n",sn,ierflg); return -1.*ierflg;}

     arglist[0]=1;
     arglist[1]=*best_xinv;
     gMinuit->mnexcm("SET PAR", arglist ,2,ierflg);
     if(ierflg) {printf("%s: command SET PAR failed %d \n",sn,ierflg); return -1.*ierflg;}

     // now the minimization
     arglist[0] = call_limit;
     arglist[1] = 1.e-5;
     gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
     if(ierflg !=0 ) {
       //switch to SIMPLEX
       printf("%s: MIGRAD5 return code %d, ",sn,ierflg);
       arglist[0] = call_limit;
       arglist[1] = 1.e-3;
       gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
       if(ierflg ==0 ) printf("  SIMPLEX5 succeeded with 10^-3 precision \n");
       if(ierflg !=0 ) { 
	 // switch to lower precision
	 arglist[0] = call_limit;
	 arglist[1] = 1.e-2;
	 gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
	 if(ierflg ==0 ) printf("  SIMPLEX5 succeeded with 10^-2 prec.\n");
	 if (ierflg !=0) /*printf("  SIMPLEX succeeded with less precision \n");*/
	   // fit failed, return with error code
	   printf("%s: last SIMPLEX5 return code %d \n",sn,ierflg); return -1.*ierflg; 
       }
     }

     // get fcn at minimum
     gMinuit->mnstat(best_amin,edm,errdef,nvpar,nparx,icstat);
     Int_t p=0;
     Double_t bestx;
     gMinuit->mnpout(p,nam,bestx,best_err,best_lo,best_hi,best_i);
     best_fcn = best_amin;
     *best_xinv = bestx;

     if(debug) {
       printf(" best xinv %f  like %f \n",bestx, best_amin);
     }
     
     // subtract the minimum fcn from the test statistics of all xinv values
     for(i=0; i< *Nstep; i++)
       *(*fcns+i)-=best_fcn;


   // now search for 90% limit
   // first check if xinv=0 is outside 90% interval, determine 90% lower limit
   
   Int_t k=0;
   
   if(*(*fcns) < *(*fcns_tlimit) ) {
     *limit90low = 0.;
   }
   else {
     for(k=1; k< *(Nstep); k++) {
       if(*(*fcns+k) < *(*fcns_tlimit+k) ) {
	 *limit90low = k*step - 
	   step*( *(*fcns_tlimit+k) - *(*fcns+k) ) /
	   (*(*fcns+k-1) - *(*fcns+k) );
	 break; 
       }
     }
   }

   // now search for 90% upper limit
   
   // now search for 90% upper limit

   for(i=*(Nstep) - 1; i>=k; i--){

     // increasing FCN after minimum, determine the 90% limit
     if( *(*fcns+i-1)  < *(*fcns_tlimit+i-1) && i>1) {
       // xinv of 90% C.L.  FCN increases (delta1),
       //  threshold also changes within step (delta2)
       double delta1 =  *(*fcns+i) - *(*fcns+i-1);
       double delta2 =   *(*fcns_tlimit+i) - *(*fcns_tlimit+i-1);
       *limit90 = (*(*fcns+i-1) -  *(*fcns_tlimit+i-1)) * step /(delta2- delta1);
       *limit90 +=  *(*xi + i - 1);//(i-1)*step;
    //   *limit90 += (i-1)*step;
       break;
     }
     *limit90=10; // no 90% limit in scaned interval
   }

/*


   for(i=k; i< *(Nstep); i++){

     // increasing FCN after minimum, determine the 90% limit
     if( *(*fcns+i)  > *(*fcns_tlimit+i) && i>1) {
       // xinv of 90% C.L.  FCN increases (delta1),
       //  threshold also changes within step (delta2)
       double delta1 =  *(*fcns+i) - *(*fcns+i-1);
       double delta2 =   *(*fcns_tlimit+i) - *(*fcns_tlimit+i-1);
       *limit90 = (*(*fcns+i-1) -  *(*fcns_tlimit+i-1)) * step /(delta2- delta1);
       *limit90 += (i-1)*step;
       break;
     }
     *limit90=10; // no 90% limit in scaned interval
   }

   if(debug) printf(" limit between x %f and %f,  fcn %f %f  thresh %f %f\n",
		    (i-1)*step,i*step,*(*fcns+i-1),*(*fcns+i),
		    *(*fcns_tlimit+i-1),*(*fcns_tlimit+i) );
  */  
   return 0;




} // end

// -------------------------------------------------------------------------
Double_t profile(Double_t xinv_true, Double_t &best_xinv)
// returns  teststatistic =
//               -2 log( L(xinv_true, theta_hat2) / L(xinv_hat, theta_hat))
//          theta_hat = best fit background
//          xinv_hat = best fit xinv
//          theta_hat2 = best fit background for xinv=xinv_true
//
//          input:  xinv_true
//          output:  best_xinv = best fit 1/T = xinv_hat
//          return value = 2-sided test statistic
{
  Double_t arglist[10];
  Int_t ierflg = 0;
  Int_t call_limit=6000;
  char sn[20]={"profile"};
  Int_t loopcnt=0;
  
   if(debug<2) gMinuit->SetPrintLevel(-1); //  -1: no print out
   else        gMinuit->SetPrintLevel(0); //  -1: no print out

   // release invT12 = parameter 1, to be sure
   arglist[0]=1;
   gMinuit->mnexcm("REL", arglist ,1,ierflg);
   if(ierflg) {printf("%s: command REL failed %d\n",sn,ierflg); return -1.*ierflg;}
   arglist[0]=1;
   arglist[1]=xinv_true;
   gMinuit->mnexcm("SET PAR", arglist ,2,ierflg);
   if(ierflg) {printf("%s: command SET PAR failed %d\n",sn,ierflg); return -1.*ierflg;}

   int irc = fixsystematic();
   if(irc!=0) return irc;

   int sysres   = 1 + Nbkgsets + 1;
   int sysshift = 1 + Nbkgsets + 2;
   int syseff   = 1 + Nbkgsets + 3;

// Now the minimization step for  bkg and invT floating, no system errors 
   arglist[0] = call_limit;
   arglist[1] = 1.e-4;
   ierflg=0;
   if(sys_sigma || sys_pos || sys_eff)
     gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
   if(ierflg !=0 ) { 
     printf("%s: pref-fit SIMPLEX return code %d \n",sn,ierflg); 
     arglist[0] = call_limit;
     arglist[1] = 1.e-3;
     gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
     if(debug && ierflg==0) printf(" but ... converged with 10^-3 precision\n"); 
     else if(ierflg!=0){
       printf("%s: pref-fit SIMPLEX return code %d \n",sn,ierflg); 
       arglist[0] = call_limit;
       arglist[1] = 1.e-1;
       gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg); 
       if(ierflg==0) printf(" but ... converged with 10^-1 precision\n"); 
       else          {return -1.*ierflg;}
     }
   }
    
   //        ------------------------------------------
   // now realease system error parameters and do the real fit

   irc = relsystematic();
   if(irc!=0) return irc;
	
   // Now ready for minimization step, first start with SIMPLEX
   arglist[0] = call_limit;
   arglist[1] = 1.e-3;
   gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
   if(ierflg !=0 ) { 
     if(debug) printf("%s: SIMPLEX RC %d ",sn,ierflg); 
     arglist[0] = call_limit;
     arglist[1] = 1.e-1;
     gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
     if(debug && ierflg==0) printf(" ... SIMPLEX converged with 10^-1 precision\n"); 
     else if(ierflg!=0)
       {
	 printf(" SIMPLEX did not even converge with 10^-1 precision\n");
	 gMinuit->SetPrintLevel(1);
	 gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
	 return -1.*ierflg;
       }
   }

 repeat:
   // repeat minimization with MIGRAD
   arglist[0] = call_limit;
   arglist[1] = 1.e-5;
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   if(debug && ierflg == 0) printf("           MIGRAD converged at 1E-4 precision all floating\n");
   if(ierflg !=0 ) { 
     if(debug)printf("%s: MIGRAD0 return code %d try again SIMPLEX \n",sn,ierflg);
     // go back to SIMPLEX fit
     arglist[0] = call_limit;
     arglist[1] = 1.e-3;
     gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
     if(ierflg !=0 ) {
       printf("%s: final fit with SIMPLEX did not converge, error %d\n",sn,ierflg);
       return -1.*ierflg;
     }
     debug=0;
   }



// get best invT12 = invT12^hat and log likelihood at min = best_amin
   Double_t  best_err,best_lo,best_hi;
   Double_t best_amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
   Int_t best_i;
   TString   nam;
   gMinuit->mnstat(best_amin,edm,errdef,nvpar,nparx,icstat);
   gMinuit->mnpout(0,nam,best_xinv,best_err,best_lo,best_hi,best_i);
   if(debug) {
     //printf("%s: out %f %f %f %d %d %d \n",sn,best_amin,edm,errdef,nvpar,nparx,icstat);
     printf("%s: best invT12= %e  best Log likelihood=%f\n",sn,best_xinv,best_amin);
   }


     //printf("Before Replacing: %f \n",best_xinv);
   
// check if invT12=0 has better likelihood in case best_xinv is close to 0
   if(best_xinv<1.e-3 && best_xinv>0) {

     arglist[0]=1;
     gMinuit->mnexcm("FIX", arglist ,1,ierflg);
     if(debug&&ierflg) {printf("%s: command FIX failed %d\n",sn,ierflg); return -1*ierflg;}

     arglist[0]=1;
     arglist[1]=0;
     gMinuit->mnexcm("SET PAR", arglist ,2,ierflg);
     if(debug&&ierflg) {printf("%s: command SET PAR failed %d \n",sn,ierflg); return -1*ierflg;}


     arglist[0] = call_limit;
     arglist[1] = 1.e-5;
     gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
     if(ierflg !=0 ) { 
       if(debug)printf("%s: fixing xinv=0, MIGRAD2 return code %d ",sn,ierflg);
       arglist[0] = call_limit;
       arglist[1] = 1.e-3;
       gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
       if(debug&&ierflg != 0) printf("     SIMPLEX also failed at 1E-3  precision \n");
       if(debug&&ierflg == 0) printf("     SIMPLEX converged at 1E-3  precision \n");
     }
     // what is FCN at xinv=0?
     Double_t current_amin=1.e8;
     if(ierflg==0) gMinuit->mnstat(current_amin,edm,errdef,nvpar,nparx,icstat);

     if(current_amin<best_amin) {
         //printf("Replacing best fit %f with 0 \n",best_xinv);
       if(debug)printf(" xinv=0 has FCN=%f  while best_xinv=%e had FCN=%f\n",
	      current_amin,best_xinv,best_amin);
       best_amin=current_amin;
       best_xinv=0;
     }

   } // best_invT12 < 1E-3 --> check for 1/T=0


   if(debug)printf("%s  last FCN at opt %f %f %f %f xinv %f  mu %f distQbb %f\n",sn,last_fcn[0],
		   last_fcn[1],last_fcn[2],last_fcn[3],best_xinv,last_fcn[4],distQbb);

   
   // calculate log likelihood at true_xinv, fix 1/T to true value

   arglist[0]=1;
   gMinuit->mnexcm("FIX", arglist ,1,ierflg);
   if(ierflg) {printf("%s: command FIX failed %d\n",sn,ierflg); return -1.*ierflg;}

   arglist[0]=1;
   arglist[1]=xinv_true;
   gMinuit->mnexcm("SET PAR", arglist ,2,ierflg);
   if(ierflg) {printf("%s: command SET PAR failed %d \n",sn,ierflg); return -1.*ierflg;}

   arglist[0] = call_limit;
   arglist[1] = 1.e-5;
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

   if(ierflg !=0 ) {
     // fit failed --> repeat fit 
     
     if(debug)printf("%s: MIGRAD3 RC %d ",sn,ierflg); 
     // try SIMPLEX at reduced precision
     arglist[0] = call_limit;
     arglist[1] = 1.e-3;

     gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
     if(debug && ierflg==0) printf(" ... but SIMPLEX converged \n");
     else if(ierflg!=0) { 
       printf("        also SIMPLEX  failed\n"); 
       return -1.*ierflg;
     }
     Double_t current_amin;
     gMinuit->mnstat(current_amin,edm,errdef,nvpar,nparx,icstat);

     // check if FCN is better than from MIGRAD call at "best xinv"
     // or if xinv are the same but FCN differ from SIMPLEX and MIGRAD
     if(current_amin < best_amin - 1.e-5 ||
	(fabs(xinv_true- best_xinv)<1.e-5 && fabs(current_amin-best_amin)>1.e-5))
       {
       arglist[0] = 10*call_limit;
       arglist[1] = 1.e-3;
       gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

       if(debug&&ierflg==0) printf(" ... and now MIGRAD as well \n");
       else if(ierflg!=0) {
	 printf("     2nd attempt with MIGRAD failed again \n");
	 return -1.*ierflg;
       }
     }// end current_amin determination
   }


   if(debug)printf("%s  last FCN contributions true xinv %f %f %f %f xinv %f mu %f\n",sn,
	  last_fcn[0], last_fcn[1],last_fcn[2],last_fcn[3],xinv_true,last_fcn[4]);

   // FCN at xinv_true fixed
   Double_t current_amin;
   gMinuit->mnstat(current_amin,edm,errdef,nvpar,nparx,icstat);
   if(debug) printf("%s: best log like for true %f is %f t=%f\n",
			sn,xinv_true,current_amin,current_amin-best_amin);

   // now check if fit with all parameters floating really found minimum
   // if not, try again to find the minimum
   
   if(current_amin < best_amin - 1.e-5) {
     loopcnt++;
     if(debug)printf("%s: negative t: %f - %f, xinv_best %e xinv_true %e\n",
	    sn,current_amin,best_amin, best_xinv,xinv_true);
     arglist[0]=1;
     gMinuit->mnexcm("REL", arglist ,1,ierflg); // releast 1/T fit parameter
     if(ierflg) {printf("%s: command REL4 failed %d\n",sn,ierflg); return -1.*ierflg;   }
     if(loopcnt<3) goto repeat;
   }


   // extract fit parameters of systematic errors
   
   Double_t xval;
   gMinuit->mnpout(sysres-1,nam,xval,best_err,best_lo,best_hi,best_i);
   last_fcn[10] = xval;  // resolution
   gMinuit->mnpout(sysshift-1,nam,xval,best_err,best_lo,best_hi,best_i);
   last_fcn[11] = xval;  // energy shift
   gMinuit->mnpout(syseff-1,nam,xval,best_err,best_lo,best_hi,best_i);
   last_fcn[12] = xval;  // efficiency


   return current_amin-best_amin;


} // end

// -------------------------------------------------------------------------
Double_t profile1(Double_t xinv_true, Double_t &best_xinv)
// returns  teststatistic =
//               -2 log( L(xinv_true, theta_hat2) / L(xinv_hat, theta_hat))
//                    for xinv_true > xinv_hat, = 0  otherwise
//          theta_hat = best fit background
//          xinv_hat = best fit xinv
//          theta_hat2 = best fit background for xinv=xinv_true

//          return 1-sided test statistic:
//                 

{


  Double_t t = profile(xinv_true, best_xinv);

   if(best_xinv > xinv_true) return 0;

   else   return t;


} // end


 
//---------------------------------------------------------
int fixsystematic(){

  // fix systematic error MINUIT parameter
  Double_t arglist[10];
  Int_t ierflg = 0;
  char sn[20]={"fixsystemaric"};

  
  int sysres   = 1 + Nbkgsets + 1;
  int sysshift = 1 + Nbkgsets + 2;
  int syseff   = 1 + Nbkgsets + 3;

  for(int i=0; i<ksys; i++) {
  
    arglist[0]=sysres + i*3;
    gMinuit->mnexcm("FIX", arglist ,1,ierflg);
    if(ierflg) {printf("%s: command FIX failed %d\n",sn,ierflg); return -1.*ierflg;}

    arglist[0]=sysres + i*3;
    arglist[1]=0.;
    gMinuit->mnexcm("SET PAR", arglist ,2,ierflg);
    if(ierflg) {printf("%s: command SET PAR failed %d\n",sn,ierflg); return -1.*ierflg;}

    // central value of peak
    arglist[0]=sysshift + i*3;
    gMinuit->mnexcm("FIX", arglist ,1,ierflg);
    if(ierflg) {printf("%s: command FIX failed %d\n",sn,ierflg); return -1.*ierflg;}
    arglist[0]=sysshift + i*3;
    arglist[1]=0;
    gMinuit->mnexcm("SET PAR", arglist ,2,ierflg);
    if(ierflg) {printf("%s: command SET PAR failed %d\n",sn,ierflg); return -1.*ierflg;}

    // efficiency
    arglist[0]=syseff + i*3;
    gMinuit->mnexcm("FIX", arglist ,1,ierflg);
    if(ierflg) {printf("%s: command FIX failed %d\n",sn,ierflg); return -1.*ierflg;}
    arglist[0]=syseff + i*3;
    arglist[1]=0;
    gMinuit->mnexcm("SET PAR", arglist ,2,ierflg);
    if(ierflg) {printf("%s: command SET PAR failed %d\n",sn,ierflg); return -1.*ierflg;}
			

  }
  return 0;
}
//------------------------------------------------------------
int relsystematic(){

  // release systematic error MINUIT parameter for fit
  Double_t arglist[10];
  Int_t ierflg = 0;
  char sn[20]={"relsystemaric"};
  
  int sysres   = 1 + Nbkgsets + 1;
  int sysshift = 1 + Nbkgsets + 2;
  int syseff   = 1 + Nbkgsets + 3;
// depending on the bit pattern release system. error parameters

  for(int i=0; i<ksys; i++) {
  
    arglist[0]=sysres + i*3;
    if(sys_sigma) gMinuit->mnexcm("REL", arglist ,1,ierflg);
    if(ierflg) {printf("%s: command REL failed %d\n",sn,ierflg); return -1.*ierflg;}
       
    arglist[0]=sysshift + i*3;
    if(sys_pos) gMinuit->mnexcm("REL", arglist ,1,ierflg);
    if(ierflg) {printf("%s: command REL failed %d\n",sn,ierflg); return -1.*ierflg;}
    
    arglist[0]=syseff + i*3;
    if(sys_eff) gMinuit->mnexcm("REL", arglist ,1,ierflg);
    if(ierflg) {printf("%s: command REL failed %d\n",sn,ierflg); return -1.*ierflg;}

  }
  return 0;
}
