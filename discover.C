// profile likelihood fit of phase I data 
// 0.27% for 3 sigma
#include "TROOT.h"
#include "TCanvas.h"
#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <time.h>
#include "TMinuit.h"
#include <math.h>
#include <stdlib.h>
#include "TMath.h"
#include "TRandom3.h"
#include "pf.h"

int fcmp(const void *in1, const void *in2);


//______________________________________________________________________________
 main(int argc, char **argv)
{
// The z values
  Int_t i,k,j;
  Double_t best_invT12,tim,limit90,limit90low;
  TRandom3 *xr = new TRandom3();

  Int_t maxgen2=20000;
  Double_t  binwidth=get_binwidth(); 
  Double_t  maxinv = 0.2000;
  
  xr->SetSeed(0);

  char *outname=NULL;
  char *dirname=NULL;
  int parser=0;

  for(parser=1; parser<argc; parser++)
    {
      if(strcmp(argv[parser],"-o")==0) outname=argv[++parser];
      if(strcmp(argv[parser],"-d")==0) dirname=argv[++parser];
      else break;
    }
  
  if(outname==NULL){ 
    printf("option -o missing, using default\n"); 
    outname = (char *) malloc(20);
    sprintf(outname,"%s","discovertmp.root");
  }
  if(dirname==NULL){ 
    printf("option -d missing, using default\n"); 
    dirname = (char *) malloc(20);
    sprintf(dirname,"%s","./");
  }
 

  char fparameter[200]; // parameter for partitions
  sprintf(fparameter,"%s/parameter_all.txt",dirname);
  char bparameter[200];  // background index
  sprintf(bparameter,"%s/bkg_all.txt",dirname);
  char rparameter[200];  // input root file for 90% cut histogram
  sprintf(rparameter,"%s/pval-0.root",dirname);
  
  // initialize MINIUT
  i = initfit(0,fparameter,bparameter,7,rparameter);  // arguments: debug, file_name, bit_pattern_systematic
  
  if(i) exit(1);

  TFile *fout = new TFile(Form("%s",outname),"RECREATE");

  Double_t xinv,t0;
  TTree outTree("discover_toyMC","discover_toyMC");
        outTree.Branch("invT12",&xinv);
        outTree.Branch("t",&t0);


  //----------------------------------------------------
    
  Int_t nbins = (maxinv+0.000001)/binwidth;

  TH2F  *ht0  = new TH2F("t0","t0",nbins+1,-binwidth/2.,maxinv+0.5*binwidth,20000,0.,40.);

  
  // throw  some toy MC for every value of 1/T  to calculate
  // the test statistics at 1/T=0, needed to estimate the discovery sensitivity
  // store result in ntuple,  analysis done in pval.C
  for(i=0;i<nbins;i++){
    xinv = i*binwidth+binwidth/5.;
    for(k=0;k<maxgen2;k++) {
      generate_sample(xr,0,xinv);
      Double_t best_xinv2;
      t0  = profile(0, best_xinv2); // what test statistics for 1/T=0
      printf("for xinv=0   best_xinv %f   t %e\n",best_xinv2,t0);
      if(t0<-1.e-5){printf(" fit failed \n"); continue ;}
      if(t0<0)  t0=0;
      ht0->Fill(xinv,t0);
      outTree.Fill();
    }
    printf("finish %f\n",xinv);
  }    
    printf("finish \n");

  fout->Write();
//
  fout->Close();



} // end main()
//-------------------------------------------------------------------------------
int fcmp(const void *in1, const void *in2) {
  double *d1 = (double *) in1;
  double *d2 = (double *) in2;
  if(*(d1) < * (d2) ) return -1;
  else if (*(d1) > * (d2) ) return 1;
  return 0;
}
