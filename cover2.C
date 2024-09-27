// profile likelihood fit of phase I data 
// 0.27% for 3 sigma
#include "TROOT.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLine.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
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
// calculates the profile likelihoods, fills an ntuple with values
// scans over a certain range in 1/T

// input options
//    -o  file_name for output ROOT file
//    -s  dir_name  for directory with input ASCI file names
  
  
  Int_t i,k,j;
  Double_t best_invT12,tim,limit90,limit90low;
  TRandom3 *xr = new TRandom3();
  Int_t     maxgen=100;  // number of toy MC data samples for every 1/T
  Double_t  binwidth;      //  incremental step in 1/T
  Double_t  maxinv = 0.151; // max value of scan in 1/T
  Double_t tval[maxgen];    // array of profile likelihoods 

  char *outname=NULL;
  char *dirname=NULL;
  int parser=0;

  xr->SetSeed(0);
    printf("%i",xr->GetSeed());

  binwidth = get_binwidth();
  
  for(parser=1; parser<argc; parser++)
    {
      if(strcmp(argv[parser],"-o")==0) outname=argv[++parser];
      if(strcmp(argv[parser],"-d")==0) dirname=argv[++parser];
//      else break;
    }
  
  if(outname==NULL){ 
    printf("option -o missing, using default\n"); 
    outname = (char *) malloc(20);
    sprintf(outname,"%s","cover2tmp.root");
  }
  if(dirname==NULL){ 
    printf("option -d missing, using default\n"); 
    dirname = (char *) malloc(20);
    sprintf(dirname,"%s","./");
  }
 

  
  Int_t nbins = (maxinv+0.000001)/binwidth;
  TFile *fout = new TFile(Form("%s",outname),"RECREATE");
  TH2F  *hh  = new TH2F("invt12_fit","invt12_fit",nbins+1,-binwidth/2.,
			maxinv+binwidth/2.,2001,-0.001,2.001);
  TH2F  *ht  = new TH2F("t","t",nbins+1,-binwidth/2.,maxinv+0.5*binwidth,2000,0.,20.);
  TH2F  *hn  = new TH2F("n","n",nbins+1,-binwidth/2.,maxinv+0.5*binwidth,150,-0.5,149.5);
  TH1F  *hlim  = new TH1F("lim","lim",nbins+1,-binwidth/2.,maxinv+0.5*binwidth);
  TH1F *hsys1 = new TH1F("resolution","resolution",100,-2.,2.);
  
  Double_t invT,profilet;
  Double_t best_xinv;
  Double_t xsig;
  int      fail_count=0;

  int n1sig,n2sig;
  TTree outTree("toyMC","toyMC");
  outTree.Branch("invT12",&invT);
  outTree.Branch("t",&profilet);
  outTree.Branch("bestfit",&best_xinv);
  outTree.Branch("nsignal",&xsig);
  outTree.Branch("n1sig",&n1sig);
  outTree.Branch("n2sig",&n2sig);
	
  char fparameter[200]; // parameter for partitions
  sprintf(fparameter,"%s/parameter.txt",dirname);
  char bparameter[200];  // background index
  sprintf(bparameter,"%s/bkg.txt",dirname);
  char rparameter[200];  // input root file for 90% cut histogram
  sprintf(rparameter,"%s/pval-0.root",dirname);
  
  i = initfit(0,fparameter,bparameter,7,rparameter);  // arguments: debug, file_name, bit_pattern_systematic
  //  i = initfit(0,"parameter_ph2.txt","bkg_ph2.txt",0,"cover2_ph2.root");  // arguments: debug, file_name, bit_pattern_systematic
  if(i) exit(1);

  double *fcnd = fcndebug();
  
  Double_t xinv=0.;

  for(xinv=binwidth/10.; xinv<=maxinv;  ) {
    
       Int_t cover=0;
    Double_t avg_inv=0;

    for(k=0;k<maxgen;k++) {

      int n;

      n=generate_sample(xr,0,xinv);
      if(n<0) return 2;
      n1sig = get_n1sig();
      n2sig = get_n2sig();

      Double_t t = profile(xinv, best_xinv);
    //printf("Teststatistic: %f \n",t);
      tval[k]=t;
      if(t<-0.001) {
	printf("profile returned with %f \n",t);
	fail_count++;
	if(fail_count>100) exit(1);
	else   continue;
      }
      if(t<0)  t=0.000001;
      int nn = nsignal();

      hh->Fill(xinv,best_xinv);
      ht->Fill(xinv,t);
      hn->Fill(xinv,(float) n);
      avg_inv+=best_xinv;

      profilet=t;
      invT=xinv;
      if(best_xinv<xinv) invT*=-1.;
      xsig = (double) nsignal();

      //printf("true_xinv %f  best_xinv %f   t %f n1sig %d n2sig %d true_signal %d\n",	     xinv,best_xinv,t,n1sig,n2sig,nn);
      
      hsys1->Fill(*(fcnd+10));
      outTree.Fill();  // fill the pair  xinv, profile_likelihood t   in the  output tree

    }

    
    qsort((void *)tval, maxgen, sizeof(double), fcmp);
    Int_t ib = 1+9*maxgen/10;
    Int_t ib2 = 1+997*maxgen/1000;
    Double_t l = tval[ib]*1.005; // factor 1.005 added 
    Double_t l2 = tval[ib2]*1.005; // factor 1.005 added 
    j = xinv/binwidth;
    hlim->SetBinContent(j+1,l);
    printf(" 0.9 limit for xinv = %f at %f, avg fitted xinv %f t(3sigma) %e \n",
	   xinv,l,avg_inv/maxgen,l2);

    
    xinv+=binwidth;
    
    //if(xinv<0.015) xinv += binwidth; // small offset from 0
    //else if(xinv<0.15) xinv += 2*binwidth; // small offset from 0
    //else xinv += 5.*binwidth;

  } // end loop xinv


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
