// profile likelihood fit of phase I data 
// 
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
#include "TGraph.h"
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


//______________________________________________________________________________
main(int argc, char **argv)
{
// The z values
  Int_t i,k,j;
  Double_t best_invT12,tim,limit90,limit90low;
  TRandom3 *xr = new TRandom3();

  FILE *floglike;
  FILE *ft12;
  char line[200];

  char fn[40];

  double maxinv=0.15;
  double binwidth=get_binwidth();
  int nbins=(maxinv+1e-7)/binwidth;
  
  TFile *fout = new TFile("../../ana/legend_gerda_mjd/t12limit-tmp.root","RECREATE");

  TH1F  *hlim  = new TH1F("dlim","dlim",nbins+1,-binwidth/2.,(nbins+0.5)*binwidth);

  Double_t invT,t0;
  TTree outTree("toyMC","toyMC");
        outTree.Branch("invT12",&invT);
        outTree.Branch("t",&t0);
  


  i = initfit(2,"../../data_sets/legend_gerda_mjd/parameter.txt","../../data_sets/legend_gerda_mjd/bkg.txt",7,"../../ana/legend_gerda_mjd/pval-ph2-tmp.root");  // init MINUIT

  if(i) exit(1);

  sprintf(fn,"../../data_sets/legend_gerda_mjd/events.txt");
  // sprintf(fn,"events_ph2.txt");
  i = read_data(fn);  // read in event list
  if(i) exit(1);

  Double_t *xinv=0;
  Double_t *t=0;
  Double_t *limit=0;
  Int_t  Nstep;


  if(1) {
    // determine likelihood as fct of  1/T=xinv,  "t" = array of likelihood values
    i = scanfit(&best_invT12, &limit90, &limit90low, &xinv, &t, &limit, &Nstep);
    printf("\n best_invT12 %f  low limit %f  high limit %f",
	   best_invT12,limit90low,limit90);
    if(limit90>0) printf(" T12lower %f \n",1/limit90);
  }
  //double t1=profile(0.0546,best_invT12);
  //return 0;
  

  printf("   t at %f is %f \n",xinv[62],t[62]); 

  // store entire profile dist in tree
  for(Int_t j=0; j < Nstep; j++) {
    t0 = t[j];
    invT = xinv[j];
    hlim->Fill(invT,t0);
    if(t0>=0) outTree.Fill();

       	//if(t0>0) printf("xinv %f t=%f \n",xinv[j],t[j]);
  }

  fout->cd();
      
  TGraph *tf = new TGraph(Nstep,xinv,t);
  tf->SetTitle("data_t");
  tf->Write();
  fout->Write();
  fout->Close();


      
} // end main()
//-------------------------------------------------------------------------------
