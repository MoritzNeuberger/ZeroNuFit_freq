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
#include "TH1D.h"
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
  xr->SetSeed(0);  // argument =0: different seed with every start, !=0: same seed every start
  Int_t nbin=1000;

  char *outname=NULL;
  char *dirname=NULL;
  int parser=0;

  
  for(parser=1; parser<argc; parser++)
    {
      if(strcmp(argv[parser],"-o")==0) outname=argv[++parser];
      if(strcmp(argv[parser],"-d")==0) dirname=argv[++parser];
      //1else break;
    }

  
  if(outname==NULL){ 
    printf("option -o missing, using default\n"); 
    outname = (char *) malloc(20);
    sprintf(outname,"%s","pvalue.root");
  }
 
  if(dirname==NULL){ 
    printf("option -d missing, using default\n"); 
    dirname = (char *) malloc(20);
    sprintf(dirname,"%s","./");
  }
 

  
  //TFile *f = new TFile("../../ana/legend_gerda_mjd/phase2abc_sys.root");
  //TFile *f = new TFile("../../ana/legend_gerda_mjd/pval-ph2-tmp.root");
  //TH2F* th = (TH2F*)f->Get("t"); 
  //  if(!th) {printf("test statistic hist not in ROOT file\n"); return 1;}


  char fparameter[200]; // parameter for partitions
  sprintf(fparameter,"%s/parameter.txt",dirname);
  char bparameter[200];  // background index
  sprintf(bparameter,"%s/bkg.txt",dirname);
  char rparameter[200];  // input root file for 90% cut histogram
  sprintf(rparameter,"./out/output_cover2.root");
//  rparameter="/mnt/atlas01/users/neuberger/L200_stat_analysis/bernhards_stat_tool/ana/legend_gerda_mjd/output_cover2.root"  

  i = initfit(0,fparameter,bparameter,7,rparameter);  // arguments: debug, file_name, bit_pattern_systematic
  //  i = initfit(0,"parameter_ph2.txt","bkg_ph2.txt",0,"cover2-ph2-0.root");  // arguments: debug, file_name, bit_pattern_systematic
  if(i) exit(1);

  TFile *fout = new TFile(Form("%s",outname),"RECREATE");
  TH1F  *hh  = new TH1F("invt12","invt12",nbin,-0.0005,1.9995);
  TH2F  *ht  = new TH2F("t12","t12",12,-0.5,11.5,150,0.,15.);
  TH1F  *hb  = new TH1F("t12best","t12best",150,0.,3.);

  Double_t invT,profile;
  TTree outTree("toyMC","toyMC");
        outTree.Branch("invT12",&invT);
        outTree.Branch("t",&profile);
  
  TH1D* hist=0;

  int  n1sig, n2sig;
  double distQbb;
  
   
  // generate toy MC without signal and calculate limit
  // the median of the limits is the sensitivity for setting a limit
  for(k=0;k<100;k++) {

      Double_t *xinv=0;
      Double_t *t=0;
      Double_t *limit=0;
      Int_t  Nstep;
      Double_t tim=1;
      
      int n=generate_sample(xr,tim,0); // true invT12=0, 2nd argument of fct ignored
      if(n<0) exit(10);
      int ns = nsignal();


      // calculate likelihood fct(1/T) and 90% limit, 2-sided
      i = scanfit(&best_invT12, &limit90, &limit90low,
		  &xinv,&t,&limit,&Nstep);
      if(i!=0) continue;
      n1sig = get_n1sig();
      n2sig = get_n2sig();
      distQbb = get_distQbb();
      printf("%d: best fit invT12 %e limit %e lowside_limit %e ntot %d nsig %d n1sig %d n2sig %d distQbb %4.2f\n",
	     k,best_invT12,limit90,limit90low,n,ns,n1sig,n2sig,distQbb);
	if(n1sig==0) printf("   t at %f is %f, at %f is %f \n",
				xinv[28],t[28],xinv[50],t[50]); 
	
      hb->Fill(best_invT12);
      hh->Fill(limit90);     // fill histogram with 90% limit for median sensitivity
      if(limit90>0) ht->Fill((Double_t)j,1./limit90);

      // store entire likelihood  dist in tree
      for(Int_t j=0; j < Nstep; j++) {
	profile = t[j];
	invT = xinv[j];
	if(profile>=0) outTree.Fill();

	//	if(profile>0) printf("xinv %f t=%f pvalue=%f\n",xinv[j],t[j],profile);
      }
      //    tt->Fill();

  } // end loop over toy MC samples
   
  // now find median of limits, k=number of toy MCs
  Double_t median=0;  
    for(i=0;i<nbin;i++) {
      median += hh->GetBinContent(i+1);
      if(median>(k+1.)/2.) break;
    }
    const TAxis *xax = hh->GetXaxis();
    median = xax->GetBinCenter(i);
    printf("time %f median T12 %f \n",tim,1/median);
  

  gROOT->ls();
  //  outTree.Write();
  fout->Write();
//
  fout->Close();



} // end main()
//-------------------------------------------------------------------------------
