
void LoadGERDAStyle () 
{
  cout << "*** loading GERDA rootlogon" << endl;

  // define and load gerda plot style
  TStyle *gerdaStyle  = new TStyle("gerdaStyle"," gerda specific root style");

  // use plain black on white colors
  gerdaStyle->SetFrameBorderMode(0);
  gerdaStyle->SetCanvasBorderMode(0);
  gerdaStyle->SetPadBorderMode(0);
  gerdaStyle->SetPadColor(0);
  gerdaStyle->SetCanvasColor(0);
  gerdaStyle->SetStatColor(0);
  //gerdaStyle->SetFillColor(0);        
  gerdaStyle->SetPalette(kGreyScale); 

  // set the paper & margin sizes
  gerdaStyle->SetPaperSize(20,26);
  gerdaStyle->SetPadTopMargin(0.0125);
  gerdaStyle->SetPadRightMargin(0.0125);
  gerdaStyle->SetPadBottomMargin(0.1);
  gerdaStyle->SetPadLeftMargin(0.1);                

  // default canvas positioning                                                 
  gerdaStyle->SetCanvasDefX(20);
  gerdaStyle->SetCanvasDefY(20);
  gerdaStyle->SetCanvasDefH(476);
  gerdaStyle->SetCanvasDefW(698);

  //fonts (42 -> sans serif, 132 -> serif)
  //43, 133  for absolute size in pts...
  gerdaStyle->SetTextFont(42);
  gerdaStyle->SetTextSize(0.08);
  gerdaStyle->SetLabelFont(42,"x");
  gerdaStyle->SetLabelFont(42,"y");
  gerdaStyle->SetLabelFont(42,"z"); 
  gerdaStyle->SetLabelSize(0.05,"x");
  gerdaStyle->SetTitleSize(0.05,"x");
  gerdaStyle->SetLabelSize(0.05,"y");
  gerdaStyle->SetTitleSize(0.05,"y");
  gerdaStyle->SetLabelSize(0.05,"z");
  gerdaStyle->SetTitleSize(0.05,"z");

  // do not display any of the standard histogram decorations
  gerdaStyle->SetOptStat(false);
  gerdaStyle->SetOptTitle(false);
  gerdaStyle->SetOptFit(0);
  gerdaStyle->SetLineScalePS(2);

  // do other things
  //...

  printf("- default style set to \"gerdaStyle\"\n");
  gROOT->SetStyle("gerdaStyle"); //uncomment to set this style                                 
}

void drawpValue () {
  
  LoadGERDAStyle ();

  // read in TGraph  generate by pval.C and plots them
  
  //  TFile *f = new TFile("pval_sys.root");
  TFile *f = new TFile("pval-0-19.root");
  TGraph *pvalueMedian = (TGraph*)f->Get("g50");
  pvalueMedian->SetLineWidth(2);
  pvalueMedian->SetLineStyle(9);

  
  TGraph *pvalLow68_in = (TGraph*)f->Get("g16");
  TGraph *pvalHigh68_in = (TGraph*)f->Get("g84");
  vector<double> x68;
  vector<double> y68;
  vector<double> ErrLow68;
  vector<double> ErrHigh68;
  double x, y, median;
  for (int i=0; i<pvalLow68_in->GetN(); i++) {
    pvalueMedian->GetPoint(i,x,median);
    y68.push_back(median);
    x68.push_back(x);
    pvalLow68_in->GetPoint(i,x,y);
    ErrLow68.push_back(median-y);
    pvalHigh68_in->GetPoint(i,x,y);
    ErrHigh68.push_back(y-median);
  }
  TGraphAsymmErrors  *pval68 = new TGraphAsymmErrors(x68.size(),&(x68[0]),&(y68[0]),0,0,&(ErrLow68[0]),&(ErrHigh68[0]));
  pval68->SetFillColor(kGreen);
  pval68->SetLineColor(kBlack);
  pval68->SetLineWidth(2);
  pval68->SetLineStyle(9);
//   pval68->RemovePoint(0);
  
  TGraph *pvalLow90_in = (TGraph*)f->Get("g05");
  TGraph *pvalHigh90_in = (TGraph*)f->Get("g95");
  vector<double> x90;
  vector<double> y90;
  vector<double> ErrLow90;
  vector<double> ErrHigh90;
  for (int i=0; i<pvalLow90_in->GetN(); i++) {
    pvalueMedian->GetPoint(i,x,median);
    y90.push_back(median);
    x90.push_back(x);
    pvalLow90_in->GetPoint(i,x,y);
    ErrLow90.push_back(median-y);
    pvalHigh90_in->GetPoint(i,x,y);
    ErrHigh90.push_back(y-median);
  }
  TGraphAsymmErrors  *pval90 = new TGraphAsymmErrors(x90.size(),&(x90[0]),&(y90[0]),0,0,&(ErrLow90[0]),&(ErrHigh90[0]));
  pval90->SetFillColor(kYellow);
  pval90->SetLineColor(kBlack);
  pval90->SetLineWidth(2);
  pval90->SetLineStyle(9);
//   pval90->RemovePoint(0);
  
//  TFile *fData = new TFile ("t12limit-all.root");
  TGraph* pvalData = (TGraph*)f->Get("data");
  pvalData->SetLineWidth(2);
  pvalData->SetLineStyle(1);

  
  //TFile *fd= new TFile("discover_phase2abc.root");
  TFile *fd= new TFile("discover-0-9.root");
  TGraph *discover = (TGraph*) fd->Get("discover_prob");
  discover->SetLineWidth(2);
  discover->SetLineStyle(4);
  
  TGraph *discover2 = (TGraph*) fd->Get("discover_pval");
  discover2->SetLineWidth(2);
  discover2->SetLineStyle(2);
  discover2->SetLineColor(kBlue);
  // set points to 0.5 instead of 1
  double binwidth=0.00075;
  //for(int ii=0; ii<5; ii++)
    //    discover2->SetPoint(ii,(ii+0.1)*binwidth,0.5); // smoothen graph here!!!

  
  //  pvalueMedian->RemovePoint(0);
  //   pvalData->RemovePoint(0);
  
  TCanvas *c1 = new TCanvas("c1", "c1",330,233,1012,500);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c1->Range(-0.1132232,-1.102507,1.014009,4.015495);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetLogy();
  c1->SetRightMargin(0.01287129);
  c1->SetTopMargin(0.01118568);
  c1->SetBottomMargin(0.1565996);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(pval90,"3 C");
  mg->SetTitle(" ;1/T [10^{-25} 1/yr];p-value");
  mg->Add(pval68,"3 C");
  mg->Add(pvalueMedian,"C");
  mg->Add(pvalData,"C");
  mg->Add(discover,"C");
  mg->Add(discover2,"C");
  
  mg->Draw("a");
  mg->GetYaxis()->SetRangeUser(0.0009,4.2);
  mg->GetXaxis()->SetRangeUser(0.001,0.122);
  mg->Draw("aC");
  
  TLine *limit = new TLine(0.001,0.1,0.0548,0.1);
  limit->SetLineWidth(2);
  limit->SetLineColor(kRed);
  TLine *limit2 = new TLine(0.0548,0.1,0.0548,0.0009);
  limit2->SetLineWidth(2);
  limit2->SetLineColor(kRed);
  
  TLine *sensit = new TLine(0.001,0.1,0.0555,0.1);
  sensit->SetLineWidth(2);
  sensit->SetLineColor(kRed);
  sensit->SetLineStyle(2);
  TLine *sensit2 = new TLine(0.0555,0.1,0.0555,0.0009);
  sensit2->SetLineWidth(2);
  sensit2->SetLineColor(kRed);
  sensit2->SetLineStyle(2);
  
  TLine *disco = new TLine(0.001,0.5,0.0805,0.5);
  disco->SetLineWidth(1);
  disco->SetLineColor(1);
  disco->SetLineStyle(4);
  TLine *disco2 = new TLine(0.0805,0.5,0.0805,0.0009);
  disco2->SetLineWidth(1);
  disco2->SetLineColor(1);
  disco2->SetLineStyle(4);
  TLine *disco3 = new TLine(0.001,0.003,0.0805,0.003);
  disco3->SetLineWidth(1);
  disco3->SetLineColor(1);
  disco3->SetLineStyle(4);
  
  limit->Draw("same");
  sensit->Draw("same");
  limit2->Draw("same");
  sensit2->Draw("same");
  disco->Draw("same");
  disco2->Draw("same");
  disco3->Draw("same");
  
  gPad->RedrawAxis();
  
  //  TLegend *leg = new TLegend(0.11,0.86,0.97,0.99,NULL,"brNDC");
  TLegend *leg = new TLegend(0.11,0.85,0.97,0.98,NULL,"TLNDC");
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);
	leg->SetBorderSize(1);
	leg->SetNColumns(2);
  leg->AddEntry(pvalData,"observed","l");
  leg->AddEntry(discover,"3#sigma discovery prob.","l");
  leg->AddEntry(pval68,"expected for no signal (median + 68% interval) ","fl");
  leg->AddEntry(discover2,"median discovery p-value","l");
  leg->AddEntry(pval90,"expected for no signal (median + 90% interval)","fl");
  //leg->AddEntry((TObject*)0, "(median + 68% prob interval)", "");
  //leg->AddEntry((TObject*)0, "(median + 90% prob interval)", "");
  leg->Draw();
  
}

