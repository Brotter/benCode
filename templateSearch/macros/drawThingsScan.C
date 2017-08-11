/*

  plotThingsScan got overwhelmingly long, so for drawing the created histograms come here!

*/




TH1* makeNormCumulative(TH1* inHist) {
  /*
    Makes a "normalized cumulative cut fraction" graph I think

    Should let me set cuts on a specific amount of "reduction"
   */


  TH1* copyHist = (TH1*)inHist->Clone();
  
  double integral = 0;
  for (int i=0; i<copyHist->GetNbinsX(); i++) {
    double value = copyHist->GetBinContent(i);
    integral += value;
  }

  copyHist->Scale(1./integral);
  TH1* outHist = (TH1*)copyHist->GetCumulative();

  for (int i=0; i<copyHist->GetNbinsX(); i++) {
    double value = outHist->GetBinContent(i);
    outHist->SetBinContent(i,1.-value);
  }


  delete copyHist;

  outHist->GetYaxis()->SetTitle("Surviving Fraction");
  return outHist;
}


TH1* makeNormHist(TH1* inHist) {
  /*

    Make a normalized histogram

  */

  TH1* copyHist = (TH1*)inHist->Clone();
  copyHist->GetYaxis()->SetTitle("Occupancy");

  
  double integral = 0;
  for (int i=0; i<copyHist->GetNbinsX(); i++) {
    double value = copyHist->GetBinContent(i);
    integral += value;
  }

  copyHist->Scale(1./integral);

  return copyHist;
}

void drawThingsScan(string inFileName="plotThingsScan.root") {
  
  stringstream name;

  TFile *inFile = TFile::Open(inFileName.c_str());

  const int numLines = 5;
  string subStrings[numLines] = {"thermal","WAIS","LDB","basePointed","minbias"};

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);


  const int numGraphs = 31;
  string graphNames[numGraphs] = {"peakPhi","peakTheta","mapPeak","mapSNR","peakRatio",
			   "linPolFrac","linPolAng","waveSNR","wavePeakVal","wavePeakHilb","impulsivity",
			   "linPolFrac_F","linPolAng_F","waveSNR_F","wavePeakVal_F","wavePeakHilb_F","impulsivity_F",
			   "linPolFrac_D","linPolAng_D","waveSNR_D","wavePeakVal_D","wavePeakHilb_D","impulsivity_D",
			   "linPolFrac_DF","linPolAng_DF","waveSNR_DF","wavePeakVal_DF","wavePeakHilb_DF","impulsivity_DF",
			   "template_Wais","template_cRay"};

  for (int graph=0; graph<numGraphs; graph++) {
    c1->Clear();
    c1->SetGrid(1);
    c1->SetLogy();
    TLegend *leg = new TLegend(0.8,0.7,0.95,0.95);
    for (int i=0; i<numLines; i++) {
      name.str("");
      name << graphNames[graph] << "_" << subStrings[i];
      TH1 *hist = makeNormCumulative((TH1D*)inFile->Get(name.str().c_str()));
      hist->SetLineColor(i+1);
      hist->GetYaxis()->SetRangeUser(1e-7,1);
      hist->SetStats(0);
      if (i==0) hist->Draw();
      else      hist->Draw("same");
      if (strcmp(subStrings[i].c_str(),"thermal")==0) leg->AddEntry(hist,"events","l");
      else leg->AddEntry(hist,subStrings[i].c_str(),"l");
    }

    leg->Draw();
    c1->SaveAs(TString("plots/")+graphNames[graph]+TString("_cumulative.png"));
  }


  for (int graph=0; graph<numGraphs; graph++) {
    c1->Clear();
    c1->SetGrid(1);
    c1->SetLogy();
    TLegend *leg = new TLegend(0.8,0.7,0.95,0.95);
    for (int i=0; i<numLines; i++) {
      name.str("");
      name << graphNames[graph] << "_" << subStrings[i];
      TH1 *hist = makeNormHist((TH1D*)inFile->Get(name.str().c_str()));
      hist->SetLineColor(i+1);
      hist->GetYaxis()->SetRangeUser(1e-8,1);
      hist->SetStats(0);
      if (i==0) hist->Draw();
      else      hist->Draw("same");
      if (strcmp(subStrings[i].c_str(),"thermal")==0) leg->AddEntry(hist,"events","l");
      else leg->AddEntry(hist,subStrings[i].c_str(),"l");
    }

    leg->Draw();
    c1->SaveAs(TString("plots/")+graphNames[graph]+TString("_hist.png"));
  }

  return;
}
  

TH1 *dividedCumulative(TH1D *h1, TH1D *h2) {

  TH1* h1Cum = makeNormCumulative(h1);
  TH1* h2Cum = makeNormCumulative(h2);

  h1Cum->Divide(h2Cum);
  delete h2Cum;
  
  return h1Cum;
}
