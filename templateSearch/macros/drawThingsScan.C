/*

  plotThingsScan got overwhelmingly long, so for drawing the created histograms come here!

*/


TH1* makeNormCumulative(TH1* inHist,Bool_t forward=true) {
  /*
    Makes a "normalized cumulative cut fraction" graph I think

    Should let me set cuts on a specific amount of "reduction"

    forward specifies the direction of the cumulative sum
   */


  TH1* copyHist = (TH1*)inHist->Clone();
  
  double integral = 0;
  for (int i=0; i<copyHist->GetNbinsX(); i++) {
    double value = copyHist->GetBinContent(i);
    integral += value;
  }

  copyHist->Scale(1./integral);
  TH1* outHist = (TH1*)copyHist->GetCumulative(forward);

  for (int i=0; i<copyHist->GetNbinsX(); i++) {
    double value = outHist->GetBinContent(i);
    outHist->SetBinContent(i,1.-value);

  }
  delete copyHist;

  if (forward) outHist->GetYaxis()->SetTitle("Surviving Fraction");
  else         outHist->GetYaxis()->SetTitle("Fraction of Events Cut");
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

void saveAllImages(string inFileName="plotThingsScan.root") {
  
  stringstream name;

  TFile *inFile = TFile::Open(inFileName.c_str());

  const int numLines = 5;
  string subStrings[numLines] = {"events","WAIS","LDB","basePointed","minbias"};

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

void getCutsFromValue(double waisCutFraction = 1e-3) {
  
  TFile *inFile = TFile::Open("plotThingsScan.root");
  if (!inFile->IsOpen()) {
    cout << "Couldn't open file" << endl;
    return;
  }

  const int numInteresting = 4;
  string interestingValues[numInteresting] = {"template_Wais","wavePeakHilb_DF","mapPeak","mapSNR"};

  for (int value=0; value<numInteresting; value++) {
    string name;
    name = interestingValues[value] + "_WAIS";
    TH1D *hist = (TH1D*)inFile->Get(name.c_str());
    TH1 *histCut = makeNormCumulative(hist,false);

    name = interestingValues[value]+"_minbias";
    TH1D *minBiasHist = (TH1D*)inFile->Get(name.c_str());
    TH1* minBiasRemain = makeNormCumulative(minBiasHist);

    double cutLoc = 0;
    double minBiasValue = 0;
    for (int bin=0; bin<histCut->GetNbinsX(); bin++) {
      if (histCut->GetBinContent(bin+1) > waisCutFraction) {
	cutLoc = histCut->GetBinCenter(bin+1);
	minBiasValue = minBiasRemain->GetBinContent(bin+1);
	break;
    }
  }

    cout << interestingValues[value] << " "  << cutLoc << " " << minBiasValue << endl;
  }


  return;
}



void drawTwoReversedHists(TH1* inHist1, TH1* inHist2, string saveName="") {

  TH1* drawHist1 = makeNormCumulative(inHist1);
  TH1* drawHist2 = makeNormCumulative(inHist2,false);
  drawHist2->SetTitle("");

  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry(drawHist1,"Possible Signal","l");
  leg->AddEntry(drawHist2,"WAIS Pulses","l");

  drawHist1->SetStats(false);
  drawHist2->SetStats(false);

  drawHist1->GetYaxis()->SetRangeUser(1e-7,1);
  drawHist2->GetYaxis()->SetRangeUser(1e-7,1);

  drawHist2->GetYaxis()->SetTitleColor(kRed);
  drawHist2->SetLineColor(kRed);


  TCanvas *c1 = new TCanvas("c1","transparent pad",1000,600);
  TPad *p1 = new TPad("p1","",0,0,1,1);
  p1->SetLogy();
  TPad *p2 = new TPad("p2","",0,0,1,1);
  p2->SetLogy();
  p2->SetFillStyle(4000);
  p2->SetFillColor(0);
  p2->SetFrameFillStyle(0);
  p1->Draw();
  p1->cd();
  drawHist1->Draw();
  p2->Draw();
  p2->cd();
  drawHist2->Draw("Y+");
  leg->Draw("");

  if (saveName != "") c1->SaveAs(saveName.c_str());

  return;

}
