#include "AnitaConventions.h"

void drawAvgMaps() {

  stringstream name;

  TChain *summaryTree = new TChain("summaryTree","summaryTree");

  for (int core=0; core<256; core++) {
    name.str("");
    name << core << ".root";
    summaryTree->Add(name.str().c_str());
  }

  AnitaNoiseSummary *noiseSum = NULL;
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);

  AnitaEventSummary *eventSummary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&eventSummary);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  TGraph *gSun = new TGraph();
  gSun->SetMarkerStyle(kStar);


  int cnt = 0;

  for (int entry=0; entry<summaryTree->GetEntries(); entry++) {
    summaryTree->GetEntry(entry);
    if (noiseSum->isMinBias && noiseSum->mapFifoFillFlag) {
      cnt++;
      c1->cd();
      c1->Clear();
      noiseSum->avgMap[0]->SetStats(0);
      noiseSum->avgMap[0]->Draw("colz");

      double sunTheta = eventSummary->sun.theta;
      double sunPhi = eventSummary->sun.phi;
      if (sunPhi < 0) sunPhi += 360;
      //      cout << sunPhi << " " << sunTheta << endl;
      gSun->SetPoint(0,sunPhi,-sunTheta);
      gSun->Draw("pSame");
      

      name.str("");
      name << "avgMaps/" << setfill('0') << setw(3) << cnt << ".png";
      c1->SaveAs(name.str().c_str());
    }
  }
    


  return;
}
