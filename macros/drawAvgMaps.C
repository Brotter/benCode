#include "AnitaConventions.h"

void getDate(int realTime,char* date) {

  time_t t = realTime;
  struct tm *tm = localtime(&t);
  strftime(date, sizeof(date), "%m.%d  %H:%M", tm);

  return;
}


void drawAvgMaps(string date="07.05.17_22h/") {

  stringstream name;

  char* resultsDir = getenv("ANITA3_RESULTSDIR");

  TChain *summaryTree = new TChain("summaryTree","summaryTree");

  for (int core=0; core<256; core++) {
    name.str("");
    name << resultsDir << date << core << ".root";
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

      noiseSum->avgMapProf[0]->SetStats(0);
      noiseSum->avgMapProf[0]->GetZaxis()->SetRangeUser(0,0.025);
      name.str("");
      char currTime[64];
      getDate(eventSummary->realTime,currTime);
      name << "Average Interferometric Map - " << currTime;
      noiseSum->avgMapProf[0]->SetTitle(name.str().c_str());

      noiseSum->avgMapProf[0]->Draw("colz");

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
