#include "AnitaEventSummary.h"

void plotSomething(string date = "06.05.17_14h"){

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  
  //  string baseDir = "/Volumes/ANITA3data/bigAnalysisFiles/templateSearch/";
  string baseDir = "results/";

  stringstream name;
  for (int run=130; run<433; run++) {
    if ( (run==130) || (run==144) || (run==150) || (run==186) || (run==198) )continue;
    name.str("");
    //    name << "/home/brotter/nfsShared/results/templateSearch/" << date << "/" << run << ".root";
    name << baseDir << date << "/" << run << ".root";
    summaryTree->Add(name.str().c_str());

  }

  AnitaEventSummary *summary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&summary);
  double templateValueV,templateValueH;
  summaryTree->SetBranchAddress("templateValueV",&templateValueV);
  summaryTree->SetBranchAddress("templateValueH",&templateValueH);


  TH2D *hist = new TH2D("hist","hist",100,0,0.2,100,0,1);
  TGraph *wais = new TGraph();
  wais->SetName("wais");
  TGraph *ldb = new TGraph();
  ldb->SetName("ldb");
  
  for (int entry=0; entry<summaryTree->GetEntries(); entry++) {
    if (entry%10000 == 0) cout << entry << "/" << summaryTree->GetEntries() << endl;
    summaryTree->GetEntry(entry);

    if (summary->flags.pulser == 0) {
      hist->Fill(summary->peak[0][0].value,templateValueH);
    }
    else if (summary->flags.pulser == 1) {
      wais->SetPoint(wais->GetN(),summary->peak[0][0].value,templateValueH);
    }
    else if (summary->flags.pulser == 2) {
      ldb->SetPoint(ldb->GetN(),summary->peak[0][0].value,templateValueH);
    }

  }


  hist->Draw("colz");
  wais->SetMarkerColor(kRed);
  if (wais->GetN() != 0) wais->Draw("pSame");
  ldb->SetMarkerColor(kBlue);
  if (ldb->GetN() != 0) ldb->Draw("pSame");

  return;

}


