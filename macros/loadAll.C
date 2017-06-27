

TChain* loadAll(string date = "06.21.17_23h"){
  //crab
  //kwaabz
  //crab
  //crabcrab
  //kiki
  //popo
  //bert
  //kona
  //tate

  TChain *summaryTree = new TChain("summaryTree");
  
  gROOT->ProcessLine(".x setupProof.C");

  char* resultsDir = getenv("ANITA3_RESULTSDIR");

  stringstream name;
  for (int run=0; run<256; run++) {

    name.str("");
    name << resultsDir << date << "/" << run << ".root";

    summaryTree->Add(name.str().c_str());

  }

    summaryTree->SetProof();

  return summaryTree;

}


TH1* makeCutStrengthPlot(TH1* inHist) {

  TH1* copyHist = (TH1*)inHist->Clone();
  
  copyHist->Scale(1./copyHist->GetIntegral()[copyHist->GetNbinsX()]);

  TH1* outHist = copyHist->GetCumulative();
  delete copyHist;

  return outHist;
}

  
