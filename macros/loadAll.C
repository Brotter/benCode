

TChain* loadAll(string date = "06.18.17_15h"){

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  
  gROOT->ProcessLine(".x setupProof.C");

  char* resultsDir = getenv("ANITA3_RESULTSDIR");

  stringstream name;
  for (int run=130; run<433; run++) {
    if ( (run==130) || (run==144) || (run==150) || (run==186) || (run==198) )continue;
    name.str("");
    name << resultsDir << date << "/" << run << ".root";

    summaryTree->Add(name.str().c_str());

  }

  summaryTree->SetProof();

  return summaryTree;

}


