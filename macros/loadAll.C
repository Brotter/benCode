

TChain* loadAll(string date = "06.18.17_22h"){

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  
  gROOT->ProcessLine(".x setupProof.C");

  char* resultsDir = getenv("ANITA3_RESULTSDIR");

  stringstream name;
  for (int run=130; run<440; run++) {
    //dead runs for whatever reason
    if ( ( run >= 214 && run <=222) || (run >= 257 && run<=263) ) continue;

    name.str("");
    name << resultsDir << date << "/" << run << ".root";

    summaryTree->Add(name.str().c_str());

  }

  summaryTree->SetProof();

  return summaryTree;

}


