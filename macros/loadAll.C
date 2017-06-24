

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

