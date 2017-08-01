

TChain* loadAll(string date = "07.28.17_17h/"){
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
    name << resultsDir << "templateSearch/" << date << "/" << run << ".root";

    summaryTree->Add(name.str().c_str());

  }

    summaryTree->SetProof();

  return summaryTree;

}


TChain* loadWais(string date = "07.28.17_17h/"){
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
  for (int core=145; core<181; core++) {

    name.str("");
    name << resultsDir << "templateSearch/" << date << "/" << core << ".root";

    summaryTree->Add(name.str().c_str());

  }

    summaryTree->SetProof();

  return summaryTree;

}
