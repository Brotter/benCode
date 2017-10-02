

TChain* loadAll(string date = "09.27.17_19h/",bool doProof = true){
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
  
  if (doProof) gROOT->ProcessLine(".x setupProof.C");

  char* resultsDir = getenv("ANITA3_RESULTSDIR");

  stringstream name;
  for (int run=0; run<256; run++) {

    //runs that are totally broken
    if (run==220) continue;

    name.str("");
    name << resultsDir << "templateSearch/" << date << "/" << run << ".root";


    cout << "core:" << run << endl;
    summaryTree->Add(name.str().c_str());
    cout << summaryTree->GetEntries() << endl;
  }

  if (doProof) summaryTree->SetProof();

  return summaryTree;

}


TChain* loadWais(string date = "08.04.17_17h"){

  TChain *summaryTree = new TChain("summaryTree");
  

  char* resultsDir = getenv("ANITA3_RESULTSDIR");

  stringstream name;
  for (int core=162; core<182; core++) {

    name.str("");
    name << resultsDir << "templateSearch/" << date << "/" << core << ".root";

    summaryTree->Add(name.str().c_str());

  }

  return summaryTree;

}
