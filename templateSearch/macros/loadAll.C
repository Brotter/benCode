

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


TChain *loadPseudoBases(string date = "10.03.17_22h",bool doProof=false) {

  /*

    Generating the list of events that cluster with the impulsive events requires the cluster
    so this imports them all (once you save them to disk!)

  */

  TChain *summaryTree = new TChain("summaryTree","summaryTree");

  if (doProof) gROOT->ProcessLine(".x setupProof.C");

  stringstream name;
  string resultsDir = getenv("ANITA3_RESULTSDIR");
    for (int i=0; i<64; i++) {
    name.str("");
    name << resultsDir << "cluster/" << date << "/pseudoBaseCluster_" << i << ".root";
    summaryTree->Add(name.str().c_str());
  }
  int lenEntries = summaryTree->GetEntries();
  cout << "loadPseudoBases(): Found " << lenEntries << " entries" << endl;

  if (doProof) summaryTree->SetProof();

  return summaryTree;
}
