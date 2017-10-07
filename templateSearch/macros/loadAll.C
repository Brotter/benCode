

TChain* loadAll(string date,bool doProof = true){
  //crab
  //kwaabz
  //crab
  //crabcrab
  //kiki
  //popo
  //bert
  //kona
  //tate

  cout << "loadAll(): Loading " << date;
  if (doProof) cout << " with PROOF enabled" << endl;
  else         cout << " without PROOF" << endl;

  TChain *summaryTree = new TChain("summaryTree");
  
  if (doProof) gROOT->ProcessLine(".x setupProof.C");

  char* resultsDir = getenv("ANITA3_RESULTSDIR");

  stringstream name;
  for (int run=0; run<256; run++) {
    name.str("");
    name << resultsDir << "templateSearch/" << date << "/" << run << ".root";
    summaryTree->Add(name.str().c_str());
  }
  cout << "loadAll(): Found " << summaryTree->GetEntries() << " entries" << endl;

  if (doProof) summaryTree->SetProof();

  return summaryTree;

}

TChain *loadAllDefault() { return loadAll("09.27.17_19h",true); }
TChain *loadAllDefault_noproof() { return loadAll("09.27.17_19h",false); }


TChain *loadPseudoBases(string date = "10.05.17_14h",bool doProof=false) {

  /*

    Generating the list of events that cluster with the impulsive events requires the cluster
    so this imports them all (once you save them to disk!)

  */
  cout << "loadPseudoBases(): Loading " << date;
  if (doProof) cout << " with PROOF enabled" << endl;
  else         cout << " without PROOF" << endl;

  TChain *summaryTree = new TChain("summaryTree","summaryTree");

  if (doProof) gROOT->ProcessLine(".x setupProof.C");

  stringstream name;
  string resultsDir = getenv("ANITA3_RESULTSDIR");
    for (int i=0; i<192; i++) {
    name.str("");
    name << resultsDir << "cluster/" << date << "/pseudoBaseCluster_" << i << ".root";
    summaryTree->Add(name.str().c_str());
  }
  int lenEntries = summaryTree->GetEntries();
  cout << "loadPseudoBases(): Found " << lenEntries << " entries" << endl;

  if (doProof) summaryTree->SetProof();

  return summaryTree;
}



/*  Default macro if called from bash command line */
void loadAll() {
  cout << "loaded loadAll.C" << endl;
  return;
}
