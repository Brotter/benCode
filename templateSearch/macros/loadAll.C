
TChain *loadWhatever(string basePrefix, int numFiles, bool doProof=false) {
  /*

    A generalized form of loading things

  */


  TChain *summaryTree = new TChain("summaryTree","summaryTree");

  if (doProof) gROOT->ProcessLine(".x setupProof.C");

  stringstream name;
    for (int i=0; i<numFiles; i++) {
    name.str("");
    name << basePrefix << "_" << i << ".root";
    summaryTree->Add(name.str().c_str());
  }
  int lenEntries = summaryTree->GetEntries();
  cout << "loadWhatever(): Found " << lenEntries << " entries" << endl;

  if (doProof) summaryTree->SetProof();

  return summaryTree;
}



TChain *loadReKey(bool doProof = true) {
  /*
    
    load the reprocessed key values for the final set

  */


  cout << "loadReKey(): Loading recalculated key values ";
  if (doProof) cout << " with PROOF enabled" << endl;
  else         cout << " without PROOF" << endl;

  TChain *summaryTree = new TChain("summaryTree");
  
  if (doProof) gROOT->ProcessLine(".x setupProof.C");

  char* resultsDir = getenv("ANITA3_RESULTSDIR");

  stringstream name;
  for (int run=0; run<64; run++) {
    name.str("");
    name << resultsDir << "templateSearch/09.27.17_19h/reCalcKeyValues/goodEvents_corr_" << run << ".root";
    summaryTree->Add(name.str().c_str());
  }
  cout << "loadReKey(): Found " << summaryTree->GetEntries() << " entries" << endl;

  if (doProof) summaryTree->SetProof();

  return summaryTree;

}

TChain *loadLabeled() {
  TChain *reKey = loadReKey(false);

  TChain *labeled = new TChain("summaryTree","summaryTree");
  labeled->Add("labelEvents_All.root");
  
  reKey->AddFriend(labeled);

  return reKey;
}
  

TChain *loadGeoAssociated() {
  TChain *outTree = new TChain("summaryTree","summaryTree");
  outTree->Add("geoAssociated/geoAssociated_ev9097075_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev11116669_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev11989349_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev16952229_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev19459851_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev23695286_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev27142546_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev32907848_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev33484995_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev39599205_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev41529195_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev58592863_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev62273732_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev66313844_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev68298837_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev70013898_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev73726742_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev75277769_geoMag.root");
  outTree->Add("geoAssociated/geoAssociated_ev83877990_geoMag.root");

  cout << "Found " << outTree->GetEntries() << " entries" << endl;
  return outTree;
}

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

    Generating the list of events that geoassociate with the impulsive events requires the cluster
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



TChain *loadBkgdCluster(string date="10.31.17_11h59m", int numFiles=64, bool doProof=false) {
  /*
    Just a quick thing I call all the time
   */

  stringstream name;
  char* dataDir = getenv("ANITA3_RESULTSDIR");
  name.str("");
  name << dataDir << "/cluster/" << date << "/clusterBackground";
  return loadWhatever(name.str(),numFiles,doProof);
}



  


/*  Default macro if called from bash command line */
void loadAll() {
  cout << "loaded loadAll.C" << endl;
  return;
}
