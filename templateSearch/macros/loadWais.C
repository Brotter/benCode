


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

void separateWais() {

  TChain *summaryTree = loadWais();

  int lenEntries = summaryTree->GetEntries();
  cout << lenEntries << " total entries found" << endl;

  AnitaEventSummary *eventSummary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&eventSummary);

  AnitaTemplateSummary *templateSummary = NULL;
  summaryTree->SetBranchAddress("template",&templateSummary);

  AnitaNoiseSummary *noiseSummary = NULL;
  summaryTree->SetBranchAddress("noiseSummary",&noiseSummary);


  TFile *outFile = TFile::Open("waisEvents.root","recreate");
  
  TTree *waisTree = new TTree("waisSummary","waisSummary");
  waisTree->Branch("eventSummary",&eventSummary);
  waisTree->Branch("template",&templateSummary);
  waisTree->Branch("noiseSummary",&noiseSummary);
  
  int waisCnt = 0;
  
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%1000 == 0) {
      cout << entry << "/" << lenEntries;
    }

    summaryTree->GetEntry(entry);

    if (eventSummary->flags.pulser == 1) {
      outFile->cd();
      waisTree->Fill();
      waisCnt++;
    }
  }
  
  outFile->cd();
  waisTree->Write();
  outFile->Close();


  cout << "Done!" << endl;

  return;

}
