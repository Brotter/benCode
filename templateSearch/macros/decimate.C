/*

  Going through the 100% dataset takes a really long time, so this is a method for decimating it I think

  Super simple, just event numbers that end in zero

 */







void decimate(int downfactor=10,string date="07.28.17_17h") {



  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  char* resultsDir = getenv("ANITA3_RESULTSDIR");
  stringstream name;
  for (int core=0; core<256; core++) {
    name.str("");
    name << resultsDir << "templateSearch/" << date << "/" << core << ".root";
    summaryTree->Add(name.str().c_str());
  }

  
  TFile *outFile = new TFile(TString(date)+"_decimated.root","recreate");
  TTree *outTree = new TTree("summaryTree","summaryTree");

  AnitaEventSummary *evSummary = NULL;
  AnitaTemplateSummary *tmpSummary = NULL;
  AnitaNoiseSummary *noiseSummary = NULL;
  Adu5Pat *gps = NULL;

  summaryTree->SetBranchAddress("eventSummary",&evSummary);
  summaryTree->SetBranchAddress("template",&tmpSummary);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSummary);
  summaryTree->SetBranchAddress("gpsEvent",&gps);

  outTree->Branch("eventSummary",&evSummary);
  outTree->Branch("template",&tmpSummary);
  outTree->Branch("noiseSummary",&noiseSummary);
  outTree->Branch("gpsEvent",&gps);

  int lenEntries = summaryTree->GetEntries();

  cout << "found " << lenEntries << " entries" << endl;

  int cnt=0;
  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    
    if (entry%10000 == 0) cout << entry << "/" << lenEntries << " (" << cnt << ")" << endl;

    if (entry%downfactor == 0) {
      outTree->Fill();
      cnt++;
    }

  }

  outTree->Write();
  outFile->Close();

  return;

}
