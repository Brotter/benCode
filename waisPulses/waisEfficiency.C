/*


  Lets calculate the efficiency of wais stuff as a function of minute or whatever


 */


void waisEfficiency() {
  

  //set up input
  TFile *inFile  = TFile::Open(TString(getenv("HOME"))+"/anita16/benCode/templateSearch/macros/waisEvents.root");
  if (!inFile->IsOpen()) {
    cout << "Couldn't find file" << endl;
    return;
  }
  TTree* summaryTree = (TTree*)inFile->Get("waisSummary");
  if (!summaryTree) {
    cout << "Couldn't find tree" << endl;
    return;
  }
  int lenEntries = summaryTree->GetEntries();
  cout << "Opened tree with " << lenEntries << " entries" << endl;

  AnitaEventSummary *evSummary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSummary);


  //set up output
  TFile *outFile = TFile::Open("waisEfficiency.root","recreate");
  TTree *outTree = new TTree("waisEfficiency","waisEfficiency");
  int count = 0;
  bool fWaisEffFilled = false;
  outTree->Branch("count",&count);
  outTree->Branch("fWaisEffFilled",&fWaisEffFilled);


  summaryTree->GetEntry(0);
  long firstTime = evSummary->realTime;
  summaryTree->GetEntry(lenEntries-1);
  long finalTime = evSummary->realTime;
  int numSeconds = finalTime-firstTime;
  summaryTree->GetEntry(0);

  const int timeWindow = 1000; //seconds
  std::deque<bool> dCount(timeWindow,false);
  std::queue<bool> qCount(dCount);

  int entry = 0;

  for (int second=0 ;second<numSeconds; second++) {
    if (second%1000 == 0) cout << second << "/" << numSeconds << endl;

    if (second >= timeWindow) fWaisEffFilled = true;
    
    while (evSummary->realTime < second+firstTime) {
      outTree->Fill();
      entry++;
      summaryTree->GetEntry(entry);
    }
    
    if (evSummary->realTime == second+firstTime) {
      count++;
      qCount.push(true);
    }
    else {
      qCount.push(false);
    }

    if (qCount.front()) count--;
    qCount.pop();

  }


  outTree->Write();
  outFile->Close();
  
  return;
}
