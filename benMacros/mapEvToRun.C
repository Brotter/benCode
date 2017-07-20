//for magic display, it is important to be able to map the event number to run (if you don't have the run)



void mapEvToRun(){

  char* dataDir = getenv("ANITA3_DATA");
  stringstream name;
  for (int run=130; run<440; run++) {
    if (run >= 257 && run <= 263) continue; //missing runs

    name.str("");
    name << dataDir << "run" << run << "/headFile" << run << ".root";
    TFile *headFile = TFile::Open(name.str().c_str());
    if (headFile == 0) continue;
    TTree *headTree = (TTree*)headFile->Get("headTree");
    RawAnitaHeader *head = NULL;
    headTree->SetBranchAddress("header",&head);

    int lenEntries = headTree->GetEntries();
    headTree->GetEntry(0);
    int firstEntry = head->eventNumber;
    headTree->GetEntry(lenEntries-1);
    int lastEntry = head->eventNumber;

    cout << run << " " << firstEntry << " " << lastEntry << endl;

    headFile->Close();

  }

  return;


}
