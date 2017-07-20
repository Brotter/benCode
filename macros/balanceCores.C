/*

  I need to generate some lookup table for balancing the cores.

  78,630,542 total events from run 130 to 440

  307,151 events for each of the 256 cores

  figure out which runs those are in

 */



void balanceCores() {

  char* dataDir = getenv("ANITA3_DATA");

  const int numCores = 256;

  stringstream name;
  
  TChain *headTree = new TChain("headTree","headtree");

  for (int run=130; run<440; run++) {
    name.str("");
    name << dataDir << "run" << run << "/headFile" << run << ".root";
    headTree->Add(name.str().c_str());

  }

  RawAnitaHeader *head = NULL;
  headTree->SetBranchAddress("header",&head);

  int lenEntries = headTree->GetEntries();
  int evPerCore = (lenEntries/numCores) + 1;

  cout << "Found " << lenEntries << " total entries, " << evPerCore << " for each of the " << numCores << " cores"  << endl;


  cout << "#core startEvNum startRun endEvNum endRun" << endl;
  for (int core=0; core<numCores; core++) {
    int startEntry = core*evPerCore;
    int endEntry = (core+1)*evPerCore - 1;

    headTree->GetEntry(startEntry);
    int startEventNumber = head->eventNumber;
    int startRun = head->run;

    headTree->GetEntry(endEntry);
    int endEventNumber = head->eventNumber;
    int endRun = head->run;

    cout << core << " " << startEventNumber << " " << startRun << " " << endEventNumber << " " << endRun;

    if (startRun != endRun) cout << " **** ";
    
    cout << endl;
  }

  TH1D *hEvsPerRun = new TH1D("hEvsPerRun","Events Per Run",400,0,50);

  ofstream outFile("entriesPerRun.txt");

  int run=130;
  int lastStart = 0;
  outFile << run << " 0";
  for (int entry=0; entry<lenEntries; entry++) {
    headTree->GetEntry(entry);
    if (head->run != run) {
      run = head->run;
      outFile << " " << entry-1 << " " << entry - lastStart << endl;
      outFile << run << " " << entry;

      hEvsPerRun->Fill((entry - lastStart)/(3.5*3600));

      lastStart = entry;
    }
  }

  outFile << lenEntries-1 << " " << lenEntries - lastStart << endl;
  hEvsPerRun->Fill((lenEntries - lastStart)/(3.5*3600));
      
  hEvsPerRun->Draw();


  outFile.close();
    

    
    
    
  return;

}

      
