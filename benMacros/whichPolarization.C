void whichPolarization(int runNum) {

  stringstream name;

  TChain *headTree = new TChain("headTree","headTree");
  name.str("");
  name << "/Volumes/ANITA3Data/root/run" << runNum << "/headFile" << runNum << ".root";
  headTree->Add(name.str().c_str());

  RawAnitaHeader *head = NULL;
  headTree->SetBranchAddress("header",&head);

  int numEntries = headTree->GetEntries();

  TChain *tdataEvent = new TChain("tdataEvent","tdataEvent");
  name.str("");
  name << "/Volumes/ANITA3Data/rootOutputs/output" << runNum <<"_1.root";
  tdataEvent->Add(name.str().c_str());
  tdataEvent->BuildIndex("eventNumber");


  int numBoth = 0;
  int numVertical = 0;
  int numHorizontal = 0;
  int noMatch = 0;
  int noTrig = 0;


  int index = -1;
  for (int entry=0; entry<numEntries; entry++) {
    headTree->GetEntry(entry);
    index = tdataEvent->GetEntryNumberWithIndex(head->eventNumber);
    if (index == -1) {
      noMatch++; }
    else {
      if (head->l3TrigPattern) {
	if (head->l3TrigPatternH) numBoth++;
	else numVertical++; }
      else if (head->l3TrigPatternH) numHorizontal++;      
      else noTrig++; 
    }
  }

  cout << "total Entries in head file: " << numEntries << endl;
  cout << "number of Veritcal trigs: " << numVertical << endl;
  cout << "number of Horizontal trigs: " << numHorizontal << endl;
  cout << "number of Both triggered: " << numBoth << endl;
  cout << "number with no trig flags: " << noTrig << endl;
  cout << "number with no eventNum match: " << noMatch << endl;
  cout << "total: " << numVertical+numHorizontal+numBoth+noTrig+noMatch << endl;


  return;

}
      
