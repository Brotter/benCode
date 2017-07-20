{

  TChain *headTree = new TChain("headTree","headTree");

  stringstream name;
  for (int i=150; i<440; i++) {
    //these runs are missing in outputs
    if ((i>=257 && i<=263) || (i==368) || (i==393)) continue;
    name.str("");
    name << "/Volumes/ANITA3Data/root/run" << i << "/headFile" << i << ".root";
    headTree->Add(name.str().c_str());
    
  }    
  
  cout << headTree->GetEntries() << " total events (in headTree)" << endl;

    
  return;
}
