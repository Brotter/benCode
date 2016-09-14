{

  TChain *newCorrelatorTree = new TChain("newCorrelatorTree","newCorrelatorTree");

  stringstream name;
  for (int i=0; i<256; i++) {
    name.str("");
    name << "rootFiles/waisNewCorrelator_" << i << ".root";
    newCorrelatorTree->Add(name.str().c_str());
  }


  cout << "number of entries: " << newCorrelatorTree->GetEntries() << endl;

  newCorrelatorTree->Draw("peakThetaDeg:peakPhiDeg-heading","","colz");

  return;

}
