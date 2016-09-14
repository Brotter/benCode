{

  TChain *newCorrelatorTree = new TChain("newCorrelatorTree","newCorrelatorTree");

  stringstream name;
  for (int i=128; i<192; i++) {
    name.str("");
    name << "rootFiles/waisNewCorrelator_" << i << ".root";
    newCorrelatorTree->Add(name.str().c_str());
  }


}
