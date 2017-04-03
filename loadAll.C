

void loadAll(string date = "03.28.17_18h"){

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  

  stringstream name;
  for (int run=130; run<213; run++) {
    name.str("");
    name << "/home/brotter/nfsShared/results/templateSearch/" << date << "/" << run << ".root";
    summaryTree->Add(name.str().c_str());

  }

  return;

}
