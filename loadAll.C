

void loadAll(string date = "03.28.17_18h"){

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  
  gROOT->ProcessLine(".x setupProof.C");

  stringstream name;
  for (int run=130; run<213; run++) {
    if ( (run==130) || (run==144) || (run==150) || (run==186) || (run==198) )continue;
    name.str("");
    //    name << "/home/brotter/nfsShared/results/templateSearch/" << date << "/" << run << ".root";
    name << "/Volumes/ANITA3data/bigAnalysisFiles/templateSearch/" << date << "/" << run << ".root";
    summaryTree->Add(name.str().c_str());

  }

  summaryTree->SetProof();

  return;

}
