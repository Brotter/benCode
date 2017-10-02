/*

  The redo runs are split into 12, this should combine them

 */



void combineRedos(int core) {

  TChain *summaryTree = new TChain("summaryTree");

  stringstream name;
  name.str("");
  char* homeDir = getenv("HOME");
  string dataDir = "/anita16/benCode/templateSearch/09.27.17_19h/";
  for (int i=0; i<12; i++) {
    name.str("");
    name << homeDir << dataDir << core << "_" << i+1 << ".root";
    summaryTree->Add(name.str().c_str());
  }

  name.str("");
  name << homeDir << dataDir << core << ".root";
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
  summaryTree->CloneTree(-1,"fast");
  outFile->Write();
  outFile->Close();


  return;
}

void combineRedos() {
  cout << "combineRedos(int core): You didn't give me a core number!  Loading instead..." << endl;
  return;
}
