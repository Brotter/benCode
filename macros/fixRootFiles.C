void CopyDir(TDirectory *source) {
  //copy all objects and subdirs of directory source as a subdir of the current directory
  source->ls();
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir->mkdir(source->GetName());
  adir->cd();
  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir);
      adir->cd();
    } else if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)source->Get(key->GetName());
      adir->cd();
      TTree *newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      adir->cd();
      obj->Write();
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}

void fixRootFiles(string date="06.28.17_15h") {

  stringstream name;
  char* resultsDir = getenv("ANITA3_RESULTSDIR");


  for (int core=0; core<256; core++) {
    cout << core << endl;
    name.str("");
    name << resultsDir << date << "/" << core << ".root";
    TFile *inFile = TFile::Open(name.str().c_str());

    if (!inFile || inFile->IsZombie()) {
      cout << "File " << name.str() << " doesn't exist!" << endl;
      continue;
    }

    name.str("");
    name << resultsDir << date << "/" << core << "_fixed.root";
    TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
    
    outFile->cd();
    CopyDir(inFile);
    inFile->Close();
    outFile->Close();
  }


  return;
}
