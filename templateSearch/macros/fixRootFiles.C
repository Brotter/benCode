void CopyDir(TDirectory *source) {
  //copy all objects and subdirs of directory source as a subdir of the current directory
  source->ls();
  TDirectory *savdir = gDirectory;
  cout << "CopyDir(): " << source->GetName() << " " << savdir->GetName() << endl;
  
  savdir->cd();
  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {    
    const char *classname = key->GetClassName();
    cout << "classname:" << classname << endl;
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      savdir->cd();
      CopyDir(subdir);
      savdir->cd();
    } else if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)source->Get(key->GetName());
      savdir->cd();
      TTree *newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      savdir->cd();
      obj->Write();
      delete obj;
    }
  }
  savdir->SaveSelf(kTRUE);
  savdir->cd();
}

void fixRootFiles(string date="06.28.17_15h",int start=0, int end=256) {

  stringstream inName,outName;
  char* resultsDir = getenv("ANITA3_RESULTSDIR");


  for (int core=start; core<end; core++) {
    cout << core << endl;


    inName.str("");
    inName << resultsDir << date << "/" << core << ".root";
    cout << "inFile = " << inName.str() << endl;
    TFile *inFile = TFile::Open(inName.str().c_str());
    if (!inFile || inFile->IsZombie()) {
      cout << "File " << inName.str() << " doesn't exist!" << endl;
      continue;
    }
    cout << "got inFile" << endl;


    outName.str("");
    outName << resultsDir << date << "/" << core << "_fixed.root";
    cout << "outFile = " << outName.str() << endl;
    TFile *outFile = TFile::Open(outName.str().c_str(),"recreate");
    if (!outFile || outFile->IsZombie()) {
      cout << "File " << outName.str() << " didn't open!" << endl;
      continue;
    }
    cout << "got outFile" << endl;


    outFile->cd();
    TDirectory *target = gDirectory;
    CopyDir(inFile);
    inFile->Close();
    outFile->Close();

    remove(inName.str().c_str());
    rename(outName.str().c_str(),inName.str().c_str());

  }


  return;
}
