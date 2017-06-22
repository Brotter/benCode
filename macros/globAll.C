/*

  The cluster returns a bunch of (possibly incorrectly closed) files that I want in ONE GIANT FILE >:D

 */

void globAll(string date = "06.11.17_19h"){

  string basedir = "/home/brotter/nfsShared/results/templateSearch/";
    //string basedir = "/Volumes/ANITA3data/bigAnalysisFiles/templateSearch/";

  stringstream name,outFile;


  outFile.str("");
  //  outFile << basedir << date << "/all.root";
  outFile << "all.root";


  TChain *summaryTree = new TChain("summaryTree","summaryTree");

  for (int run=130; run<433; run++) {
    //    if ( (run==130) || (run==144) || (run==150) || (run==186) || (run==198) )continue;
    name.str("");
    name << basedir << date << "/" << run << ".root";
    
    summaryTree->Add(name.str().c_str());

    cout << name.str() << endl;

  }

  summaryTree->Merge(outFile.str().c_str());
  
  return;


}
