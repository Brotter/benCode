/*
  Should read in a bunch of timefresnel-root.dat files and then maybe plot them or 
  find interesting things about them


  Ben Rotter - Feb 2017 - University of Hawaii at Manoa
 */


void ZHAiresParser(int energy=195) {


  //  string baseDir = "/storageb/harmscho/ANITA_CR_PAPER/";
  //  string subDir = "exposure_events/";

  string baseDir = "/Volumes/ANITA3Data/ZHAires_results/";

  string subDir = "measured_events/Event-649637/";

  stringstream fullDir;
  fullDir.str("");
  fullDir << baseDir << subDir << "En" << energy << "_off";

  TGraph *outGraphX;
  TGraph *outGraphY;
  TGraph *outGraphZ;

  stringstream filename,name;
  filename.str("");
  filename << "ZHAiresParser_En" << energy << ".root";
  TFile *outFile = TFile::Open(filename.str().c_str(),"recreate");


  for (int off=0; off<4; off++) {
    for (int sub=0; sub<4; sub++) {
      cout << "off" << off << " sub" << sub << endl;
      filename.str("");
      filename << fullDir.str() << off+1 << "_sub" << sub+1 << "/timefresnel-root.dat";


      TTree *inTree = new TTree("inTree","inTree");
      inTree->ReadFile(filename.str().c_str(),"showerNum:antNum:antX:antY:antZ:time:A:Ax:Ay:Az:E:Ex:Ey:Ez");
      Float_t Ex,Ey,Ez,antx,anty,antz,antNum,time;
      inTree->SetBranchAddress("Ex",&Ex);
      inTree->SetBranchAddress("Ey",&Ey);
      inTree->SetBranchAddress("Ez",&Ez);
      inTree->SetBranchAddress("time",&time);
      inTree->SetBranchAddress("antX",&antx);
      inTree->SetBranchAddress("antY",&anty);
      inTree->SetBranchAddress("antZ",&antz);
      inTree->SetBranchAddress("antNum",&antNum);
      int prevAntNum = -1;


      for (int entry=0; entry<inTree->GetEntries(); entry++) {
	inTree->GetEntry(entry);

	if (prevAntNum != antNum) { 
	  cout << "antNum " << antNum << endl;
	  prevAntNum = antNum;

	  if (antNum != 1) {
	    outFile->cd();
	    outGraphX->Write();
	    outGraphY->Write();
	    outGraphZ->Write();
	    
	    delete outGraphX;
	    delete outGraphY;
	    delete outGraphZ;
	  }

	  outGraphX = new TGraph();
	  name.str("");
	  name << "off" << off << "sub" << sub << "ant" << antNum << "_X";
	  outGraphX->SetTitle(name.str().c_str());
	  outGraphX->SetName(name.str().c_str());
	  
	  outGraphY = new TGraph();
	  name.str("");
	  name << "off" << off << "sub" << sub << "ant" << antNum << "_Y";
	  outGraphY->SetTitle(name.str().c_str());
	  outGraphY->SetName(name.str().c_str());
	  
	  outGraphZ = new TGraph();
	  name.str("");
	  name << "off" << off << "sub" << sub << "ant" << antNum << "_Z";
	  outGraphZ->SetTitle(name.str().c_str());
	  outGraphZ->SetName(name.str().c_str());
	  
	  
	}



	int pt = entry-20000*(antNum-1);
	if (pt >= 2000) continue;

	outGraphX->SetPoint(pt,pt*0.3,Ex);
	outGraphY->SetPoint(pt,pt*0.3,Ey);
	outGraphZ->SetPoint(pt,pt*0.3,Ez);
      }



      delete inTree;
    }
  }

  outFile->Close();


  return;
}
