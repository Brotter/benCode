/*
  Should read in a bunch of timefresnel-root.dat files and then maybe plot them or 
  find interesting things about them


  Ben Rotter - Feb 2017 - University of Hawaii at Manoa
 */


void ZHAiresReader(int zen=8) {


  string baseDir = "/storageb/harmscho/ANITA_CR_PAPER/";

  string subDir = "exposure_events/";

  stringstream fullDir;
  fullDir.str("");
  fullDir << baseDir << subDir << "RefCR_zen8__arm";

  const int waveSize = 2000;
  TH2D *hWavesX = new TH2D("hWavesX","hWavesX",waveSize,-0.5,waveSize-0.5, 16*9,-0.5,16*9 - 0.5);
  TH2D *hWavesY = new TH2D("hWavesY","hWavesY",waveSize,-0.5,waveSize-0.5, 16*9,-0.5,16*9 - 0.5);
  TH2D *hWavesZ = new TH2D("hWavesZ","hWavesZ",waveSize,-0.5,waveSize-0.5, 16*9,-0.5,16*9 - 0.5);

  TGraph2D *antPos = new TGraph2D();
  antPos->SetName("antPos");
  antPos->SetTitle("antPos");

  TGraph *outGraphX;
  TGraph *outGraphY;
  TGraph *outGraphZ;

  stringstream filename;
  for (int arm=0; arm<4; arm++) {
    for (int sub=0; sub<4; sub++) {
      cout << "arm" << arm << " sub" << sub << endl;
      filename.str("");
      filename << fullDir.str() << arm+1 << "_sub" << sub+1 << "/timefresnel-root.dat";
      TTree *inTree = new TTree("inTree","inTree");
      inTree->ReadFile(filename.str().c_str(),"showerNum:antNum:antX:antY:antZ:time:A:Ax:Ay:Az:E:Ex:Ey:Ez");
      Float_t Ex,Ey,Ez,antx,anty,antz,antNum;
      inTree->SetBranchAddress("Ex",&Ex);
      inTree->SetBranchAddress("Ey",&Ey);
      inTree->SetBranchAddress("Ez",&Ez);
      inTree->SetBranchAddress("antX",&antx);
      inTree->SetBranchAddress("antY",&anty);
      inTree->SetBranchAddress("antZ",&antz);
      inTree->SetBranchAddress("antNum",&antNum);
      int prevAntNum = -1;



      name.str("");
      name << "arm" << arm << "sub" << sub;
      outGraph->SetTitle(name.str().c_str())
      outGraph->SetName(name.str().c_str())

      for (int entry=0; entry<inTree->GetEntries(); entry++) {
	inTree->GetEntry(entry);

	if (prevAntNum != antNum) { 
	  cout << "antNum " << antNum << endl;
	  antPos->SetPoint(antPos->GetN(),antx,anty,antz);
	  prevAntNum = antNum; }

	int pt = entry-20000*(antNum-1);
	outGraph->SetPoint(entry,pt,Ex);
	hWavesX->Fill(pt,(arm*4+sub)*9+antNum-1,Ex);
	hWavesY->Fill(pt,(arm*4+sub)*9+antNum-1,Ey);
	hWavesZ->Fill(pt,(arm*4+sub)*9+antNum-1,Ez);
      }
      delete inTree;
    }
  }

  filename.str("");
  filename << "ZHAiresReader_zen" << zen << ".root";
  TFile *outFile = TFile::Open(filename.str().c_str(),"recreate");
  hWavesX->Write();
  hWavesY->Write();
  hWavesZ->Write();
  antPos->Write();
  outFile->Close();


  return;
}
