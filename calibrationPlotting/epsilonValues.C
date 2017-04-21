{


  TTree *inTree = new TTree("inTree","inTree");

  inTree->ReadFile("/Users/brotter/anita16/local/share/anitaCalib/epsilonFromBenS.dat",
		   "surf/I:lab/I:rco/I:epsilon/D");

  int surf,lab,rco;
  double epsilon;

  inTree->SetBranchAddress("surf",&surf);
  inTree->SetBranchAddress("lab",&lab);
  inTree->SetBranchAddress("rco",&rco);
  inTree->SetBranchAddress("epsilon",&epsilon);

 
  TH2D *hist = new TH2D("hist","Eplison Values;Surf;Lab*2 + RCO; Epsilon (ns)",12,0,12,8,0,8);

  for (int entry=0; entry<inTree->GetEntries(); entry++) {
    inTree->GetEntry(entry);
    hist->Fill(surf,lab*2 + rco,epsilon);
  }

  gStyle->SetPaintTextFormat("4.2f");
  hist->SetMarkerSize(1);
  hist->GetYaxis()->SetTitleOffset(0.5);
  hist->SetStats(0);
  hist->Draw("text colz");

}
