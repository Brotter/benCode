{

  TTree *inTree = new TTree("inTree","inTree");
  inTree->ReadFile("/Users/brotter/anita16/local/share/anitaCalib/relativeCableDelays.dat",
		   "surf/I:chan/I:lab/I:delay/D");
  int surf,chan,lab;
  double delay;
  inTree->SetBranchAddress("surf",&surf);
  inTree->SetBranchAddress("chan",&chan);
  inTree->SetBranchAddress("lab",&lab);
  inTree->SetBranchAddress("delay",&delay);


  TProfile2D *relativeCableDelays = new TProfile2D("hDelays","Channel average invarient delay with respect to 16BH;Surf;Channel;Relative Channel Delay (ns)",
						   12,0,12, 8,0,8);

  for (int entry=0; entry<inTree->GetEntries(); entry++) {
    inTree->GetEntry(entry);
    //    if ( lab != 0) continue;

    relativeCableDelays->Fill(surf+0.5,chan+0.5,-delay);

  }

  relativeCableDelays->SetMarkerSize(1);
  gStyle->SetPaintTextFormat("4.2f");
  relativeCableDelays->GetYaxis()->SetTitleOffset(0.5);
  relativeCableDelays->SetStats(0);
  relativeCableDelays->Draw("text colz");

}
