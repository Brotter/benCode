void sunPointing(){

  cout << "Starting physics!" << endl;

  stringstream name;

  //link in the GPS trees too
  TChain *patTree = new TChain("adu5PatTree","adu5PatTree");
  for (int i=150; i<439; i++) {
    name.str("");
    name << "/home/brotter/anita16/ANITA3_rootData/root/run" << i << "/gpsFile" << i << ".root";
    patTree->Add(name.str().c_str());
  }
  Adu5Pat *patPtr = NULL;
  patTree->SetBranchAddress("pat",&patPtr);
  patTree->BuildIndex("realTime");


  TGraph *gSunPhi = new TGraph();
  gSunPhi->SetTitle("Phi location of Sun");
  gSunPhi->SetName("gSunPhi");
  TGraph *gSunTheta = new TGraph();
  gSunTheta->SetTitle("Theta location of Sun");
  gSunTheta->SetName("gSunTheta");


  Double_t sunPhi,sunTheta;

  int numEntries = patTree->GetEntries();

  for (int entry=0; entry<numEntries; entry++) {
    if (entry%1000!=0) continue;

    cout << entry << "/" << numEntries << endl;
    
    patTree->GetEntry(entry);
    UsefulAdu5Pat *useful = new UsefulAdu5Pat(patPtr);
    useful->getSunPosition(sunPhi,sunTheta);
    delete useful;

    sunPhi = patPtr->heading - sunPhi;
    while(sunPhi < 0) sunPhi += 360;
    while(sunPhi >= 0) sunPhi -= 360;


    gSunPhi->SetPoint(gSunPhi->GetN(),entry,sunPhi);
    gSunTheta->SetPoint(gSunTheta->GetN(),entry,sunTheta);
    

  }//end entry loop


  TCanvas *c1 = new TCanvas();
  gSunPhi->Draw("alp");
  TCanvas *c2 = new TCanvas();
  gSunTheta->Draw("alp");
    
    
  return;
}

