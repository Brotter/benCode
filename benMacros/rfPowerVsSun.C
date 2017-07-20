{


  cout << "Starting physics!" << endl;

  AnitaGeomTool *geom = new AnitaGeomTool();

  stringstream name;

  //to make things run faster
  int startRun = 150;
  int stopRun = 440;
  int maxEntries = -1;

  //link in the GPS trees
  TChain *patTree = new TChain("adu5PatTree","adu5PatTree");
  patTree->Add("/Volumes/ANITA3Data/gpsFileAll.root");

  Adu5Pat *patPtr = NULL;
  patTree->SetBranchAddress("pat",&patPtr);
  patTree->BuildIndex("realTime");

  TChain *surfHkTree = new TChain("surfHkTree","surfHkTree");
  for (int i=startRun; i<stopRun; i++) {
    name.str("");
    name << "/Volumes/ANITA3Data/root/run" << i << "/surfHkFile" << i << ".root";
    surfHkTree->Add(name.str().c_str());
  }
  SurfHk *surfPtr = NULL;
  surfHkTree->SetBranchAddress("surf",&surfPtr);
  surfHkTree->BuildIndex("realTime");
  
  int numEntries = 0;
  if (maxEntries == -1)  numEntries = surfHkTree->GetEntries();
  else numEntries = maxEntries;
  cout << "Number of entries: " << numEntries << endl;


  surfHkTree->GetEntry(0);
  int firstTime = surfPtr->realTime;
  surfHkTree->GetEntry(numEntries-1);
  int lastTime = surfPtr->realTime;

  double sunPhi,sunTheta;

  int surf,chan,ant;
  int ring,pol;

  TProfile2D *rfPowVsSunVsTime = new TProfile2D("rfPowVsSunVsTime","rfPowVsSunVsTime",
						numEntries/1000,firstTime,lastTime,
						16,0,360);

  TProfile2D *rfPowVsSunVsTimeV = new TProfile2D("rfPowVsSunVsTimeV","rfPowVsSunVsTimeV",
						numEntries/1000,firstTime,lastTime,
						16,0,360);

  TProfile2D *rfPowVsSunVsTimeH = new TProfile2D("rfPowVsSunVsTimeH","rfPowVsSunVsTimeH",
						numEntries/1000,firstTime,lastTime,
						16,0,360);


  TProfile2D *rfPowVsChan = new TProfile2D("rfPowVsChan","rfPowVsChan",
					   numEntries/1000,firstTime,lastTime,
					   96,-0.5,95.5);

  for (int entry=0; entry<numEntries; entry++) {
    if (entry%1000 == 0) cout << entry << " / " << numEntries << "\r";

    surfHkTree->GetEntry(entry);
    int gpsEntry = patTree->GetEntryNumberWithBestIndex(surfPtr->realTime);
    patTree->GetEntry(gpsEntry);

//    UsefulAdu5Pat *useful = new UsefulAdu5Pat(patPtr);
//    useful->getSunPosition(sunPhi,sunTheta);
//    delete useful;

    for (int phi=0; phi<16; phi++) {
      double phiDiff = phi*22.5 - patPtr->heading;
      if (phiDiff > 360) phiDiff -= 360;
      else if (phiDiff < 0 ) phiDiff += 360;
      for (ring=0; (AnitaRing::AnitaRing_t)ring != AnitaRing::kNotARing; ring++) {
	for (pol=0; (AnitaPol::AnitaPol_t)pol != AnitaPol::kNotAPol; pol++) {
	  int index = phi*6+ring*2+pol;
	  AnitaRing::AnitaRing_t ringRing = (AnitaRing::AnitaRing_t)ring;
	  AnitaPol::AnitaPol_t polPol = (AnitaPol::AnitaPol_t)pol;
	  geom->getSurfChanAntFromRingPhiPol(ringRing,phi,polPol,surf,chan,ant);

	  rfPowVsChan->Fill(surfPtr->realTime,index,surfPtr->getRFPowerInK(surf,chan));

	  if ( (index == 95) || (index == 94) ||
	       (index == 47) || (index == 43) ||
	       (index == 24) || (index == 22)) {
	    if (entry==0) {
	      cout << "skipping: index=" << index << " AKA " << phi+1  << AnitaRing::ringAsChar((AnitaRing::AnitaRing_t)ring) << AnitaPol::polAsChar((AnitaPol::AnitaPol_t)pol) << endl;
	    }
	    continue;
	  }
	  rfPowVsSunVsTime->Fill(surfPtr->realTime,phiDiff,surfPtr->getRFPowerInK(surf,chan));

	  if (polPol == AnitaPol::kVertical) {
	    rfPowVsSunVsTimeV->Fill(surfPtr->realTime,phiDiff,surfPtr->getRFPowerInK(surf,chan)); }
	  else if (polPol == AnitaPol::kHorizontal) {
	    rfPowVsSunVsTimeH->Fill(surfPtr->realTime,phiDiff,surfPtr->getRFPowerInK(surf,chan)); }
	}
      }
    
    }
  }

  TCanvas *c1 = new TCanvas();
  rfPowVsSunVsTime->Draw("colz");
  TCanvas *c2 = new TCanvas();
  rfPowVsChan->Draw("colz");

  TFile *outFile = TFile::Open("rfPowerVsSun.root","recreate");
  rfPowVsSunVsTime->Write();
  rfPowVsSunVsTimeV->Write();
  rfPowVsSunVsTimeH->Write();
  rfPowVsChan->Write();
  outFile->Close();


}
