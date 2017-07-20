{
  cout << "Starting Physics! :) " << endl;


  //import all the surfHk data
  cout << "Loading TChains: ";
  TChain *surfHkTree = new TChain("surfHkTree","surfHkTree");
  stringstream name;
  for (int i=150; i<440; i++) {
    name.str("");
    name << "/Volumes/ANITA3Data/surfHk/surfHkFile" << i << ".root";
    surfHkTree->Add(name.str().c_str());
  }
  SurfHk *surf = NULL;
  surfHkTree->SetBranchAddress("surf",&surf);
  cout << "surfHk, ";

  //link in the GPS trees too
  TChain *patTree = new TChain("adu5PatTree","adu5PatTree");
  patTree->Add("/Volumes/ANITA3Data/gpsFileAll.root");
  Adu5Pat *patPtr = NULL;
  patTree->SetBranchAddress("pat",&patPtr);
  cout << "gps." << endl;
  patTree->BuildIndex("realTime");


  //how many entries?
  int numEntries = surfHkTree->GetEntries();
  cout << "TChains loaded, " << numEntries << " entries" << endl;


  //times
  surfHkTree->GetEntry(0);
  int startTime = surf->realTime;
  surfHkTree->GetEntry(numEntries-1);
  int endTime = surf->realTime;

  //storage objects
  TProfile2D *thresholdsV = new TProfile2D("thresholdsV","thresholdsV",1000,startTime,endTime,32,-360,360);
  TProfile2D *thresholdsH = new TProfile2D("thresholdsH","thresholdsH",1000,startTime,endTime,16,-360,360);


  //since the diodes aren't identical, we should find the MEAN threshold for each chan
  cout << "Calculating means..." << endl;
  double meansV[16] = {0};
  double meansH[16] = {0};

  //this is the same calculation every time, so I saved the output as thresholdMeans.txt
  FILE *meanFile;
  meanFile = fopen("thresholdMeans.txt", "r");
  //however, if it doesn't exist, remake it
  if(meanFile == NULL) {
    cout << "Can't find mean file!" << endl;
    for (int entry=0; entry<numEntries; entry++) {
      if (entry%10000 == 0) {
	cout << entry << " / " << numEntries << "\r";
	fflush(stdout);
      }
      surfHkTree->GetEntry(entry);
      for (int phi=0; phi<16; phi++) {
	for (int ring_i=0; (AnitaRing::AnitaRing_t)ring_i != AnitaRing::kNotARing; ring_i++) {
	  AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t) ring_i;
	  meansV[phi] += surf->getThreshold(phi,ring,AnitaPol::kVertical); //remember to average later!
	  meansH[phi] += surf->getThreshold(phi,ring,AnitaPol::kHorizontal);
	}
      }
    } 
    cout << "Recreating mean file and writing" << endl;
    meanFile = fopen("thresholdMeans.txt","w");
    for (int phi=0; phi<16; phi++) {
      meansH[phi] /= numEntries*3; //remember that you've just been adding them!
      meansV[phi] /= numEntries*3;
      fprintf(meanFile,"%lf %lf\n",meansV[phi],meansH[phi]);
    }
  }
  //if the threshold file DOES exist, just import it
  else {
    cout << "Using existing mean file" << endl;
    for (int phi=0; phi<16; phi++) {
      fscanf(meanFile,"%lf %lf",&meansV[phi],&meansH[phi]);
    }
  }

  fclose(meanFile);
  cout << "means:" << endl;
  for (int phi=0; phi<16; phi++) {
    cout << "phi " << phi << ": H=" << meansH[phi] << " V=" << meansV[phi] << endl;
  }




  //  for (int entry=0; entry<numEntries; entry++) {
  for (int entry=0; entry<numEntries; entry++) {
    if (entry%10000 == 0) {
      cout << entry << " / " << numEntries << "\r";
      fflush(stdout);
    }
    surfHkTree->GetEntry(entry);

    patTree->GetEntry(patTree->GetEntryNumberWithBestIndex(surf->realTime));
    
    UsefulAdu5Pat *useful = new UsefulAdu5Pat(patPtr);


    for (int phi=0; phi<16; phi++) {
      double phiAngle = phi*22.5;
      double sunAngle = useful->getDifferencePointingToSun(phiAngle);
      for (int ring_i=0; (AnitaRing::AnitaRing_t)ring_i != AnitaRing::kNotARing; ring_i++) {
	AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t) ring_i;
	thresholdsV->Fill(surf->realTime,sunAngle,
			  surf->getThreshold(phi,ring,AnitaPol::kVertical)/meansV[phi]);
	thresholdsH->Fill(surf->realTime,sunAngle,
			  surf->getThreshold(phi,ring,AnitaPol::kHorizontal)/meansH[phi]);
      }
    }
    delete useful;
  }




  return;
}
  
