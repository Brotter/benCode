{

  // define the TChains where abby stores the data (most of these are actually TNtuples but whatever)
  TChain *ndata = new TChain("ndata","ndata");
  TChain *ndata2 = new TChain("ndata2","ndata2");
  TChain *ndata3 = new TChain("ndata3","ndata3");
  TChain *ndata4 = new TChain("ndata4","ndata4");
  TChain *tdataPointed = new TChain("tdataPointed","tdataPointed");
  TChain *tdataEvent = new TChain("tdataEvent","tdataEvent");
  TChain *tdataTracing = new TChain("tdataTracing","tdataTracing");
  int startRun,endRun;



  //link all the processed runs to those chains
  // I linked each output directory on the 4 different servers via nfs
  stringstream name;
  for (int i=150; i<440; i++) {
    name.str("");
    name << "/Volumes/NO NAME/ANITA3/rootOutputs/output" << i << "_1.root";
    ndata->Add(name.str().c_str());
    ndata2->Add(name.str().c_str());
    ndata3->Add(name.str().c_str());
    ndata4->Add(name.str().c_str());
    tdataPointed->Add(name.str().c_str());
    tdataEvent->Add(name.str().c_str());
    tdataTracing->Add(name.str().c_str());
  }

  cout << tdataPointed->GetEntries() << " passing events out of " << tdataEvent->GetEntries();
  cout << " (" << 100.*tdataPointed->GetEntries()/float(tdataEvent->GetEntries()) << "%)" << endl;
  

  //link all the variables to the branches in the chains
  //tdataTracing tree
  int eventNumber,eventTracedFlag,quietBaseFlagIndex,realTime;
  double peakThetaFinal,peakPhiFinal,sourceLon,sourceLat,sourceAlt,anitaLatitude,anitaLongitude,anitaAltitude;

  //just the event number
  tdataTracing->SetBranchAddress("eventNumber",&eventNumber);
  //peak interferometry bins in payload coordinates IN RADIANS
  tdataTracing->SetBranchAddress("peakThetaFinal",&peakThetaFinal);
  tdataTracing->SetBranchAddress("peakPhiFinal",&peakPhiFinal);
  //the max direction pointed back to the continent (0 if eventTracedFlag==0)
  tdataTracing->SetBranchAddress("sourceLon",&sourceLon);
  tdataTracing->SetBranchAddress("sourceLat",&sourceLat);
  tdataTracing->SetBranchAddress("sourceHeight",&sourceAlt);
  //whether the event traced back to the continent
  tdataTracing->SetBranchAddress("eventTracedFlag",&eventTracedFlag);
  //where ANITA was at event
  tdataTracing->SetBranchAddress("anitaLatitude",&anitaLatitude);
  tdataTracing->SetBranchAddress("anitaLongitude",&anitaLongitude);
  tdataTracing->SetBranchAddress("anitaAltitude",&anitaAltitude);
  //seems to always be -1 and doesn't do anyting?
  tdataTracing->SetBranchAddress("quietBaseFlagIndex",&quietBaseFlagIndex);
  //time of event
  tdataTracing->SetBranchAddress("realTime",&realTime);

  //ndata nTuple
  Float_t deltaTheta,deltaPhi,deltamcmTheta,deltamcmPhi,thetaMap,phiMap,mapSNR,peakVal,ratioFirstToSecondPeak,snrCoherent,snrPeakAnt,maxSignalPeak,distanceTD,peakHilbertCoherent;
  //Taylor Dome angular separation
  ndata->SetBranchAddress("deltaTheta",&deltaTheta);
  ndata->SetBranchAddress("deltaPhi",&deltaPhi);
  //McMurdo angular separation
  ndata->SetBranchAddress("deltamcmTheta",&deltamcmTheta);
  ndata->SetBranchAddress("deltamcmPhi",&deltamcmPhi);
  //peak interferometry bins in payload coordinates IN DEGREES
  ndata->SetBranchAddress("thetaMap",&thetaMap);
  ndata->SetBranchAddress("phiMap",&phiMap);
  //height of interferometry peak vs RMS of map
  ndata->SetBranchAddress("mapSNR",&mapSNR);
  //peak value of interferometry map
  ndata->SetBranchAddress("peakVal",&peakVal);
  //ratio of first to second peak in interferometry map
  ndata->SetBranchAddress("ratioFirstToSecondPeak",&ratioFirstToSecondPeak);
  //snr = peak to peak of waveform divided by the rms noise of second half of waveform (maybe not great, Peter did this)
  //snr of the coherently summed waveform (12 antennas, 4 phi sectors?)
  ndata->SetBranchAddress("snrCoherent",&snrCoherent);
  //snr of the peak antenna (whole payload)
  ndata->SetBranchAddress("snrPeakAnt",&snrPeakAnt);
  //maximum peak of signal (p2p/2)
  ndata->SetBranchAddress("maxSignalPeak",&maxSignalPeak);
  //distance to taylor dome
  ndata->SetBranchAddress("distanceTD",&distanceTD);
  //peak of the hilbert transform
  ndata->SetBranchAddress("peakHilbertCoherent",&peakHilbertCoherent);

  //ndata2 nTuple
  Float_t deltaTTD,snrPeakAfterFilter,didIFilter,triggerOrPhiMaskFlag,thetaTD,phiTD,thetaWilly,phiWilly,hwTriggerAngle,thisPhiMask,distanceMcM,secondTheta,secondPhi,strongCWFlag;
  //time difference between trigTime and time of flight of signal from taylor dome
  ndata2->SetBranchAddress("deltaTTD",&deltaTTD);
  //200->1200MHz band pass filter
  ndata2->SetBranchAddress("snrPeakAfterFilter",&snrPeakAfterFilter);
  //how many filtered freq bands (5 total)
  ndata2->SetBranchAddress("didIFilter",&didIFilter);
  ndata2->SetBranchAddress("triggerOrPhiMaskFlag",&triggerOrPhiMaskFlag);
  ndata2->SetBranchAddress("thetaTD",&thetaTD);
  ndata2->SetBranchAddress("phiTD",&phiTD);
  ndata2->SetBranchAddress("thetaWilly",&thetaWilly);
  ndata2->SetBranchAddress("phiWilly",&phiWilly);
  ndata2->SetBranchAddress("hwTriggerAngle",&hwTriggerAngle);
  ndata2->SetBranchAddress("thisPhiMask",&thisPhiMask);
  ndata2->SetBranchAddress("distanceMcM",&distanceMcM);
  ndata2->SetBranchAddress("secondTheta",&secondTheta);
  ndata2->SetBranchAddress("secondPhi",&secondPhi);
  ndata2->SetBranchAddress("strongCWFlag",&strongCWFlag);

  Float_t headingOfThisEvent,nadirFlag,thirdTheta,thirdPhi,varnerFlag,varnerFlag2,pitch,roll,heading,phiMaskFlag,hwTriggerFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,xCorPassFlag;

  ndata3->SetBranchAddress("headingOfThisEvent",&headingOfThisEvent);
  ndata3->SetBranchAddress("nadirFlag",&nadirFlag);
  ndata3->SetBranchAddress("thirdTheta",&thirdTheta);
  ndata3->SetBranchAddress("thirdPhi",&thirdPhi);
  ndata3->SetBranchAddress("varnerFlag",&varnerFlag);
  ndata3->SetBranchAddress("varnerFlag2",&varnerFlag2);
  ndata3->SetBranchAddress("pitch",&pitch);
  ndata3->SetBranchAddress("roll",&roll);
  ndata3->SetBranchAddress("heading",&heading);
  ndata3->SetBranchAddress("phiMaskFlag",&phiMaskFlag);
  ndata3->SetBranchAddress("hwTriggerFlag",&hwTriggerFlag);
  ndata3->SetBranchAddress("ratioOfPeaksPassFlag",&ratioOfPeaksPassFlag);
  ndata3->SetBranchAddress("elevationAnglePassFlag",&elevationAnglePassFlag);
  ndata3->SetBranchAddress("xCorPassFlag",&xCorPassFlag);

  Float_t payloadBlastFlag,polAngleCoherent,polFractionCoherent,didIFilterAboveSatellite,didIFilterHoriz,didIFilterAboveSatelliteHoriz,meanFreqVert,meanFreqHoriz,deltaWAISTheta,deltaWAISPhi,distanceWAIS,uhecrCorrMax,uhecrCorrMin,trigTns;
  ndata4->SetBranchAddress("payloadBlastFlag",&payloadBlastFlag);
  ndata4->SetBranchAddress("polAngleCoherent",&polAngleCoherent);
  ndata4->SetBranchAddress("polFractionCoherent",&polFractionCoherent);
  ndata4->SetBranchAddress("didIFilterAboveSatellite",&didIFilterAboveSatellite);
  ndata4->SetBranchAddress("didIFilterHoriz",&didIFilterHoriz);
  ndata4->SetBranchAddress("didIFilterAboveSatelliteHoriz",&didIFilterAboveSatelliteHoriz);
  ndata4->SetBranchAddress("meanFreqVert",&meanFreqVert);
  ndata4->SetBranchAddress("meanFreqHoriz",&meanFreqHoriz);
  ndata4->SetBranchAddress("deltaWAISTheta",&deltaWAISTheta);
  ndata4->SetBranchAddress("deltaWAISPhi",&deltaWAISPhi);
  ndata4->SetBranchAddress("distanceWAIS",&distanceWAIS);
  ndata4->SetBranchAddress("uhecrCorrMax",&uhecrCorrMax);
  ndata4->SetBranchAddress("uhecrCorrMin",&uhecrCorrMin);
  ndata4->SetBranchAddress("trigTns",&trigTns);



  int numEntries = tdataTracing->GetEntries();

  int startTime,endTime;
  tdataTracing->GetEntry(0);
  startTime = realTime;
  tdataTracing->GetEntry(numEntries-1);
  endTime = realTime;


  //  tdataPointed->Draw("sourceLat:sourceLon:eventNumber","","colz");
  gStyle->SetOptStat(0);
  tdataTracing->SetMarkerStyle(0);
  tdataPointed->SetMarkerStyle(0);
  tdataPointed->SetMarkerColor(kRed);

  AntarcticaMapPlotter *aMap = new AntarcticaMapPlotter("mapPlotter","mapPlotter",1000,1000);
  aMap->addHistogram("pointMap","pointMap",1000,1000);
  aMap->addProfile("PeakHilbert","PeakHilbert",1000,1000);
  TProfile2D* pPeakHilbert = (TProfile2D*)aMap->getCurrentHistogram();


  Double_t sunPhi,sunTheta;  
    
  Double_t sunDiff;


  for (int entry=0; entry<numEntries; entry++) {
  //  for (int entry=1000000; entry<1200000; entry++) {
    if (entry%10000==0) cout << entry << "/" << numEntries << endl;
    tdataTracing->GetEntry(entry);
    ndata->GetEntry(entry);
    ndata2->GetEntry(entry);
    ndata3->GetEntry(entry);
    ndata4->GetEntry(entry);
    aMap->setCurrentHistogram("pointMap");
    aMap->Fill(sourceLat,sourceLon);
    aMap->setCurrentHistogram("PeakHilbert");
    aMap->Fill(sourceLat,sourceLon,peakHilbertCoherent);
    
  }

  TCanvas *c4 = new TCanvas();
  aMap->DrawHist("colz");

  TFile *outFile = TFile::Open("antarcticaMap.root","recreate");
  aMap->Write();
  outFile->Close();

    
  return;
}
