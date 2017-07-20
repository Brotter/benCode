{

  cout << "Starting Physics! :) " << endl;
  // define the TChains where abby stores the data (most of these are actually TNtuples but whatever)
  TChain *ndata = new TChain("ndata","ndata");
  TChain *tdataTracing = new TChain("tdataTracing","tdataTracing");



  //link all the processed runs to those chains
  // I linked each output directory on the 4 different servers via nfs
  stringstream name;
  for (int i=150; i<440; i++) {
    //these runs are missing or didn't finish MyCorrelator (retransfer later)
    if ((i>=257 && i<=263) || (i==170) || (i==171) || (i==387) || (i==388) || (i==395)) continue;
    name.str("");
    name << "/Volumes/ANITA3Data/rootOutputs/output" << i << "_0.root";
    ndata->Add(name.str().c_str());
    tdataTracing->Add(name.str().c_str());
  }

  

  //link all the variables to the branches in the chains
  //tdataTracing tree
  Int_t eventNumber,realTime;
  Double_t peakThetaFinal,peakPhiFinal,sourceLon,sourceLat,sourceHeight,anitaLatitude,anitaLongitude,anitaAltitude;

  //just the event number
  tdataTracing->SetBranchAddress("eventNumber",&eventNumber);
  //peak interferometry bins in payload coordinates IN RADIANS
  tdataTracing->SetBranchAddress("peakThetaFinal",&peakThetaFinal);
  tdataTracing->SetBranchAddress("peakPhiFinal",&peakPhiFinal);
  //the max direction pointed back to the continent (0 if eventTracedFlag==0)
  tdataTracing->SetBranchAddress("sourceLon",&sourceLon);
  tdataTracing->SetBranchAddress("sourceLat",&sourceLat);
  tdataTracing->SetBranchAddress("sourceHeight",&sourceHeight);
  //where ANITA was at event
  tdataTracing->SetBranchAddress("anitaLatitude",&anitaLatitude);
  tdataTracing->SetBranchAddress("anitaLongitude",&anitaLongitude);
  tdataTracing->SetBranchAddress("anitaAltitude",&anitaAltitude);
  //time of event
  tdataTracing->SetBranchAddress("realTime",&realTime);

  //ndata nTuple
  Float_t thetaMap,phiMap;
  ndata->SetBranchAddress("thetaMap",&thetaMap);
  ndata->SetBranchAddress("phiMap",&phiMap);
 

  //link in the GPS trees too
  TChain *patTree = new TChain("adu5PatTree","adu5PatTree");
  patTree->Add("/Volumes/ANITA3Data/gpsFileAll.root");
  Adu5Pat *patPtr = NULL;
  patTree->SetBranchAddress("pat",&patPtr);
  patTree->BuildIndex("realTime");


  int numEntries = tdataTracing->GetEntries();

  int startTime,endTime;
  tdataTracing->GetEntry(0);
  startTime = realTime;
  tdataTracing->GetEntry(numEntries-1);
  endTime = realTime;


  gStyle->SetOptStat(0);
  tdataTracing->SetMarkerStyle(0);


  TH2D *thing1 = new TH2D("myHist","Pointing offset from sun;phi;theta",
			  1000,-360,360,
			  1000,-100,100);

  TH1D *gpsOffset = new TH1D("gpsOffset","Diff between GPS time and event time;dT;coutns",
			     1000,-.1,.1);

  TH2D *sunPointingPhi = new TH2D("sunPointingPhi","Phi Separation between event and sun",
			      1000,startTime,endTime,
			      1000,-360,360);

  TH2D *sunDiffTheta = new TH2D("sunDiffTheta","Theta Separation between event and sun",
			      1000,startTime,endTime,
			      1000,-5,5);



  TH2D *sunReflectDiffTheta = new TH2D("sunReflectDiffTheta","Theta Separation between event and sun reflection",
			      1000,startTime,endTime,
			      1000,-5,5);



  TGraph *gSunTheta = new TGraph();
  gSunTheta->SetName("gSunTheta");


  TGraph *gSunReflect = new TGraph();
  gSunReflect->SetName("gSunReflect");

  double sunPhi,sunTheta;  
    
  Double_t sunDiff;

  Double_t sunReflect = 0.0;
  Double_t temp = 0.0;
  Double_t pi = 3.14159;
  Double_t Rearth = 6.371e6; //in m
  Double_t degToRad = 3.14159/180.;

  for (int entry=0; entry<numEntries; entry++) {
  //  for (int entry=1000000; entry<3500000; entry++) {
    if (entry%10000==0) cout << entry << "/" << numEntries << "\r";
    tdataTracing->GetEntry(entry);
    ndata->GetEntry(entry);
    
    patTree->GetEntry(patTree->GetEntryNumberWithBestIndex(realTime));
    gpsOffset->Fill(patPtr->realTime-realTime);
    UsefulAdu5Pat *useful = new UsefulAdu5Pat(patPtr);
    useful->getSunPosition(sunPhi,sunTheta);
    delete useful;

    sunTheta *= -1; //This is inverted from how Abby reports theta in thetaMap(negative down)

    //    temp = 0.5*(TMath::Sin(sunTheta*degToRad)*TMath::Tan(sunTheta*degToRad)) - TMath::Cos(sunTheta*degToRad);
    //    sunReflect = TMath::ATan((1/3.)*(-1*TMath::Sqrt(pow(temp,2)+3)-temp));

    //this seems to work, though I don't know why (I spent a WEEK on this trig D: )
//    temp = TMath::Cos(sunTheta*degToRad) + TMath::Sin(sunTheta*degToRad)/2.;
//    if (sunTheta < 26.6) {
//      sunReflect = 2*TMath::ATan((1-TMath::Sqrt(5-4*pow(temp,2)))/(2*(temp+1)));
//    }
//    else {
//      sunReflect = 2*TMath::ATan((1+TMath::Sqrt(5-4*pow(temp,2)))/(2*(temp+1)));
//    }
    

    //Solved by Brian Davis (Joan's mom's friend!)
    sunReflect = TMath::ASin(Rearth/(Rearth+patPtr->altitude)*TMath::Cos((-1*sunTheta*degToRad)/2));

    sunReflect /= degToRad;
    if (entry%1000==0) {
      gSunTheta->SetPoint(gSunTheta->GetN(),realTime,sunTheta);
      gSunReflect->SetPoint(gSunReflect->GetN(),realTime,-sunReflect);
    }
    
    sunDiff = phiMap-sunPhi -360;
    
    //    cout << heading << " " << phiMap << " " << sunAz << " | ";
    //    cout << sunTheta << " " << thetaMap << endl;

    thing1->Fill(sunDiff,thetaMap-sunTheta);
    sunPointingPhi->Fill(realTime,sunDiff);
    sunDiffTheta->Fill(realTime,thetaMap-sunTheta);
    sunReflectDiffTheta->Fill(realTime,thetaMap-sunReflect);
  }

  //  TCanvas *c1 = new TCanvas("c","c",4000,2000);


  //TCanvas *c1 = new TCanvas();
  //thing1->Draw("colz");
  //TCanvas *c2 = new TCanvas();
  //sunPointingPhi->Draw("colz");
  //TCanvas *c3 = new TCanvas();
  //sunDiffTheta->Draw("colz");


  TFile *outFile = TFile::Open("sunVsPointing.root","recreate");
  thing1->Write();
  gpsOffset->Write();
  sunPointingPhi->Write();
  sunDiffTheta->Write();
  sunReflectDiffTheta->Write();
  gSunTheta->Write();
  gSunReflect->Write();
  
  outFile->Close();

  //  c1->SaveAs("cutEventPointing.tiff");
    
  return;
}
