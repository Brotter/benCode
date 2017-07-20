#include "AnitaConventions.h"
#include "BedmapReader.h"


TProfile2D* bedMapPlot(){

  BedmapReader *bedmap = new BedmapReader(0);

  const double latLow = -82;
  const double latHigh = -73;
  const double lonLow = -140;
  const double lonHigh = -90;

  const int numLatBins = 101;
  double latPerBin = (latHigh-latLow)/numLatBins;
  const int numLonBins = 101;
  double lonPerBin = (lonHigh-lonLow)/numLonBins;

  cout <<  latPerBin << " " << lonPerBin << endl;

  TProfile2D *hist = new TProfile2D("bedmap","bedmap;lat,lon",numLatBins-1,latLow,latHigh,numLonBins-1,lonLow,lonHigh);

  for (int latBin=0; latBin<numLatBins; latBin++) {
    double lat = latLow+latBin*latPerBin;

    for (int lonBin=0; lonBin<numLatBins; lonBin++) {
      double lon = lonLow+lonBin*lonPerBin;

      double value = bedmap->Surface(lon,lat);
      hist->Fill(lat,lon,value);
    }
  }

  hist->GetZaxis()->SetRangeUser(6374e3,6380e3);
  hist->SetStats(0);

  return hist;
}
    
  


void makeMovie() {
  
  
  TChain *newCorrelatorTree = new TChain("newCorrelatorTree","newCorrelatorTree");

  stringstream name;
  for (int i=0; i<256; i++) {
    name.str("");
    name << "rootFiles/waisNewCorrelator_" << i << ".root";
    newCorrelatorTree->Add(name.str().c_str());
  }

  Double_t lat,lon;
  newCorrelatorTree->SetBranchAddress("lat",&lat);
  newCorrelatorTree->SetBranchAddress("lon",&lon);

  int numEntries = newCorrelatorTree->GetEntries();
  cout << "number of entries: " << numEntries << endl;
  
  int numFrames = 250;
  int eventsPerFrame = numEntries/numFrames;


  TGraph *wais = new TGraph();
  wais->SetPoint(0,AnitaLocations::LATITUDE_WAIS,AnitaLocations::LONGITUDE_WAIS);
  wais->SetMarkerStyle(kFullStar);
  //  wais->SetMarkerColor(kRed);

  TProfile2D *bedmap = bedMapPlot();

  TH2D *hist = new TH2D("WAIS pointing","WAIS pointing;lat;lon",100,-82,-73,100,-140,-90);
  hist->SetStats(0);

  for (int frame=1; frame<numFrames; frame++) {
    hist->Reset();
    for (int entry=(frame-1)*eventsPerFrame; entry<frame*eventsPerFrame; entry++) {
      newCorrelatorTree->GetEntry(entry);
      hist->Fill(lon,lat);
    }

    TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
    hist->GetZaxis()->SetRangeUser(0,200);
    bedmap->Draw("col");
    hist->Draw("colzSame");
    wais->Draw("pSame");


    //    aMap->DrawTGraph("pSame");
    name.str(""); 
    name << "images/newWais_" << frame << ".png";
    
    c1->SaveAs(name.str().c_str(),"q");
    delete c1;
  }

  return;

}

void waisNewCorrelatorPlotter(){

  //    makeMovie();


  TChain *newCorrelatorTree = new TChain("newCorrelatorTree","newCorrelatorTree");

  stringstream name;
  for (int i=0; i<256; i++) {
    name.str("");
    name << "rootFiles/waisNewCorrelator_" << i << ".root";
    newCorrelatorTree->Add(name.str().c_str());
  }


  cout << "number of entries: " << newCorrelatorTree->GetEntries() << endl;

  newCorrelatorTree->Draw("peakThetaDeg-waisTheta:peakPhiDeg-waisPhi","","colz");

  return;

}

