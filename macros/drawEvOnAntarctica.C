#include "AntarcticaMapPlotter.h"
#include "AnitaEventSummary.h"
#include "AnitaConventions.h"
#include <Riostream.h>

/*

  I want to plot the events onto the continent, probably with the base list, so that I can see how likely they are an event


 */



int evToRun(int ev) {

  ifstream runToEv("/Users/brotter/anita16/benMacros/runToEv.txt");
  int run,evLow,evHigh;
  int runOut = -1;
  while (runToEv >> run >> evLow >> evHigh)
    if (ev >= evLow && ev <= evHigh) {
      runOut = run;
    }

  runToEv.close();
  
  return runOut;

}


  
TGraph *loadBases(Acclaim::AntarcticaMapPlotter *aMap) {


  TGraph *gBases = new TGraph();

  stringstream name;
  char* anitaInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  name.str("");
  name << anitaInstallDir << "/share/anitaCalib/baseListA3.root";
  TFile *fBaseList = TFile::Open(name.str().c_str());
  if (fBaseList == NULL) {
    cout << "Couldn't find base list in " << name.str().c_str() << ". Not drawing bases." << endl;
    return gBases;
  }

  TTree *baseCampTree = (TTree*)fBaseList->Get("baseCampTree");
  TTree *awsTree = (TTree*)fBaseList->Get("awsTree");
  TTree *fixedWingTree = (TTree*)fBaseList->Get("fixedWingTree");
  TList *treeList = new TList();
  treeList->Add(baseCampTree);
  treeList->Add(awsTree);
  treeList->Add(fixedWingTree);
  TTree *allBases = TTree::MergeTrees(treeList);
  double lat,lon;
  allBases->SetBranchAddress("fullLong",&lon);
  allBases->SetBranchAddress("fullLat",&lat);
  
  for (int entry=0; entry<allBases->GetEntries(); entry++) {
    allBases->GetEntry(entry);
    double x,y;
    aMap->getRelXYFromLatLong(lat,lon,x,y);    
    gBases->SetPoint(entry,x,y);
  }

  delete treeList;
  fBaseList->Close();

  return gBases;

}



void drawOnAntarcticaFromLatLonList() {

  ifstream inFile("passingEvs.txt");

  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();
  TH2D *evMap = aMap->addHistogram("evMap","evMap",250,250);

  TGraph *gBases = loadBases(aMap);
  gBases->SetMarkerColor(kRed);
  gBases->SetMarkerStyle(3);
  
  int row,instance,evNum;
  double lat,lon;
  while (inFile >> row >> instance >> evNum >> lat >> lon) {
    if (lat != -9999 && lon != -9999) continue;
    double x,y;
    aMap->getRelXYFromLatLong(lat,lon,x,y);
    evMap->Fill(x,y);
  }
  
  aMap->DrawHist("colz");
  aMap->DrawTGraph("pSame");


}



TProfile2D* makeCutHist() {  

  TChain *summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");

  TProfile2D *cutLatLon = new TProfile2D("cutLatLon","cutLatLon",250,-90,-65,360,-180,180);

  summaryTree->Draw("latitude:longitude:templateCRayH[5] >> bothCuts","flags.pulser == 0","colz");


  return cutLatLon;
}


void drawOnAntarcticaFromLatLonFile(string inFileName="passingLocations") {

  TFile* inFile = TFile::Open("passingLatLons.root");

  TH2D* latLons = (TH2D*)inFile->Get("passingLatLons");


  return;
}



TH2D* drawOnAntarcticaFromLatLonHist(TH2* latLons,Acclaim::AntarcticaMapPlotter *aMap) {


  TH2D *evMap = aMap->addHistogram("evMap","evMap",250,250);

  for (int latBin=0; latBin<latLons->GetNbinsX(); latBin++) {
    for (int lonBin=0; lonBin<latLons->GetNbinsY(); lonBin++) {
      cout << latBin << " " << lonBin << endl;
      int binValue = latLons->GetBinContent(latBin,lonBin);
      double lonValue = latLons->GetYaxis()->GetBinCenter(lonBin);
      double latValue = latLons->GetXaxis()->GetBinCenter(latBin);
      double x,y;
      aMap->getRelXYFromLatLong(latValue,lonValue,x,y);
      evMap->Fill(x,y,binValue);
    }
  }

  return evMap;
}



void drawOnAntarcticaFromCuts() {

  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();

  TGraph *gBases = loadBases(aMap);
  gBases->SetMarkerColor(kRed);
  gBases->SetMarkerStyle(3);
  
  TCanvas *c1 = new TCanvas("c1","c1",1024,800);

  aMap->img->Draw();
  gBases->Draw("pSame");

  cout << "loaded bases" << endl;

  TProfile2D *cutHist = makeCutHist();
  cout << "made cut histogram" << endl;
  TH2D* antCutHist = drawOnAntarcticaFromLatLonHist(cutHist,aMap);
  cout << "projected it onto antarctica" << endl;

  antCutHist->Draw("colzSame");

  c1->SaveAs("antTemplateMap.png");

  return;
}  




    


void drawHistOnAntarcticaFromCuts_old(){
  
  TChain *summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");

  const int latBins = 250;
  const int lonBins = 360;
  TH2D *cutLatLon = new TH2D("cutLatLon","cutLatLon",latBins,-90,-65,lonBins,-180,180);

  stringstream cuts;

  //evs not associated with a cal pulser
  cuts << "flags.pulser == 0";
  //polarization angle within 15 degs of horizontal
  cuts << " && ";
  cuts << "TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 30";
  //linear polarization fraction
  cuts << " && ";
  cuts << "TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4";
  //template
  cuts << " && ";
  cuts << "templateCRayH[5] > 0.7";
  //map peak
  cuts << " && ";
  cuts << "peak[0][0].value > 0.06";
  //not close to ldb
  cuts << " && ";
  cuts << "TMath::Abs(peak[0][0].phi - ldb.phi) > 5";
  //not close to wias
  cuts << " && ";
  cuts << "TMath::Abs(peak[0][0].phi - ldb.phi) > 5";
  //hits the continent
  cuts << " && ";
  cuts << "peak[0][0].latitude > -999";
  //not close to some super loud base
  //  cuts << " && ";
  //  cuts << "!( (peak[0][0].latitude < -78 && peak[0][0].latitude > -82) && (peak[0][0].longitude < -105 && peak[0][0].longitude > -115) )";
    

  summaryTree->Draw("peak[0][0].longitude:peak[0][0].latitude >> cutLatLon",cuts.str().c_str(),"colz");

  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();
  TH2D *evMap = aMap->addHistogram("evMap","evMap",250,250);

  TGraph *gBases = loadBases(aMap);
  gBases->SetMarkerColor(kRed);
  gBases->SetMarkerStyle(3);


  for (int latBin=0; latBin<latBins; latBin++) {
    for (int lonBin=0; lonBin<lonBins; lonBin++) {
      int binValue = cutLatLon->GetBinContent(latBin,lonBin);
      double lonValue = cutLatLon->GetYaxis()->GetBinCenter(lonBin);
      double latValue = cutLatLon->GetXaxis()->GetBinCenter(latBin);
      double x,y;
      aMap->getRelXYFromLatLong(latValue,lonValue,x,y);
      evMap->Fill(x,y,binValue);
    }
  }

  aMap->DrawHist("colz");
  gBases->Draw("psame");

}
  

  




void drawOnAntarctica_slow(string fileName="") {

  /*

    This is a stupidly slow way of drawing the events.

    TChain::Draw() is way faster, then I can just parse the resulting TH2D

   */

  char* resultsDir = getenv("ANITA3_RESULTSDIR");
  string date="06.11.17_19h/";

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  AnitaEventSummary *summary = NULL;

  stringstream name;

  TGraph *passingEvs;

  if (fileName != "") {
    for (int run=130;run<440;run++) {
      name.str("");
      name << resultsDir << date << run << ".root";
      summaryTree->Add(name.str().c_str());
    }
    cout << "There are " << summaryTree->GetEntries() << " entries in the summary tree" << endl;

    summaryTree->SetBranchAddress("eventSummary",&summary);

    cout << "Building index..."; 
    fflush(stdout);
    summaryTree->BuildIndex("eventNumber");
    cout << " done!" << endl;



    if (fileName != "") {
      passingEvs = new TGraph(fileName.c_str(),"%lg %lg*");
      cout << "found " << passingEvs->GetN() << " passing events in " << fileName << endl;
    }
    
  }
  else {
    passingEvs = new TGraph();
    cout << "not doing any events" << endl;
  }
  


  
  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();
  aMap->addTGraph("event","event");
  TGraph *gEv = aMap->getCurrentTGraph();
  gEv->SetMarkerStyle(29); //29=star
  gEv->SetMarkerColor(kGreen);
  
  TH2D *myHist = aMap->addHistogram("hist","hist",1000,1000);
  

  for (int ev=0; ev<passingEvs->GetN(); ev++) {
    int summaryEntry = summaryTree->GetEntryNumberWithBestIndex(passingEvs->GetY()[ev]);
    cout << ev << " summaryEntry=" << summaryEntry << endl;
    summaryTree->GetEntry(summaryEntry);

    double xEv,yEv;
    double latEv = summary->peak[0][0].latitude;
    double lonEv = summary->peak[0][0].longitude;
    cout << "event position: " << latEv << " , " << lonEv << endl;
    aMap->getRelXYFromLatLong(latEv,lonEv,xEv,yEv);
    gEv->SetPoint(ev,xEv,yEv);
    aMap->Fill(latEv,lonEv);
  }

  aMap->addTGraph("baseList","baseList");
  TGraph *gBaseList = aMap->getCurrentTGraph();
  gBaseList->SetMarkerStyle(20); //20=filled circle
  gBaseList->SetMarkerSize(1);
  gBaseList->SetMarkerColor(kRed);

  char* anitaInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  name.str("");
  name << anitaInstallDir << "/share/anitaCalib/baseListA3.root";
  cout << "looking for base list in " << name.str() << endl;
  TFile *fBaseList = TFile::Open(name.str().c_str());
  TTree *baseCampTree = (TTree*)fBaseList->Get("baseCampTree");
  TTree *awsCampTree = (TTree*)fBaseList->Get("awsTree");
  double fullLat,fullLong;
  double fullLataws,fullLongaws;
  baseCampTree->SetBranchAddress("fullLat",&fullLat);
  baseCampTree->SetBranchAddress("fullLong",&fullLong);

  awsCampTree->SetBranchAddress("fullLat",&fullLataws);
  awsCampTree->SetBranchAddress("fullLong",&fullLongaws);

  cout << "baseCampTree->GetEntries() = " << baseCampTree->GetEntries() << endl;
  for (int entry=0; entry<baseCampTree->GetEntries(); entry++) {
    baseCampTree->GetEntry(entry);
    double x,y;
    aMap->getRelXYFromLatLong(fullLat,fullLong,x,y);
    gBaseList->SetPoint(entry,x,y);
  }


  cout << "awsCampTree->GetEntries() = " << awsCampTree->GetEntries() << endl;
  for (int entry=0; entry<awsCampTree->GetEntries(); entry++) {
    awsCampTree->GetEntry(entry);
    double x,y;
    aMap->getRelXYFromLatLong(fullLataws,fullLongaws,x,y);
    gBaseList->SetPoint(gBaseList->GetN(),x,y);
  }


  aMap->setCurrentHistogram("hist");
  aMap->DrawHist("colz");

  new TCanvas();
  aMap->setCurrentTGraph("event");
  aMap->DrawTGraph("pSame");
  gBaseList->Draw("pSame");




  return;
}


