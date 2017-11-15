/*

  
  I need to determine what the stokes parameters, and therefore the linear polarization angle, is for all these CR events

  Ben Strutt and John Russell wrote it, so I will go ahead and use it


 */
#include "AnitaEventSummary.h"
#include "GeoMagnetic.h"


double getGeomagExpected(AnitaEventSummary *evSum, Adu5Pat *gps, bool upgoing=false) {
  /*
    For integration into other code:
    
    Inputs:
    evSum - an event summary that gives stokes info and pointing info
    gps - an Adu5Pat that gives location and heading info
    
    Outputs:
    expected - from GeoMagnetic::getExpected
    
    optional:
    upgoing - maybe you want to see an upgoing shower instead
    
  */
  
  UsefulAdu5Pat usefulGPS(gps);
  
  double phi = TMath::DegToRad() * evSum->peak[0][0].phi;
  double theta = TMath::DegToRad() * evSum->peak[0][0].theta;
  
  double exPol;
  if (!upgoing) exPol = TMath::RadToDeg() * GeoMagnetic::getExpectedPolarisation(usefulGPS,phi,theta);
  else          exPol = TMath::RadToDeg() * GeoMagnetic::getExpectedPolarisationUpgoing(usefulGPS,phi,theta,1000);
  
  
  return exPol;
  
}



void saveGeoMagnetic(string inFileName="cutsClust_oct14.root") {
  /*

    Read in a summaryTree file, then add an expGeoPol field to it

    Outputs the same filename with "_geo" appended before the .root

   */


  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add(inFileName.c_str());

  int lenEntries = summaryTree->GetEntries();
  cout << " Found " << lenEntries << " events" << endl;


  size_t pos = inFileName.find(".root");
  string baseName = inFileName.substr(0,pos);
  stringstream outFileName;
  outFileName << baseName << "_geo.root";

  TFile *outFile = TFile::Open(outFileName.str().c_str(),"recreate");
  TTree *outTree = new TTree("summaryTree","summaryTree");

  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *tempSum = NULL;
  AnitaNoiseSummary *noiseSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  outTree->Branch("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);
  outTree->Branch("template",&tempSum);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  outTree->Branch("noiseSummary",&noiseSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  outTree->Branch("gpsEvent",&gps);
  double clusterValue;
  summaryTree->SetBranchAddress("clusterValue",&clusterValue);
  outTree->Branch("clusterValue",&clusterValue);
  double exGeoPol;
  outTree->Branch("exGeoPol",&exGeoPol);


  for (int entry=0; entry<lenEntries; entry++) {
    if (!(entry%100)) cout << entry << "/" << lenEntries << endl;

    summaryTree->GetEntry(entry);

    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);
    double phi = TMath::DegToRad() * evSum->peak[0][0].phi;
    double theta = TMath::DegToRad() * evSum->peak[0][0].theta;
    exGeoPol = TMath::RadToDeg() * GeoMagnetic::getExpectedPolarisation(*usefulGPS,phi,theta);

    outTree->Fill();

    delete usefulGPS;

  }
    
  outFile->cd();
  outTree->Write();
  outFile->Close();

  return;

}


void plotGeoMagneticLabeled(string inFilename="cutsClust_oct14_geo_labeled.root") {
/*

  Take in the 5997 events that passed my "impulsivity" cuts, then make the geomagnetic graphs, appropriately labeled, for them

  needs to have a "_geo" and a "_labeled"

*/
  stringstream name;

  //  GeoMagnetic::setDebug(true);

  TFile *inFile = TFile::Open(inFilename.c_str());
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  int lenEntries = summaryTree->GetEntries();
  cout << "Found " << lenEntries << " entries in " << inFilename << endl;

  AnitaEventSummary *evSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  TString *labelString = NULL;
  summaryTree->SetBranchAddress("label",&labelString);
  double exGeoPol;
  summaryTree->SetBranchAddress("exGeoPol",&exGeoPol);


  //the collection of TGraphErrors for the candidates
  vector<TGraphErrors*> vCandidateGeoMap;

  //the collection of TGraphErrors for the candidates and the dirty dozen
  vector<TGraphErrors*> vStrongGeoMap;




  //A histogram for the measured values for strong events
  TH1D *hStrongMeas = new TH1D("hStrongMeas","Measured Polarization Angle: Strong events;degrees;count",46,-45,45);
  TH1D *hStrongExpe = new TH1D("hStrongExpe","Expected Polarization Angle: Strong events;degrees;count",46,-45,45);
  TH1D *hStrongDiff = new TH1D("hStrongDiff","Difference between Measured and Expected : Strong events;degrees;count",91,-45,45);
  TH2D *h2Strong = new TH2D("h2Strong","Strong Events; Expected Polarization Angle (degrees);Measured Polarization Angle (degrees)",
			    91,-45,45, 46,-45,45);

  //A histogram for the measured values for all impulsive events
  TH1D *hImpulsMeas = new TH1D("hImpulsMeas","Measured Polarization Angle: Impulsive events;degrees;count",46,-45,45);
  TH1D *hImpulsExpe = new TH1D("hImpulsExpe","Expected Polarization Angle: Impulsive events;degrees;count",46,-45,45);
  TH1D *hImpulsDiff = new TH1D("hImpulsDiff","Difference between Measured and Expected : Impulsive events;degrees;count",91,-45,45);
  TH2D *h2Impuls = new TH2D("h2Impuls","Impuls Events;Expected Polarization Angle (degrees); Measured Polarization Angle (degrees)",
			    91,-45,45, 46,-45,45);


  for (int entry=0; entry<lenEntries; entry++) {
    if (!(entry%10)) cout << entry << "/" << lenEntries << endl;

    summaryTree->GetEntry(entry);


    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);

    double phi = TMath::DegToRad() * evSum->peak[0][0].phi;
    double theta = TMath::DegToRad() * evSum->peak[0][0].theta;

    double measPol = evSum->coherent_filtered[0][0].linearPolAngle();

    //candidates
    if (strstr(labelString->Data(),"Below Horizon Candidate") || strstr(labelString->Data(),"Above Horizon Candidate")) {
      TGraphErrors *geCurr = new TGraphErrors();
      geCurr->SetPoint(0,exGeoPol,measPol);
      geCurr->SetPointError(0,2.,6.);
      if (strstr(labelString->Data(),"Below Horizon Candidate")) { geCurr->SetMarkerColor(kBlue);geCurr->SetLineColor(kBlue); }
      if (strstr(labelString->Data(),"Above Horizon Candidate")) { geCurr->SetMarkerColor(kRed); geCurr->SetLineColor(kRed); }
      string sName = "ev"+to_string(evSum->eventNumber);
      geCurr->SetName(sName.c_str());
      geCurr->SetTitle(labelString->Data());

      geCurr->GetXaxis()->SetRangeUser(-30,30);
      geCurr->GetYaxis()->SetRangeUser(-30,30);

      vCandidateGeoMap.push_back(geCurr);
      vStrongGeoMap.push_back(geCurr);
    }

    //dirty dozen
    if (strstr(labelString->Data(),"Dirty Dozen")) {
      TGraphErrors *geCurr = new TGraphErrors();
      geCurr->SetPoint(0,exGeoPol,measPol);
      geCurr->SetPointError(0,2.,6.);
      string sName = "ev"+to_string(evSum->eventNumber);
      geCurr->SetName(sName.c_str());
      geCurr->SetTitle(labelString->Data());

      geCurr->GetXaxis()->SetRangeUser(-30,30);
      geCurr->GetYaxis()->SetRangeUser(-30,30);

      vStrongGeoMap.push_back(geCurr);
    }
      

    //strong non-candidates (clustered or above horizon)
    if (strstr(labelString->Data(),"Above Horizon Passing") || strstr(labelString->Data(),"Clustered Passing")) {
      hStrongMeas->Fill(measPol);
      hStrongExpe->Fill(exGeoPol);
      hStrongDiff->Fill(measPol-exGeoPol);
      h2Strong->Fill(exGeoPol,measPol);

    }


    //all impulsive non-candidates (clustered or above horizon)
    if (strstr(labelString->Data(),"Above Horizon Failing") || strstr(labelString->Data(),"Clustered Failing")) {
      hImpulsMeas->Fill(measPol);
      hImpulsExpe->Fill(exGeoPol);
      hImpulsDiff->Fill(measPol-exGeoPol);
      h2Impuls->Fill(exGeoPol,measPol);      

    }

    delete usefulGPS;
  }

  cout << "Number of Candidates " << vCandidateGeoMap.size() << endl;
  cout << "Number of Isolated Events" << vStrongGeoMap.size() << endl;



  TF1 *line = new TF1("line","x",-45,45);
  TF1 *lineP1 = new TF1("line","x+6",-45,45);
  TF1 *lineM1 = new TF1("line","x-6",-45,45);
  lineP1->SetLineStyle(2);
  lineM1->SetLineStyle(2);
  TF1 *lineP2 = new TF1("line","x+12",-45,45);
  TF1 *lineM2 = new TF1("line","x-12",-45,45);
  lineP2->SetLineStyle(3);
  lineM2->SetLineStyle(3);
  lineP2->SetLineWidth(1);
  lineM2->SetLineWidth(1);

  TCanvas *cCands = new TCanvas("cCands","cCands",1000,5000);

  line->Draw();
  for (int i=0; i<vCandidateGeoMap.size(); i++) {
    vCandidateGeoMap[i]->Draw("psame");

  }


  TCanvas *cDD = new TCanvas("cDD","cDD",1000,5000);
  line->Draw();
  for (int i=0; i<vStrongGeoMap.size(); i++) {
    vStrongGeoMap[i]->Draw("psame");

  }

  TCanvas *cStrong = new TCanvas("cStrong","cStrong",1000,5000);
  cStrong->Divide(2,2);
  cStrong->cd(1);
  h2Strong->Draw("colz");
  cStrong->cd(2);
  hStrongMeas->Draw();
  cStrong->cd(3);
  hStrongExpe->Draw();
  cStrong->cd(4);
  hStrongDiff->Draw();

  TCanvas *cImpuls = new TCanvas("cImpuls","cImpuls",1000,5000);
  cImpuls->Divide(2,2);
  cImpuls->cd(1);
  h2Impuls->Draw("colz");
  cImpuls->cd(2);
  hImpulsMeas->Draw();
  cImpuls->cd(3);
  hImpulsExpe->Draw();
  cImpuls->cd(4);
  hImpulsDiff->Draw();
  
    
  TFile *outRootFile = TFile::Open("findGeomagneticLabeled.root","recreate");
  hStrongMeas->Write();
  hImpulsMeas->Write();
  hStrongExpe->Write();
  hImpulsExpe->Write();
  hStrongDiff->Write();
  hImpulsDiff->Write(); 
  for (int i=0; i<vStrongGeoMap.size(); i++) {
    vStrongGeoMap[i]->Write();
  }
  for (int i=0; i<vCandidateGeoMap.size(); i++) {
    vCandidateGeoMap[i]->Write();
  }
  outRootFile->Close();


  return;
}


void findGeoMagnetic2(string inFilename="cuts_oct14.root") {
/*
  Do both upgoing and direct, and include the direct events as well.
*/
  stringstream name;

  //  GeoMagnetic::setDebug(true);

  TFile *inFile = TFile::Open(inFilename.c_str());
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  int lenEntries = summaryTree->GetEntries();
  cout << "Found " << lenEntries << " entries in " << inFilename << endl;

  AnitaEventSummary *evSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);


  TGraphErrors **gExVsMeas = (TGraphErrors**)malloc(sizeof(TGraphErrors*)*lenEntries);
  TGraphErrors **gExVsMeasUp = (TGraphErrors**)malloc(sizeof(TGraphErrors*)*lenEntries);


  TH1D *hMeas = new TH1D("hMeas","Measured Polarization Angle;degrees;count",46,-45,45);

  TH1D *hExp = new TH1D("hExp","Expected Polarization Angle;degrees;count",46,-45,45);
  TH1D *hDiff = new TH1D("hDiff","Difference between Measured and Expected;degrees;count",91,-45,45);

  TH1D *hExpUp = new TH1D("hExpUp","Expected Polarization Angle;degrees;count",46,-45,45);
  TH1D *hDiffUp = new TH1D("hDiffUp","Difference between Measured and Expected;degrees;count",91,-45,45);

  vector<TArrow*> arrows;

  for (int entry=0; entry<summaryTree->GetEntries(); entry++) {
    summaryTree->GetEntry(entry);
    //    if (evSum->eventNumber != 15717147) continue;


    

    name.str("");
    name << "ev" << evSum->eventNumber;
    gExVsMeas[entry] = new TGraphErrors();
    gExVsMeas[entry]->SetName(name.str().c_str());
    gExVsMeas[entry]->SetTitle("Reflected: Geomagnetic Polarization; Expected Polarisation Angle (degrees); Measured Polarisation Angle (degrees)");
    
    


    name << "_up";
    gExVsMeasUp[entry] = new TGraphErrors();    
    gExVsMeasUp[entry]->SetName(name.str().c_str());
    gExVsMeasUp[entry]->SetTitle("Direct (Upgoing): Geomagnetic Polarization; Expected Polarisation Angle (degrees); Measured Polarisation Angle (degrees)");
    gExVsMeasUp[entry]->SetMarkerColor(kBlue);
    gExVsMeasUp[entry]->SetLineColor(kBlue);

    cout << "--------------------------------------------------" << endl;

    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);

    //    if (evSum->eventNumber == 78011657 || evSum->eventNumber == 80299372) continue; //HiCal
    //    if (evSum->eventNumber == 56444667 || evSum->eventNumber == 57475797) continue; //Um, base things?  :(

    double phi = TMath::DegToRad() * evSum->peak[0][0].phi;
    double theta = TMath::DegToRad() * evSum->peak[0][0].theta;

    cout << "findGeomagnetic2(): EventNumber=" << evSum->eventNumber << endl;
    cout << "Phi: " << evSum->peak[0][0].phi << " Theta: " << evSum->peak[0][0].theta << endl;

    double exPol = TMath::RadToDeg() * GeoMagnetic::getExpectedPolarisation(*usefulGPS,phi,theta);

    double exPolUp = -9999;
    if (evSum->peak[0][0].latitude != -9999) {
      exPolUp = TMath::RadToDeg() * GeoMagnetic::getExpectedPolarisationUpgoing(*usefulGPS,phi,theta,1000);
      //      if (exPolUp > 90) exPolUp = 180-exPolUp;
      //      if (exPolUp < -90) exPolUp = 180+exPolUp;

    }
    /*exPol = TMath::Abs(exPol);
    if (exPol > 90) exPol = 180-exPol;
    exPol = TMath::Abs(exPol);
    if (exPol > 45) exPol = 90-exPol;
    exPol = TMath::Abs(exPol);*/

    double measPol = evSum->coherent_filtered[0][0].linearPolAngle();
    
    cout << "measPol: " << measPol << " exPol: " << exPol << " exPolUp:" << exPolUp << endl;

    hMeas->Fill(measPol);
    hExp->Fill(exPol);
    hExpUp->Fill(exPolUp);
    hDiff->Fill(measPol-exPol);
    hDiffUp->Fill(measPol-exPolUp);
    gExVsMeas[entry]->SetPoint(gExVsMeas[entry]->GetN(),exPol,measPol);
    gExVsMeas[entry]->SetPointError(gExVsMeas[entry]->GetN()-1,2.,6.);
    if (exPolUp > -9999) {
      gExVsMeasUp[entry]->SetPoint(gExVsMeasUp[entry]->GetN(),exPolUp,measPol);
      gExVsMeasUp[entry]->SetPointError(gExVsMeasUp[entry]->GetN()-1,5,4.5);
      TArrow *currArrow = new TArrow(exPol,measPol,exPolUp,measPol,0.01,"->-");
      arrows.push_back(currArrow);
    }
    delete usefulGPS;
  }

  TF1 *line = new TF1("line","x",-45,45);
  TF1 *lineP1 = new TF1("line","x+6",-45,45);
  TF1 *lineM1 = new TF1("line","x-6",-45,45);
  lineP1->SetLineStyle(2);
  lineM1->SetLineStyle(2);
  TF1 *lineP2 = new TF1("line","x+12",-45,45);
  TF1 *lineM2 = new TF1("line","x-12",-45,45);
  lineP2->SetLineStyle(3);
  lineM2->SetLineStyle(3);
  lineP2->SetLineWidth(1);
  lineM2->SetLineWidth(1);
  

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  //  c1->Divide(2);
  c1->cd(1);
  gExVsMeas[0]->Draw("ap"); 
  gExVsMeas[0]->GetXaxis()->SetLimits(-30,30);
  gExVsMeas[0]->GetYaxis()->SetRangeUser(-30,30);
  line->Draw("lSame");
  lineP1->Draw("lSame");
  lineM1->Draw("lSame");
  lineP2->Draw("lSame");
  lineM2->Draw("lSame");
  for (int i=0; i<lenEntries; i++) {
    gExVsMeas[i]->Draw("psame");
  }
  /*
  c1->cd(2);
  gExVsMeas[0]->Draw("ap"); 
  gExVsMeas[0]->GetXaxis()->SetLimits(-30,30);
  gExVsMeas[0]->GetYaxis()->SetRangeUser(-30,30);
  line->Draw("lSame");
  lineP1->Draw("lSame");
  lineM1->Draw("lSame");
  lineP2->Draw("lSame");
  lineM2->Draw("lSame");
  for (int i=0; i<lenEntries; i++) {
    gExVsMeasUp[i]->Draw("psame");
  }
  */
  //  for (int i=0; i<arrows.size(); i++) arrows[i]->Draw();
    
  return;
}






double* findExpectedError(int evNum, string inFileName="",bool plots=false) {
  /*
    Idea instead of math:
    Take a single event, then vary the angle by the angular uncertainty and the xMax by the xMax uncertainty
    -> see what the distribution is.

    Probably do it for all the events that you are testing, so have evNum as an input

    if inFileName is filled, it loads that file looking for a summaryTree to get the event from.
    Otherwise it loads the whole thing.
  */
  
  TChain *summaryTree;
  if (inFileName == "") {
    summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");
  }
  else {
    summaryTree = new TChain("summaryTree","summaryTree");
    summaryTree->Add(inFileName.c_str());
  }

  summaryTree->BuildIndex("eventNumber");
  int entry=summaryTree->GetEntryNumberWithIndex(evNum);
  if (entry==-1) {
    cout << "Couldn't find event number " << evNum << ". Quitting" << endl;
    return NULL;
  }

  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("gpsEvent",&gps);

  summaryTree->GetEntry(entry);

  //get anita location
  UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);

  //get central pointing hypothesis
  double centralPhi = evSum->peak[0][0].phi*TMath::DegToRad();
  double centralTheta = evSum->peak[0][0].theta*TMath::DegToRad();
  
  const double phiSigma = 0.5*TMath::DegToRad();
  const double thetaSigma = 0.25*TMath::DegToRad();

  double phis[3] = {centralPhi-phiSigma,centralPhi,centralPhi+phiSigma};
  double thetas[3] = {centralTheta-thetaSigma,centralTheta,centralTheta+thetaSigma};

  double xMaxHigh = GeoMagnetic::xMaxFe18;
  double xMaxLow = GeoMagnetic::xMaxP19;
  double xMaxCent = (xMaxHigh + xMaxLow)/2.;
    
  double xMaxesDn[3] = {xMaxLow,xMaxCent,xMaxHigh};

  double xMaxesUp[3] = {200,5000,10000};

  //three things to vary, sigma, phi and xMax, each with three values, leaves 3^3=27
  double expectedDn[27];
  double expectedUp[27];
  TH1D *hExpDn = new TH1D("hExpDn","Expected polarization downgoing",5001,-5,5);
  TH1D *hExpUp = new TH1D("hExpUp","Expected polarization upgoing",5001,-5,5);

  TGraph *gExpDn = new TGraph();
  TGraph *gExpUp = new TGraph();

  for (int phii=0; phii<3; phii++) {
    for (int thetai=0; thetai<3; thetai++) {
      for (int xMaxi=0; xMaxi<3; xMaxi++) {
	int index = phii*3*3 + thetai*3 + xMaxi;
	expectedUp[index] = GeoMagnetic::getExpectedPolarisationUpgoing(*usefulGPS,phis[phii],
									thetas[thetai],xMaxesUp[xMaxi]);
	hExpUp->Fill(expectedUp[index]);
	gExpUp->SetPoint(index,index,expectedUp[index]);
	expectedDn[index] = GeoMagnetic::getExpectedPolarisation(*usefulGPS,phis[phii],
								 thetas[thetai],xMaxesDn[xMaxi]);
	hExpDn->Fill(expectedDn[index]);
	gExpDn->SetPoint(index,index,expectedDn[index]);
	cout << index << " " << expectedUp[index] << " " << expectedDn[index] << endl;
      }
    }
  }

  if (plots) {
    TCanvas *cDn = new TCanvas("cDn","",1000,600);
    cDn->Divide(2);
    cDn->cd(1);
    gExpDn->Draw("alp");
    cDn->cd(2);
    hExpDn->Draw("");
    
    TCanvas *cUp = new TCanvas("cUp","",1000,600);
    cUp->Divide(2);
    cUp->cd(1);
    gExpUp->Draw("alp");
    cUp->cd(2);
    hExpUp->Draw("");
  }

  //dnMin,dnMax,upMin,upMax
  double *returnVal = new double[4];
  returnVal[0] = TMath::MinElement(gExpDn->GetN(),gExpDn->GetY());
  returnVal[1] = TMath::MaxElement(gExpDn->GetN(),gExpDn->GetY());
  returnVal[2] = TMath::MinElement(gExpUp->GetN(),gExpUp->GetY());
  returnVal[3] = TMath::MaxElement(gExpUp->GetN(),gExpUp->GetY());

  return returnVal;
}



TGraph* findGeoMagneticAngle() {

  GeoMagnetic::setDebug(true);

  TFile *inFile = TFile::Open("trueCandidates_oct14_reMasked.root");
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  AnitaEventSummary *evSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);

  

  TGraphErrors *gExVsMeas = new TGraphErrors();
  gExVsMeas->SetName("gExVsMeas");
  gExVsMeas->SetTitle("Geomagnetic Polarization, measured vs expected; Expected Polarisation Angle (degrees); Measured Polarisation Angle (degrees)");

  TH1D *hMeas = new TH1D("hMeas","Measured Polarization Angle;degrees;count",46,-45,45);
  TH1D *hExp = new TH1D("hExp","Expected Polarization Angle;degrees;count",46,-45,45);
  TH1D *hDiff = new TH1D("hDiff","Difference between Measured and Expected;degrees;count",91,-45,45);


  for (int entry=0; entry<summaryTree->GetEntries(); entry++) {
    summaryTree->GetEntry(entry);


    cout << "--------------------------------------------------" << endl;
    cout << "EventNumber: " << evSum->eventNumber << endl;
    cout << "Phi: " << evSum->peak[0][0].phi << " Theta: " << evSum->peak[0][0].theta << endl;

    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);
    double phi = TMath::DegToRad() * evSum->peak[0][0].phi;
    double theta = TMath::DegToRad() * evSum->peak[0][0].theta;
    double exPol = TMath::RadToDeg() * GeoMagnetic::getExpectedPolarisation(*usefulGPS,phi,theta);
    delete usefulGPS;

    double measPol = evSum->coherent_filtered[0][0].linearPolAngle();
    
    cout << "End: ev" << evSum->eventNumber << " exPol: " << exPol << " measPol: " << measPol << endl;

    hMeas->Fill(measPol);
    hExp->Fill(exPol);
    hDiff->Fill(measPol-exPol);
    gExVsMeas->SetPoint(gExVsMeas->GetN(),exPol,measPol);
    gExVsMeas->SetPointError(gExVsMeas->GetN()-1,2,4.5);

  }

  TF1 *line = new TF1("line","x",-45,45);
  TF1 *lineP1 = new TF1("line","x+10",-45,45);
  TF1 *lineM1 = new TF1("line","x-10",-45,45);
  lineP1->SetLineStyle(2);
  lineM1->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(2,2);
  TVirtualPad *p1 = c1->cd(1);
  p1->SetGridx();
  p1->SetGridy();
  gExVsMeas->Draw("ap");
  line->Draw("lSame");
  //  lineP1->Draw("lSame");
  //  lineM1->Draw("lSame");
  c1->cd(2);
  hMeas->Draw("");
  c1->cd(3);
  hExp->Draw("");
  c1->cd(4);
  hDiff->Draw("");

  return gExVsMeas;
}





void findGeoMagneticNoise() {
  /*
    Whats that anthropogenic sample look like?
   */
  
  TFile *inFile = TFile::Open("mergeBaseClusters.root");
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  
  TFile *outFile = TFile::Open("geomagneticNoise.root","recreate");
  TTree *geomagTree = new TTree("geomagTree","geomagTree");
  int eventNumber;
  double expectUp,expectDn,measured;
  geomagTree->Branch("eventNumber",&eventNumber);
  geomagTree->Branch("expectUp",&expectUp);
  geomagTree->Branch("expectDn",&expectDn);
  geomagTree->Branch("measured",&measured);



  AnitaEventSummary *evSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);

  TH1D *hDiffDn = new TH1D("hDiffDn","Difference between Measured and Expected - Downgoing;degrees;count",181,-90,90);
  TH1D *hDiffUp = new TH1D("hDiffUp","Difference between Measured and Expected - Upgoing;degrees;count",361,-180,180);

  int lenEntries = summaryTree->GetEntries();
  lenEntries = 10000;
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%100 == 0) cout << entry << "/" << lenEntries << endl;
    summaryTree->GetEntry(entry);
    
    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);
    
    eventNumber = evSum->eventNumber;

    measured = evSum->coherent_filtered[0][0].linearPolAngle();

    double phi = TMath::DegToRad() * evSum->peak[0][0].phi;
    double theta = TMath::DegToRad() * evSum->peak[0][0].theta;

    //If the event is assumed to be "downgoing," ie reflected
    expectDn = GeoMagnetic::getExpectedPolarisation(*usefulGPS,phi,theta);
    //returns a value that is in radians, is inverted, and needs to be wrapped
    // I don't remember why it needs to be inverted though? maybe not do that
    cout << "1:" << expectDn << " ";
    expectDn = FFTtools::wrap(TMath::RadToDeg() * expectDn,180,0);
    cout << expectDn << endl;
    hDiffDn->Fill(expectDn-measured);

    //do the upgoing estimate too, but only if it is pointed at the continent (which it will be all the time I think)
    if (evSum->peak[0][0].latitude >  -999) {
      expectUp = GeoMagnetic::getExpectedPolarisationUpgoing(*usefulGPS,phi,theta,1000);
      expectUp = FFTtools::wrap(TMath::RadToDeg() * expectUp,180,0);
      hDiffUp->Fill(expectUp-measured);
    }
    else {
      cout << evSum->peak[0][0].latitude << endl;
      cout << "skipping... theta=" << theta*TMath::RadToDeg() << endl;
      expectUp = -999;
      measured = -999;
    }

    geomagTree->Fill();

    delete usefulGPS;

  }

  outFile->Write();
  outFile->Close();


  TCanvas *c1 = new TCanvas("c1","",1000,600);
  c1->Divide(1,2);
  c1->cd(1);
  hDiffDn->Draw();
  c1->cd(2);
  hDiffUp->Draw();
  

  return;
}

  //  crab = code = pincing = potato = code again    
    

/*
  
  Default

 */

void geomagnetic() {
  
  cout << "loaded geomagnetic.C" << endl;
  return;

}
