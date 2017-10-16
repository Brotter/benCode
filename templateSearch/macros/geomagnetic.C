/*

  
  I need to determine what the stokes parameters, and therefore the linear polarization angle, is for all these CR events

  Ben Strutt and John Russell wrote it, so I will go ahead and use it


 */
#include "AnitaEventSummary.h"
#include "GeoMagnetic.h"





void findGeoMagnetic2(string inFilename="trueCandidates.root") {
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
    gExVsMeas[entry]->SetPointError(gExVsMeas[entry]->GetN()-1,2,4.5);
    if (exPolUp > -9999) {
      gExVsMeasUp[entry]->SetPoint(gExVsMeasUp[entry]->GetN(),exPolUp,measPol);
      gExVsMeasUp[entry]->SetPointError(gExVsMeasUp[entry]->GetN()-1,5,4.5);
      TArrow *currArrow = new TArrow(exPol,measPol,exPolUp,measPol,0.01,"->-");
      arrows.push_back(currArrow);
    }
    delete usefulGPS;
  }

  TF1 *line = new TF1("line","x",-45,45);
  TF1 *lineP1 = new TF1("line","x+5",-45,45);
  TF1 *lineM1 = new TF1("line","x-5",-45,45);
  lineP1->SetLineStyle(2);
  lineM1->SetLineStyle(2);
  TF1 *lineP2 = new TF1("line","x+10",-45,45);
  TF1 *lineM2 = new TF1("line","x-10",-45,45);
  lineP2->SetLineStyle(3);
  lineM2->SetLineStyle(3);
  lineP2->SetLineWidth(1);
  lineM2->SetLineWidth(1);
  

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(2);
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

  TFile *inFile = TFile::Open("candidates.root");
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
    //quality cuts
    if (TMath::Abs(evSum->peak[0][0].hwAngle) > 45 || evSum->flags.maxBottomToTopRatio[0] > 6) {
      cout << "ev" << evSum->eventNumber << " cut due to hw or blast" << endl;
      continue;
    }

    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);

    double phi = TMath::DegToRad() * evSum->peak[0][0].phi;
    double theta = TMath::DegToRad() * evSum->peak[0][0].theta;

    cout << "EventNumber: " << evSum->eventNumber << endl;
    cout << "Phi: " << evSum->peak[0][0].phi << " Theta: " << evSum->peak[0][0].theta << endl;

    double exPol = -1 * TMath::RadToDeg() * GeoMagnetic::getExpectedPolarisation(*usefulGPS,phi,theta);

    /*exPol = TMath::Abs(exPol);
    if (exPol > 90) exPol = 180-exPol;
    exPol = TMath::Abs(exPol);
    if (exPol > 45) exPol = 90-exPol;
    exPol = TMath::Abs(exPol);*/

    double measPol = evSum->coherent_filtered[0][0].linearPolAngle();
    
    cout << "eventNumber: " << evSum->eventNumber << " exPol: " << exPol << " measPol: " << measPol << endl;

    hMeas->Fill(measPol);
    hExp->Fill(exPol);
    hDiff->Fill(measPol-exPol);
    gExVsMeas->SetPoint(gExVsMeas->GetN(),exPol,measPol);
    gExVsMeas->SetPointError(gExVsMeas->GetN()-1,0.1,4.5);
    delete usefulGPS;
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
    

