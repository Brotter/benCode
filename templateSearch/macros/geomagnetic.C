/*

  
  I need to determine what the stokes parameters, and therefore the linear polarization angle, is for all these CR events

  Ben Strutt and John Russell wrote it, so I will go ahead and use it


 */
#include "AnitaEventSummary.h"
#include "GeoMagnetic.h"

TGraph* findGeomagneticAngle() {

  TFile *inFile = TFile::Open("candidates.root");
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  AnitaEventSummary *evSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);

  TGraph *gExVsMeas = new TGraph();
  gExVsMeas->SetName("gExVsMeas");
  gExVsMeas->SetTitle("Geomagnetic Polarization, measured vs expected; Expected Polarisation Angle (degrees); Measured Polarisation Angle (degrees)");

  TH1D *hMeas = new TH1D("hMeas","Measured Polarization Angle",46,-45,45);
  TH1D *hExp = new TH1D("hExp","Expected Polarization Angle",46,-45,45);
  TH1D *hDiff = new TH1D("hDiff","Difference between Measured and Expected",91,-90,90);


  for (int entry=0; entry<summaryTree->GetEntries(); entry++) {
    summaryTree->GetEntry(entry);

    //quality cuts
    if (TMath::Abs(evSum->peak[0][0].hwAngle) > 45 || evSum->flags.maxBottomToTopRatio[0] > 6) {
      cout << "ev" << evSum->eventNumber << " cut due to hw or blast" << endl;
      continue;
    }

    cout << "--------------------------------------------------" << endl;

    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);

    double phi = TMath::DegToRad() * evSum->peak[0][0].phi;
    double theta = TMath::DegToRad() * evSum->peak[0][0].theta;

    cout << "Phi: " << phi << " Theta: " << theta << endl;

    double exPol = -1*FFTtools::wrap(TMath::RadToDeg() * GeoMagnetic::getExpectedPolarisation(*usefulGPS,phi,theta),180,0);

    /*exPol = TMath::Abs(exPol);
    if (exPol > 90) exPol = 180-exPol;
    exPol = TMath::Abs(exPol);
    if (exPol > 45) exPol = 90-exPol;
    exPol = TMath::Abs(exPol);*/

    double measPol = evSum->coherent_filtered[0][0].linearPolAngle();
    
    cout << "eventNumber: " << evSum->eventNumber << " exPol: " << exPol << " measPol: " << measPol << endl;

    hMeas->Fill(measPol);
    hExp->Fill(exPol);
    hDiff->Fill(TMath::Abs(measPol-exPol));
    gExVsMeas->SetPoint(entry,exPol,measPol);

    delete usefulGPS;
  }

  TF1 *line = new TF1("line","x",-45,45);
  TF1 *lineP1 = new TF1("line","x+10",-45,45);
  TF1 *lineM1 = new TF1("line","x-10",-45,45);
  lineP1->SetLineStyle(2);
  lineM1->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(2,2);
  c1->cd(1);
  gExVsMeas->Draw("ap");
  line->Draw("lSame");
  lineP1->Draw("lSame");
  lineM1->Draw("lSame");
  c1->cd(2);
  hMeas->Draw("");
  c1->cd(3);
  hExp->Draw("");
  c1->cd(4);
  hDiff->Draw("");

  return gExVsMeas;
}



