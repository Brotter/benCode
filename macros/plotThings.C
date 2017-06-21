#include "AnitaConventions.h"


void plotThings() {

  TFile *outFile = TFile::Open("plotThings.root","recreate");

  TChain* summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");


  /*================
    Template Correlation Stuff
  */

  TH2D *noiseTmplt = new TH2D("noiseTmplt","Cosmic Ray Template (5) Correlation - No Pulsers; Interferometric Peak; Template Corr",
			      500,0,1,500,0,1);
  TH2D *waisTmplt = new TH2D("waisTmplt","Cosmic Ray Template (5) Correlation - WAIS; Interferometric Peak; Template Corr",
			     500,0,1,500,0,1);
  TH2D *ldbTmplt = new TH2D("ldbTmplt","Cosmic Ray Template (5) Correlation - LDB; Interferometric Peak; Template Corr",
			    500,0,1,500,0,1);
  TH2D *cutTmplt = new TH2D("cutTmplt","Cosmic Ray Template (5) Correlation - Simple Cut; Interferometric Peak; Template Corr",
			    500,0,1,500,0,1);

  TCanvas *cTemplate = new TCanvas("cTemplate","cTemplate",800,600);
  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> noise",
		    "flags.pulser == 0","colz");
  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> wais",
		    "flags.pulser == 1","colzSame");
  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> ldb",
		    "flags.pulser == 2","colzSame");
  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> ldb",
		    "flags.pulser == 0 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05","colz");

  outFile->cd();
  noiseTmplt->Write();
  waisTmplt->Write();
  ldbTmplt->Write();
  cutTmplt->Write();



  TH1D *h1Tmplt = noiseTmplt->ProjectionY("h1Tmplt",0,500);
  h1Tmplt->SetTitle("Cosmic Ray Template (5) Correlation - No Pulsers;Template Corr; Occupancy");
  TH1D *h1TmpltWais = waisTmplt->ProjectionY("h1TmpltWais",0,500);
  h1TmpltWais->SetTitle("Cosmic Ray Template (5) Correlation - Wais;Template Corr; Occupancy");
  TH1D *h1TmpltLDB = ldbTmplt->ProjectionY("h1TmpltLDB",0,500);
  h1TmpltLDB->SetTitle("Cosmic Ray Template (5) Correlation - LDB;Template Corr; Occupancy");

  TH1D *h1IPeak = noiseTmplt->ProjectionX("h1IPeak",0,500);
  h1IPeak->SetTitle("Interferometric Map Peak - No Pulsers; Map Peak; Occupancy");
  TH1D *h1IPeakWais = waisTmplt->ProjectionX("h1IPeakWais",0,500);
  h1IPeakWais->SetTitle("Interferometric Map Peak - WAIS; Map Peak; Occupancy");
  TH1D *h1IPeakLDB = ldbTmplt->ProjectionX("h1IPeakLDB",0,500);
  h1IPeakLDB->SetTitle("Interferometric Map Peak - LDB; Map Peak; Occupancy");

  outFile->cd();
  h1Tmplt->Write();
  h1TmpltWais->Write();
  h1TmpltLDB->Write();
  h1IPeak->Write();
  h1IPeakWais->Write();
  h1IPeakLDB->Write();

  /*
    --------------*/


  /*=================
    Polarization Plots
  */
  TCanvas *cPolariz = new TCanvas("cPolariz","Polarization Info",800,600);
  TH2D *polariz = new TH2D("polariz","Polarization Info - Noise; Linear Polarization Fraction; Linear Polarization Angle",
			   500,0,1,360,-45,45);
  TH2D *polarizWais = new TH2D("polarizWais","Polarization Info - Wais; Linear Polarization Fraction; Linear Polarization Angle",
			       500,0,1,360,-45,45);
  TH2D *polarizLDB = new TH2D("polarizLDB","Polarization Info - LDB; Linear Polarization Fraction; Linear Polarization Angle",
			      500,0,1,360,-45,45);

    


  summaryTree->Draw("TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2.:TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I >> polariz"
		    "flags.pulser == 0","colz");


  summaryTree->Draw("TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2.:TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I >> polarizWais"
		    "flags.pulser == 1","colzSame");


  summaryTree->Draw("TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2.:TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I >> polarizLDB"
		    "flags.pulser == 2","colzSame");



  outFile->cd();
  polariz->Write();
  polarizWais->Write();
  polarizLDB->Write();
  

  TH1D *linPolFrac = polariz->ProjectionY("linPolFrac",0,500);
  linPolFrac->SetTitle("Linear Polarization Fraction - No Pulsers; Linear Polarization Fraction; Occupancy");
  TH1D *linPolFracWais = polarizWais->ProjectionY("linPolFracWais",0,500);
  linPolFracWais->SetTitle("Linear Polarization Fraction - No Pulsers; Linear Polarization Fraction; Occupancy");
  TH1D *linPolFracLDB = polarizLDB->ProjectionY("linPolFracLDB",0,500);
  linPolFracLDB->SetTitle("Linear Polarization Fraction - No Pulsers; Linear Polarization Fraction; Occupancy");

  TH1D *linPolAng = polariz->ProjectionY("linPolAng",0,360);
  linPolAng->SetTitle("Linear Polarization Angtion - No Pulsers; Linear Polarization Angle; Occupancy");
  TH1D *linPolAngWais = polarizWais->ProjectionY("linPolAngWais",0,360);
  linPolAngWais->SetTitle("Linear Polarization Angtion - No Pulsers; Linear Polarization Angle; Occupancy");
  TH1D *linPolAngLDB = polarizLDB->ProjectionY("linPolAngLDB",0,360);
  linPolAngLDB->SetTitle("Linear Polarization Angtion - No Pulsers; Linear Polarization Angle; Occupancy");


  outFile->cd();
  linPolFrac->Write();
  linPolFracWais->Write();
  linPolFracLDB->Write();

  linPolAng->Write();
  linPolAngWais->Write();
  linPolAngLDB->Write();

  /*
    ---------*/

  
  TCanvas *cCut = new TCanvas("cCut","cCut",800,600);
  TH1D *bothCutsWais = new TH1D("bothCutsWais","bothCutsWais",1000,0,9e7);

  summaryTree->Draw("eventNumber >> bothCutsWais","flags.pulser == 1 && TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 15 && TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05 ","colz");


  TH1D *bothCuts = new TH1D("bothCuts","bothCuts",1000,0,9e7);

  summaryTree->Draw("eventNumber >> bothCuts","flags.pulser == 0 && TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 15 && TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05 && TMath::Abs(peak[0][0].phi - ldb.phi) > 5 && TMath::Abs(peak[0][0].phi - ldb.phi) > 5","colz");


  TH2D *cutLatLon = new TH2D("cutLatLon","cutLatLon",250,-90,-65,360,-180,180);

  summaryTree->Draw("peak[0][0].longitude:peak[0][0].latitude >> cutLatLon","flags.pulser == 0 && TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 15 && TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05 && TMath::Abs(peak[0][0].phi - ldb.phi) > 5 && TMath::Abs(peak[0][0].phi - ldb.phi) > 5 && peak[0][0].latitude > -999 && !( (peak[0][0].latitude < -78 && peak[0][0].latitude > -82) && (peak[0][0].longitude < -105 && peak[0][0].longitude > -115) )","colz");

  TFile *baseList = TFile::Open("~/anita16/local/share/anitaCalib/baseListA3.root");
  TTree *baseTree = (TTree*)baseList->Get("baseCampTree");
  TTree *awsTree  = (TTree*)baseList->Get("awsTree");

  baseTree->SetMarkerStyle(3);
  baseTree->SetMarkerColor(kRed);
  awsTree->SetMarkerStyle(3);
  awsTree->SetMarkerColor(kBlue);

  baseTree->Draw("fullLong:fullLat","","same");
  awsTree->Draw("fullLong:fullLat","","same");
  
    


  outFile->Close();

  return;
}
