#include "AnitaConventions.h"


void plotThings() {

  
  TChain* summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");
  
  TH2D *noise = new TH2D("noise","noise",500,0,1,500,0,1);
  TH2D *wais = new TH2D("wais","wais",500,0,1,500,0,1);
  TH2D *ldb = new TH2D("ldb","ldb",500,0,1,500,0,1);


  TCanvas *cTemplate = new TCanvas("cTemplate","cTemplate",800,600);
  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> noise","flags.pulser == 0","colz");
  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> wais","flags.pulser == 1","colzSame");
  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> ldb","flags.pulser == 2","colzSame");
  
  
  TH2D *cut = new TH2D("cut","cut",500,0,1,500,0,1);
  
  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> ldb",
		    "flags.pulser == 0 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05","colz");

  TCanvas *cPolariz = new TCanvas("cPolariz","cPolariz",800,600);
  TH2D *polariz = new TH2D("polariz","polairiz",500,0,1,360,-45,45);
  TH2D *polarizWais = new TH2D("polarizWais","polarizWais",500,0,1,360,-45,45);
  TH2D *polarizLDB = new TH2D("polarizLDB","polarizLDB",500,0,1,360,-45,45);

    


  summaryTree->Draw("TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2.:TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I >> polariz"
		    "flags.pulser == 0","colz");


  summaryTree->Draw("TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2.:TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I >> polarizWais"
		    "flags.pulser == 1","colzSame");


  summaryTree->Draw("TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2.:TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I >> polarizLDB"
		    "flags.pulser == 2","colzSame");



  TH2D *polarizCut = new TH2D("polarizCut","polarizCut",500,0,1,360,-45,45);

  summaryTree->Draw("TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2.:TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I >> polarizCut"
		    "flags.pulser == 2 && TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 15 && TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4","colz");



  TH1D *bothCutsWais = new TH1D("bothCutsWais","bothCutsWais",1000,0,9e7);

  summaryTree->Draw("eventNumber >> bothCutsWais","flags.pulser == 1 && TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 15 && TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05 ","colz");


  TH1D *bothCuts = new TH1D("bothCuts","bothCuts",1000,0,9e7);

  summaryTree->Draw("eventNumber >> bothCuts","flags.pulser == 0 && TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 15 && TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05 ","colz");


  return;
}
