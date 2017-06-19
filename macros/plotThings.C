#include "AnitaConventions.h"


void plotThings() {

  
  gROOT->ProcessLine(".x loadAll.C");
  
  TH2D *noise = new TH2D("noise","noise",500,0,1,500,0,1);
  TH2D *wais = new TH2D("wais","wais",500,0,1,500,0,1);
  TH2D *ldb = new TH2D("ldb","ldb",500,0,1,500,0,1);


  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> noise","flags.pulser == 0","colz");
  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> wais","flags.pulser == 1","colz");
  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> ldb","flags.pulser == 2","colz");
  
  
  TH2D *cut = new TH2D("cut","cut",500,0,1,500,0,1);
  
  summaryTree->Draw("templateCRayH[5]:peak[0][0].value >> ldb",
		    "flags.pulser == 0 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05","colz");

  
  TH2D *polariz = new TH2D("polariz","polairiz",500,0,1,360,-45,45);
  TH2D *polarizWais = new TH2D("polarizWais","polarizWais",500,0,1,360,-45,45);
  TH2D *polarizLDB = new TH2D("polarizLDB","polarizLDB",500,0,1,360,-45,45);

    


  summaryTree->Draw("TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2.:TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I >> polariz"
		    "flags.pulser == 0","colz");


  summaryTree->Draw("TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2.:TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I >> polarizWais"
		    "flags.pulser == 1","colz");


  summaryTree->Draw("TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2.:TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I >> polarizLDB"
		    "flags.pulser == 2","colz");



  TH2D *polarizCut = new TH2D("polarizCut","polarizCut",500,0,1,360,-45,45);

  summaryTree->Draw("TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2.:TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I >> polarizCut"
		    "flags.pulser == 2 && TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 15 && TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4","colz");


  TH1D *bothCuts = new TH1D("bothCuts","bothCuts",10000,0,4e6);

  summaryTree->Draw("Entry$ >> bothCuts","flags.pulser == 1 && TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 15 && TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05 ","colz");



  return;
}
