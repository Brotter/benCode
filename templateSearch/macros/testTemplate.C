/*

  I want to do a correlation of the template vs a bunch of stuff and see what the various values are.


  Like for example vs an averaged wais pulse, or vs various off angles from the coherence angle

  this script sort of just plays around with the idea and sees how various things compare

 */


#include "FFTtools.h"
#include "AnitaConventions.h"
void compareTemplates(){

  TFile *waisFile = TFile::Open("/Users/brotter/anita16/benPrograms/waisPulses/waisImpulseResponse_wAngs.root");

  TGraph *waisPulse = (TGraph*)waisFile->Get("wais01BH");
  cout << "Length: " << waisPulse->GetN() << endl;
  double waisSum = 0;
  for (int pt=0; pt<waisPulse->GetN(); pt++) waisSum += pow(waisPulse->GetY()[pt],2);
  //  waisSum /= waisPulse->GetN();
  double* waisGetCorr = FFTtools::getCorrelation(waisPulse->GetN(),waisPulse->GetY(),waisPulse->GetY());
  TGraph *waisAutoCorr = FFTtools::getCorrelationGraph(waisPulse,waisPulse);
  double waisMax = TMath::MaxElement(waisAutoCorr->GetN(),waisAutoCorr->GetY());
  for (int pt=0; pt<waisPulse->GetN(); pt++) waisPulse->GetY()[pt] /= TMath::Sqrt(waisMax);
  delete waisAutoCorr;
  waisAutoCorr = FFTtools::getCorrelationGraph(waisPulse,waisPulse);


  TGraph *waisPulse2 = (TGraph*)waisFile->Get("wais02BH");
  cout << "Length: " << waisPulse2->GetN() << endl;
  double waisSum2 = 0;
  for (int pt=0; pt<waisPulse2->GetN(); pt++) waisSum2 += pow(waisPulse2->GetY()[pt],2);
  //  waisSum2 /= waisPulse2->GetN();
  TGraph *waisAutoCorr2 = FFTtools::getCorrelationGraph(waisPulse2,waisPulse2);
  double waisMax2 = TMath::MaxElement(waisAutoCorr2->GetN(),waisAutoCorr2->GetY());
  for (int pt=0; pt<waisPulse2->GetN(); pt++) waisPulse2->GetY()[pt] /= TMath::Sqrt(waisMax2);
  delete waisAutoCorr2;
  waisAutoCorr2 = FFTtools::getCorrelationGraph(waisPulse2,waisPulse2);

  TGraph *impulseResponse = new TGraph("/Users/brotter/anita16/local/share/UCorrelator/responses/SingleBRotter/all.imp");
  TGraph *iRPadded = FFTtools::padWaveToLength(impulseResponse,waisPulse->GetN());
  double iRSum = 0;
  for (int pt=0; pt<iRPadded->GetN(); pt++) iRSum += pow(iRPadded->GetY()[pt],2);
  iRSum /= iRPadded->GetN();
  TGraph *iRAutoCorr = FFTtools::getCorrelationGraph(iRPadded,iRPadded);
  double iRMax = TMath::MaxElement(iRAutoCorr->GetN(),iRAutoCorr->GetY());
  for (int pt=0; pt<iRPadded->GetN(); pt++) iRPadded->GetY()[pt] /= TMath::Sqrt(iRMax);
  delete iRAutoCorr;
  iRAutoCorr = FFTtools::getCorrelationGraph(iRPadded,iRPadded);
  
  TGraph *iRPaddedSmall = (TGraph*)iRPadded->Clone();
  for (int pt=0; pt<iRPaddedSmall->GetN(); pt++) iRPaddedSmall->GetY()[pt] /= 1000;
  
  TGraph *corr = FFTtools::getCorrelationGraph(waisPulse,iRPadded);
  TGraph *corr2 = FFTtools::getCorrelationGraph(waisPulse2,iRPadded);
  TGraph *corrWais = FFTtools::getCorrelationGraph(waisPulse,waisPulse2);

  TGraph *corrN = FFTtools::getNormalisedCorrelationGraph(waisPulse,iRPaddedSmall);
  corrN->SetMarkerColor(kRed);

  cout << TMath::MaxElement(waisPulse->GetN(),waisGetCorr) << endl;
  cout << waisMax << " " << waisMax2 << " " << iRMax << endl;
  cout << waisSum << " " << waisSum2 << " " << iRSum << endl;

  TCanvas* c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(3,3);
  c1->cd(1);
  corr->Draw("alp");
  corr->GetXaxis()->SetRangeUser(12,22);
  corr->SetTitle("Wais 01BH vs IR correlation");
  //  corrN->Draw("lpSame");
  cout << "Wais 01BH vs IR Max|Min: ";
  cout << TMath::MaxElement(corr->GetN(),corr->GetY()) << "|" << TMath::MinElement(corr->GetN(),corr->GetY()) << endl;


  c1->cd(2);
  corr2->Draw("alp");
  corr2->GetXaxis()->SetRangeUser(-6,1);
  corr2->SetTitle("Wais 02BH vs IR correlation");
  cout << "Wais 02BH vs IR Max|Min: ";
  cout << TMath::MaxElement(corr2->GetN(),corr2->GetY()) << "|" ;
  cout << TMath::MinElement(corr2->GetN(),corr2->GetY()) << endl;

  c1->cd(3);
  corrWais->Draw("alp");
  corrWais->GetXaxis()->SetRangeUser(16,28);
  corrWais->SetTitle("Wais vs Wais corr");
  cout << "Wais vs Wais Max|Min: ";
  cout << TMath::MaxElement(corrWais->GetN(),corrWais->GetY()) << "|" ;
  cout << TMath::MinElement(corrWais->GetN(),corrWais->GetY()) << endl;
  

  c1->cd(4);
  iRAutoCorr->Draw("alp");
  iRAutoCorr->GetXaxis()->SetRangeUser(-4,4);
  iRAutoCorr->SetTitle("impulse response auto");

  c1->cd(5);
  waisAutoCorr->Draw("alp");
  waisAutoCorr->GetXaxis()->SetRangeUser(-4,4);
  waisAutoCorr->SetTitle("wais 01BH auto");

  c1->cd(6);
  waisAutoCorr2->Draw("alp");
  waisAutoCorr2->GetXaxis()->SetRangeUser(-4,4);
  waisAutoCorr2->SetTitle("wais 02BH auto");


  c1->cd(7);
  iRPadded->Draw("alp");
  iRPadded->SetTitle("impulse response");
  
  c1->cd(8);
  waisPulse->Draw("alp");

  c1->cd(9);
  waisPulse2->Draw("alp");


  TGraph *impulseResponses[96];
  stringstream fullName,name;
  for (int phi=0; phi<16; phi++) {
    for (int ring=0; ring<3; ring++) {
      for (int pol=0; pol<2; pol++) {
	int index = phi*6 + ring*2 + pol;
	name.str("");
	if (phi<9) name << "0";
	name << phi+1 << AnitaRing::ringAsChar((AnitaRing::AnitaRing_t)ring) << AnitaPol::polAsChar( (AnitaPol::AnitaPol_t) pol) << ".imp";
	cout << index << " " <<  name.str() << endl;
	fullName.str("");
	fullName << "/Users/brotter/anita16/local/share/UCorrelator/responses/IndividualBRotter/" << name.str();
	impulseResponses[index] = new TGraph(fullName.str().c_str());
      }
    }
  }

  
  TGraph *irCorrsVsEachOther[96];
  for (int i=0; i<96; i++) {
    irCorrsVsEachOther[i] = FFTtools::getCorrelationGraph(impulseResponses[0],impulseResponses[i]);
    cout << "max|min" << TMath::MaxElement(irCorrsVsEachOther[i]->GetN(),irCorrsVsEachOther[i]->GetY()) <<  "|" << TMath::MinElement(irCorrsVsEachOther[i]->GetN(),irCorrsVsEachOther[i]->GetY()) << endl;
  }
  
  
  for (int i=0; i<96; i++) {
    delete impulseResponses[i];
    delete irCorrsVsEachOther[i];
  }
  

  return;

}


double getTemplateCorrVal (TGraph *waisPulse, TGraph *finalWave) {
  //just get the max absolute template correlation value
  //make sure both these inputs have the same number of points and same dTs!
 
  int lenTemplate = waisPulse->GetN();
  if (lenTemplate != finalWave->GetN()) {
    cout << "!!!Warning in getTemplateCorrVal: waisPulse->GetN() != finalWave->GetN()" << endl; }

  double* corr = FFTtools::getCorrelation(lenTemplate,waisPulse->GetY(),finalWave->GetY());

  double max = TMath::MaxElement(lenTemplate,corr);
  double min = TMath::MinElement(lenTemplate,corr);
  
  delete[] corr;

  return TMath::Max(max,TMath::Abs(min));
}
  

void compareWaisTemplateVsRealEvents() {
  
  /*
    Another thing that is interesting is seeing how real events compare to the template
   */


  //lets get the wais tempate
  TFile *waisFile = TFile::Open("/Users/brotter/anita16/benPrograms/waisPulses/waisImpulseResponse_wAngs.root");

  //and calculate the normalization (method found and discribed in testNormalization.C)
  TGraph *waisPulse = (TGraph*)waisFile->Get("wais01BH");
  int lenTemplate = waisPulse->GetN();
  cout << "Length: " << lenTemplate << endl;
  double waisSum = 0;
  for (int pt=0; pt<lenTemplate; pt++) waisSum += pow(waisPulse->GetY()[pt],2);
  for (int pt=0; pt<lenTemplate; pt++) waisPulse->GetY()[pt] /= TMath::Sqrt(waisSum / (lenTemplate/4));

  cout << "waisPulse dT=" << waisPulse->GetX()[1] - waisPulse->GetX()[0] << endl;

  
  //Now lets load an event file from the list of wais files I generated
  TChain *headTree = new TChain("headTree","headTree");
  headTree->Add("/Users/brotter/anita16/rootFilesLocal/headFileWais.root");
  RawAnitaHeader *head = NULL;
  headTree->SetBranchAddress("header",&head);

  TChain *eventTree = new TChain("eventTree","eventTree");
  eventTree->Add("/Users/brotter/anita16/rootFilesLocal/calEventFileWais.root");
  CalibratedAnitaEvent *event = NULL;
  eventTree->SetBranchAddress("event",&event);

  //output stuff
  TGraph *gCorrAbs = new TGraph();
  gCorrAbs->SetName("gCorrAbs");


  int lenEntries = eventTree->GetEntries();  
  cout << "lenEntries = " << lenEntries << endl;
  lenEntries = 10000;
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%100 == 0) cout << "entry " << entry << "/" << lenEntries << endl;
    headTree->GetEntry(entry);
    eventTree->GetEntry(entry);

    double absMaxSum = 0;

    for (int ringi=0; (AnitaRing::AnitaRing_t)ringi != AnitaRing::kNotARing; ringi++) {
      AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t) ringi;
    
      //get waveform for that channel
      UsefulAnitaEvent *useful = new UsefulAnitaEvent(event);
      TGraph *waveform = useful->getGraph(ring,1,AnitaPol::kHorizontal);
      delete useful;
      //evenly sample
      TGraph *evenWave = FFTtools::getInterpolatedGraph(waveform,0.1);
      delete waveform;
      //zero pad so it is equal in length
      TGraph *finalWave = FFTtools::padWaveToLength(evenWave,lenTemplate);
      delete evenWave;
      
      //normalize it
      double waveSum = 0;
      for (int pt=0; pt<finalWave->GetN(); pt++) waveSum += pow(finalWave->GetY()[pt],2);
      for (int pt=0; pt<finalWave->GetN(); pt++) finalWave->GetY()[pt] /= TMath::Sqrt(waveSum / (finalWave->GetN()/4));

      absMaxSum += getTemplateCorrVal(waisPulse,finalWave);
      delete finalWave;
    }

    gCorrAbs->SetPoint(entry,entry,absMaxSum/3.);

  }

    
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gCorrAbs->Draw("alp");
  /*
  gCorrPeaks->Draw("ap");
  gCorrPeaks->GetYaxis()->SetRangeUser(-1,1);
  gCorrLows->SetMarkerColor(kRed);
  gCorrLows->Draw("pSame");
  */

  return;
}
