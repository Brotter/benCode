#include "Analyzer.h"
#include "AnalysisConfig.h"
#include "SpectrumAverage.h"

/*

  Using the impulse response and different off angles from the ZHAires simulation, make some templates!

  Uses the results from ZHAiresReader and convolves it with the stuff from transferFunctionFinal

 */




void makeCalculatedTemplate() {

  cout << "This one doesn't do anything yet" << endl;

  return;

}



void makeMeasuredTemplate() {

  /*

    Peter seems to trust the measured events most, so this is a coherent sum of BenS's events that he found

    ev#s = 32907848, 33484995, 41529195, 58592863
    from runs 248, 250, 284, 343

    The plan here is to use UCorrelator to generate the coherent sums, then coherently sum those together and use it as
    a template

   */

  
  AnitaDataset *data[4];
  
  data[0] = new AnitaDataset(248,false);
  data[0]->getEvent(32907848);
  data[1] = new AnitaDataset(250,false);
  data[1]->getEvent(33484995);
  data[2] = new AnitaDataset(284,false);
  data[2]->getEvent(41529195);
  data[3] = new AnitaDataset(343,false);
  data[3]->getEvent(58592863);

  //Make a filter strategy
  FilterStrategy *strategy = new FilterStrategy();


  //the spectrum average is used for a couple of filters to make them sort of "adaptive"
  char* specAvgDir = getenv("UCORRELATOR_SPECAVG_DIR");
  const UCorrelator::SpectrumAverageLoader *specAvgLoader = new UCorrelator::SpectrumAverageLoader(specAvgDir);

  //Sine subtract alghorithm (this is the complicated way to do it)
  //  UCorrelator::SineSubtractFilter *sineSub = new UCorrelator::SineSubtractFilter(0.05,2);
  //  sineSub->makeAdaptive(specAvgLoader);
  //  strategy->addOperation(sineSub);
  // This seems like it should work and is easier
  //add "adsinsub_2_5_13" (default in MagicDisplay)
  //  UCorrelator::fillStrategyWithKey(strategy,"sinsub_05_1_ad_1");
  

  //Brick wall filter, should be way faster
  // UCorrelator::AdaptiveBrickWallFilter(const UCorrelator::SpectrumAverageLoader * spec, double thresh=2, bool fillNotch = true);  
  // Don't fill in the noise because whats the point of that really
  UCorrelator::AdaptiveBrickWallFilter *brickWall = new UCorrelator::AdaptiveBrickWallFilter(specAvgLoader,2,false);
  strategy->addOperation(brickWall);


  //  with abby's list of filtering
  //  UCorrelator::applyAbbysFilterStrategy(&strategy);



  //and a configuration for the analysis
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig(); 
  //set the response to my "single" response
  //  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseIndividualBRotter;
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;

  //and create an analyzer object
  UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true); ;

  AnitaEventSummary *eventSummary = new AnitaEventSummary; 

  stringstream name;

  TGraph *coherent[4][2];
  for (int i=0; i<4; i++) {
    FilteredAnitaEvent *filtEv = new FilteredAnitaEvent(data[i]->useful(), strategy, data[i]->gps(), data[i]->header());
    analyzer->analyze(filtEv, eventSummary); 
    for (int poli=0; poli<2; poli++) {
      cout << "wave " << i << " pol " << poli << endl;
      const TGraphAligned *coherentAligned = analyzer->getCoherent((AnitaPol::AnitaPol_t)poli,0)->even();
      TGraph *coherentRaw = new TGraph(coherentAligned->GetN(),coherentAligned->GetX(),coherentAligned->GetY());
      //make sure it is the same length as the template
      
      coherent[i][poli] = FFTtools::padWaveToLength(coherentRaw,2048);
      name.str("");
      name << data[i]->header()->eventNumber;
      if (poli) name << "H";
      else name << "V";
      coherent[i][poli]->SetTitle(name.str().c_str());
      name.str("");
      name << "Wave" << i << "pol" << poli;
      coherent[i][poli]->SetName(name.str().c_str());
      delete coherentRaw;
    }
  }

  
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->Divide(2,4);
  c1->cd(1);
  coherent[0][0]->Draw("alp");
  c1->cd(2);
  coherent[0][1]->Draw("alp");
  c1->cd(3);
  coherent[1][0]->Draw("alp");
  c1->cd(4);
  coherent[1][1]->Draw("alp");
  c1->cd(5);
  coherent[2][0]->Draw("alp");
  c1->cd(6);
  coherent[2][1]->Draw("alp");
  c1->cd(7);
  coherent[3][0]->Draw("alp");
  c1->cd(8);
  coherent[3][1]->Draw("alp");

  TFile *outFile = TFile::Open("measuredCR.root","recreate");
  for (int i=0; i<4; i++) {
    for (int j=0; j<2; j++) {
      coherent[i][j]->Write();
    }
  }
  outFile->Close();

  return;
  }


void drawMeasured() {

  TFile *inFile = TFile::Open("measuredCR.root");
  TGraph* gIRraw = new TGraph("~/anita16/local/share/AnitaAnalysisFramework/responses/SingleBRotter/all.imp");
  TGraph *gIR = FFTtools::padWaveToLength(gIRraw,2048);
  for (int pt=0; pt<gIR->GetN(); pt++) gIR->GetY()[pt] *= 6;
  gIR->SetMarkerColor(kRed);
  gIR->SetLineColor(kRed);
  gIR->SetName("gIR");
  gIR->SetTitle("Impulse Response");

  stringstream name;

  TGraph *waves[5];
  for (int i=0; i<4; i++) {
    name.str("");
    name << "Wave" << i << "pol0";
    waves[i] = (TGraph*)inFile->Get(name.str().c_str());

    double yMax = waves[i]->GetX()[TMath::LocMax(waves[i]->GetN(),waves[i]->GetY())];
    for (int pt=0; pt<waves[i]->GetN(); pt++) {
      double xVal = waves[i]->GetX()[pt];

      //      if (xVal < yMax-10 || xVal > yMax+50) waves[i]->GetY()[pt] = 0;
    }
  }
  waves[4] = gIR;
  
  int shifts[4];
  for (int i=0; i<4; i++) {
    double *corr = FFTtools::getCorrelation(waves[0]->GetN(),waves[0]->GetY(),waves[i+1]->GetY());
    int maxLoc = TMath::LocMax(waves[0]->GetN(),corr);
    int minLoc = TMath::LocMin(waves[0]->GetN(),corr);
    int max = TMath::MaxElement(waves[0]->GetN(),corr);
    int min = TMath::MinElement(waves[0]->GetN(),corr);

    cout << maxLoc << " " << max << " " << minLoc << " " << min << endl;
    if (max > -1*min) shifts[i] = maxLoc;
    else shifts[i] = minLoc;
    if (shifts[i] > waves[0]->GetN()/2.) shifts[i] -= waves[0]->GetN();
    delete[] corr;
    double dT = waves[i+1]->GetX()[1] - waves[i+1]->GetX()[0];
    cout << shifts[i] << " " << dT << endl;

    for (int pt=0; pt<waves[i+1]->GetN(); pt++) {
      waves[i+1]->GetX()[pt] += shifts[i]*dT;

    }
  }

  


  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->Divide(2,3);

  c1->cd(6);
  waves[0]->SetMarkerColor(kBlue);
  waves[0]->SetLineColor(kBlue);
  waves[0]->Draw("alp");
  //  waves[0]->GetXaxis()->SetRangeUser(0,50);
  for (int i=1; i<5; i++) {    
    waves[i]->Draw("lpSame");
  }
  for (int i=0; i<5; i++) {
    c1->cd(i+1);
    waves[i]->Draw("alp");
    waves[i]->GetXaxis()->SetRangeUser(0,50);
  }




}
