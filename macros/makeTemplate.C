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


  return;

}




TGraph *normalizeWaveform(TGraph *inGraph) {
  // from templateSearch.cc!!!!  I need a library for this :P
  TGraph *outGraph = (TGraph*)inGraph->Clone();
  
  //normalize it ( as seen in macros/testTemplate.C )
  double waveSum = 0;
  for (int pt=0; pt<outGraph->GetN(); pt++) waveSum += pow(outGraph->GetY()[pt],2);
  for (int pt=0; pt<outGraph->GetN(); pt++) outGraph->GetY()[pt] /= TMath::Sqrt(waveSum / (outGraph->GetN()/4));

  return outGraph;

}


void compTemplatesVsMeasured() {

  stringstream name;
  
  
  int length = 10000; //length of ZHAires waveforms

  TFile *measFile = TFile::Open("measuredCR.root");
  TGraph *measWaves[4];
  for (int i=0; i<4; i++) {
    name.str("");
    name << "Wave" << i << "pol0";
    TGraph* tempWave = (TGraph*)measFile->Get(name.str().c_str());
    TGraph* paddedWave = FFTtools::padWaveToLength(tempWave,length);
    measWaves[i] = normalizeWaveform(paddedWave);
    delete paddedWave;
  }


  TFile *tempFile = TFile::Open("/Users/brotter/benCode/ZHAiresReader/convolveCRWithSigChain.root");
  TGraph *templates[10];
  for (int i=0; i<10; i++) {
    int wave = i+13; //peak seems to be at around the 13th one, then by 23 it is basically zero
    name.str("");
    name << "wave" << wave;
    TGraph* tempWave = (TGraph*)tempFile->Get(name.str().c_str());
    templates[i] = normalizeWaveform(tempWave);
    
  }


  TH2D *hComp = new TH2D("hComp","correlation values",12,-0.5,11.5, 4,0,3.5);

  for (int measi=0; measi<4; measi++) {
    cout << "Measurement " << measi << " | ";
    double peaks[10];
    for (int tempi=0; tempi<10; tempi++) {
      double *corr = FFTtools::getCorrelation(length,measWaves[measi]->GetY(),templates[tempi]->GetY());
      double max = TMath::MaxElement(length,corr);
      double min = TMath::MinElement(length,corr);
      double peak = TMath::Max(max,-1.*min);
      peaks[tempi] = peak;
      cout << peak<< " ";
      hComp->Fill(tempi,measi,peak);
    }
    int peakLoc = TMath::LocMax(10,peaks);
    double peak = TMath::MaxElement(10,peaks);
    cout << " | max: " << peak << " @ " << peakLoc << endl;

  }

  TGraph* gIRraw = new TGraph("~/anita16/local/share/AnitaAnalysisFramework/responses/SingleBRotter/all.imp");
  TGraph *gIR = FFTtools::padWaveToLength(gIRraw,length);
  gIR->SetMarkerColor(kRed);
  gIR->SetLineColor(kRed);
  gIR->SetName("gIR");
  gIR->SetTitle("Impulse Response");


  cout << "Vs Impulse Response | ";
  for (int measi=0; measi<4; measi++) {
    double peaks[10];
    double *corr = FFTtools::getCorrelation(length,measWaves[measi]->GetY(),gIR->GetY());
    double max = TMath::MaxElement(length,corr);
    double min = TMath::MinElement(length,corr);
    double peak = TMath::Max(max,-1.*min);
    hComp->Fill(11,measi,peak);
    cout << peak<< " ";
  }
  cout << endl;

  hComp->Draw("colz");

  return;

}



void measuredVsEachOther() {

  stringstream name;

  const int length = 2048; //length of saved measured waveforms

  TFile *measFile = TFile::Open("measuredCR.root");
  TGraph *measWaves[4];
  for (int i=0; i<4; i++) {
    name.str("");
    name << "Wave" << i << "pol0";
    TGraph* tempWave = (TGraph*)measFile->Get(name.str().c_str());
    TGraph* paddedWave = FFTtools::padWaveToLength(tempWave,length);
    measWaves[i] = normalizeWaveform(paddedWave);
    delete paddedWave;
  }


  for (int measi=0; measi<4; measi++) {
    cout << "Measurement " << measi << " | ";
    for (int measj=0; measj<4; measj++) {
      double *corr = FFTtools::getCorrelation(length,measWaves[measi]->GetY(),measWaves[measj]->GetY());
      double max = TMath::MaxElement(length,corr);
      double min = TMath::MinElement(length,corr);
      double peak = TMath::Max(max,-1.*min);
      cout << peak<< " ";
    }
    cout << endl;
  }

  return;

}



void chansVsEachOther() {

  int length = 4096;

  char* baseDir = getenv("ANITA_UTIL_INSTALL_DIR");

  stringstream name;
  TGraph *grChans[96];
  for (int phi=0; phi<16; phi++) {
    for (int ringi=0; ringi<3; ringi++) {
      for (int poli=0; poli<2; poli++) {
	//	int chanIndex = phi*6 + ringi*2 + poli;
	int chanIndex = poli*48 + ringi*16 + phi;

	AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t)ringi;
	AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)poli;
	name.str("");
	name << baseDir << "/share/UCorrelator/responses/IndividualBRotter/";
	if (phi<9) name << "0";
	name << phi+1 << AnitaRing::ringAsChar(ring) << AnitaPol::polAsChar(pol) << ".imp";
	cout << chanIndex << " " << name.str() << endl;
	TGraph* tempGraph = new TGraph(name.str().c_str());
	TGraph *paddedGraph = FFTtools::padWaveToLength(tempGraph,4096);
	grChans[chanIndex] = normalizeWaveform(paddedGraph);
      }
    }
  }

  TH1D *corrs = new TH1D("hCorrs","hCorrs",100,0,1);

  TH2D *h2Corrs = new TH2D("h2Corrs","h2Corrs",96,-0.5,95.5,96,-0.5,95.5);

  for (int i=0; i<96; i++) {
    for (int j=0; j<96; j++) {
      double *corr = FFTtools::getCorrelation(length,grChans[i]->GetY(),grChans[j]->GetY());
	double max = TMath::MaxElement(length,corr);
	double min = TMath::MinElement(length,corr);
	double peak = TMath::Max(max,-1.*min);	
	h2Corrs->Fill(i,j,peak);

	if (j > i) corrs->Fill(peak);


      }
  }


  h2Corrs->Draw("colz");

}
      
