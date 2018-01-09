/*

  using Correlator.cc and setting a bunch of disallowed antennas to make example maps

 */

#include "AnitaEventSummary.h"
#include "UCFilters.h"

TH2* interferometry(int eventNumber=15717147) {

  AnitaDataset *data = new AnitaDataset(342);
  data->getEvent(eventNumber);

  FilterStrategy *fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");

  FilteredAnitaEvent *filtered = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());


  UCorrelator::Correlator *corr = new UCorrelator::Correlator(360,0,360,71,-50,20);

  corr->compute(filtered,AnitaPol::kHorizontal);
  TH2D *histTot = new TH2D(*corr->getHist());
  
  int x,y,z;
  int maxBin = histTot->GetMaximumBin();
  histTot->GetBinXYZ(maxBin,x,y,z);
  cout << x << endl;
  TH2 *histRot = UCorrelator::rotateHistogram(histTot,x+180);
  histRot->SetName("rotatedMap");
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  histRot->Draw("colz");

  

  return histRot;

}
  

void interferometry_example(int eventNumber=15717147) {

  AnitaDataset *data = new AnitaDataset(342);
  data->getEvent(eventNumber);

  FilterStrategy *fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");

  FilteredAnitaEvent *filtered = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());


  UCorrelator::Correlator *corr = new UCorrelator::Correlator(140,40,180,41,-20,20);


  const UCorrelator::AntennaPositions * ap = UCorrelator::AntennaPositions::instance(); 

  int closest[15];
  int nant = ap->getClosestAntennas(110, 15, closest, 0, AnitaPol::kHorizontal);

  for (int i=0; i<15; i++) {
    cout << closest[i] << endl; }


  ULong64_t allowed = 0;  
  TH2D* hists[6];
  //set two of the antennas to 1 iteratively,
  //0,1 0,2 0,3 1,2 1,3 2,3 
  int cnt=0;
  for (int i=0; i<4; i++) {
    for (int j=i+1; j<4; j++) {
      cout << closest[i] << "," << closest[j] << endl;
      allowed = 0;
      allowed |= (1ul << closest[i]);
      allowed |= (1ul << closest[j]);
      corr->setAllowedAntennas(allowed);
      corr->compute(filtered,AnitaPol::kHorizontal);
      hists[cnt] = new TH2D(*corr->getHist());
      cnt++;
    }
  }
  allowed=0;
  for (int i=0; i<4; i++) {
    allowed |= (1ul << closest[i]);
  }
  cout << allowed << endl;
  corr->setAllowedAntennas(allowed);
  corr->compute(filtered,AnitaPol::kHorizontal);
  TH2D *histTot = new TH2D(*corr->getHist());
    
  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.0, 0.8,1.0);
  pad1->Divide(3,2);
  TPad *pad2 = new TPad("pad2", "The pad 80% of the height",0.8,0.0, 1.0,1.0);
  pad1->Draw();
  pad2->Draw();
  for (int i=0; i<6; i++) {
    pad1->cd(i+1);
    hists[i]->Draw("colz");
  }
  pad2->cd();
  histTot->Draw("colz");

  

}
  
