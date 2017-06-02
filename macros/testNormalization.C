/*

  I want to do a correlation of the template vs a bunch of stuff and see what the various values are.


  Like for example vs an averaged wais pulse, or vs various off angles from the coherence angle

  this script sort of just plays around with the idea and the normalization

  Note that there are SIGNIFICANT differences between FFTtools::getCorrelationGraph and FFTtools::getCorrelation
  chief amongst them is that getCorrelationGraph zero pads the waveforms so that they are much longer first, which
  will change the normalization a lot!  However, it also will sort of help with making them the same length

  I should probably just do zero padding and length matching by hand, because then it is clear.
  then use FFTtools::getCorrelation

 */


#include "FFTtools.h"
#include "AnitaConventions.h"
void testNormalization(){

  //get a waveform (this is from correlating and aligning wais pulses)
  TFile *waisFile = TFile::Open("/Users/brotter/anita16/benPrograms/waisPulses/waisImpulseResponse_wAngs.root");
  TGraph *waisPulse = (TGraph*)waisFile->Get("wais01BH");
  cout << "Wais Length: " << waisPulse->GetN() << endl;

  //calculate the sum of squares (should be a fast normalization constant)
  double waisSum = 0;
  for (int pt=0; pt<waisPulse->GetN(); pt++) waisSum += pow(waisPulse->GetY()[pt],2);
  cout << "Wais sum: " << waisSum << endl;

  //calculate the autocorrelation using getCorrelationGraph
  TGraph *waisAutoCorr = FFTtools::getCorrelationGraph(waisPulse,waisPulse);
  cout << "waisAutoCorr->GetN():" << waisAutoCorr->GetN() << endl;
  double waisMax = TMath::MaxElement(waisAutoCorr->GetN(),waisAutoCorr->GetY());

  //calculate the autocorrelation using getCorrelation
  cout << "waisAutoCorr Max:" << waisMax << endl;
  double* waisGetCorr = FFTtools::getCorrelation(waisPulse->GetN(),waisPulse->GetY(),waisPulse->GetY());
  cout << "waisGetCorr->GetN():" << waisPulse->GetN() << endl;
  double corrMax = TMath::MaxElement(waisPulse->GetN(),waisGetCorr);
  cout << "waisGetCorr Max:" << corrMax << endl;

  //determine how each of those is different from the sum method (hopefully something we can replicate!)
  cout << "sum/autoMax " << (float)waisSum/(float)waisMax << endl;
  cout << "sum/corrMax " << (float)waisSum/(float)corrMax << endl;

  //turns out the above are both N/4 (though they have different Ns) for whatever reason
  //so lets try normalizing them and seeing if it works


  //previously dividing the waveform by sqrt(autoMax) had the correct normalization
  // but then you have to calculate the autoCorrelation max for each waveform!  Lots of FFTs!
  //  This is why I wanted to try the sum of squares method, since it should get the same thing but faster
  TGraph *waisCorrGraphNorm1 = new TGraph();
  TGraph *waisCorrGraphNorm2 = new TGraph();
  TGraph *waisCorrNorm1 = new TGraph();
  TGraph *waisCorrNorm2 = new TGraph();
  for (int pt=0; pt<waisPulse->GetN(); pt++) {

    //this isn't quite going to work, because you'll have to calculate what you think the getCorrelationGraph()
    // routine is going to decide on for the final length.  This is not a great way to do it!
    waisCorrGraphNorm1->SetPoint(pt,waisPulse->GetX()[pt],waisPulse->GetY()[pt]/TMath::Sqrt(waisMax));
    waisCorrGraphNorm2->SetPoint(pt,waisPulse->GetX()[pt],waisPulse->GetY()[pt]/TMath::Sqrt(waisSum/(waisPulse->GetN()/4)));
    
    //this will work better, though getCorrelation() is sort of a more clumbsy things to have to use
    waisCorrNorm1->SetPoint(pt,waisPulse->GetX()[pt],waisPulse->GetY()[pt]/TMath::Sqrt(corrMax));
    waisCorrNorm2->SetPoint(pt,waisPulse->GetX()[pt],waisPulse->GetY()[pt]/TMath::Sqrt(waisSum/(waisPulse->GetN()/4)));
  }


    
  TGraph *autoCorrGraphNorm1 = FFTtools::getCorrelationGraph(waisCorrGraphNorm1,waisCorrGraphNorm1);
  TGraph *autoCorrGraphNorm2 = FFTtools::getCorrelationGraph(waisCorrGraphNorm2,waisCorrGraphNorm2);
  TGraph *autoCorrNorm1 = new TGraph(waisPulse->GetN(),waisPulse->GetX(), FFTtools::getCorrelation(waisCorrNorm1->GetN(),waisCorrNorm1->GetY(),waisCorrNorm1->GetY()));
  TGraph *autoCorrNorm2 = new TGraph(waisPulse->GetN(),waisPulse->GetX(),FFTtools::getCorrelation(waisCorrNorm2->GetN(),waisCorrNorm2->GetY(),waisCorrNorm2->GetY()));


				     
  TCanvas *c1 = new TCanvas("hi","hi",800,600);
  c1->Divide(2);
  c1->cd(1);
  autoCorrGraphNorm1->Draw("al");
  autoCorrGraphNorm2->SetMarkerColor(kRed);
  autoCorrGraphNorm2->Draw("pSame");
  c1->cd(2);
  autoCorrNorm1->Draw("al");
  autoCorrNorm2->SetMarkerColor(kRed);
  autoCorrNorm2->Draw("pSame");


  return;

}



