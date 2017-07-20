#include "FFTtools.h"

//just making sure this code in templateSearch.cc worked.  Turns out the bug was in a V/H mismatch typo


TGraph *normalizeWaveform(TGraph *inGraph) {
  
  TGraph *outGraph = (TGraph*)inGraph->Clone();
  
  //normalize it ( as seen in macros/testTemplate.C )
  double waveSum = 0;
  for (int pt=0; pt<outGraph->GetN(); pt++) waveSum += pow(outGraph->GetY()[pt],2);
  for (int pt=0; pt<outGraph->GetN(); pt++) outGraph->GetY()[pt] /= TMath::Sqrt(waveSum / (outGraph->GetN()/4));

  return outGraph;

}



void waisTemplateTest() {

  int length = 2048;

  //and get the "averaged" impulse response as the template"                                                             
  TFile *inFile = TFile::Open("/Users/brotter/anita16/benPrograms/waisPulses/waisImpulseResponse_wAngs.root");
  TGraph *grTemplateRaw = (TGraph*)inFile->Get("wais01TH");
  //the wais waveform is like N=2832, but most of it is dumb, so cut off the beginning                                   
  TGraph *grTemplateCut = new TGraph();
  for (int pt=0; pt<grTemplateRaw->GetN(); pt++) {
    if (grTemplateRaw->GetX()[pt] > 0) grTemplateCut->SetPoint(grTemplateCut->GetN(),grTemplateRaw->GetX()[pt],grTemplateRaw->GetY()[pt]);
  }
  inFile->Close();
  //of course this way of doing it probably makes it too short :P                                                        
  TGraph *grTemplatePadded = FFTtools::padWaveToLength(grTemplateCut,length);
  delete grTemplateCut;
  //and then normalize it                                                                                                
  TGraph *grTemplate = normalizeWaveform(grTemplatePadded);
  delete grTemplatePadded;

  //and get the FFT of it as well, since we don't want to do this every single event                                     
  grTemplate->Draw("alp");

}


