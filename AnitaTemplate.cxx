#include "AnitaTemplate.h"


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
AnitaTemplate::AnitaTemplate(){
  zeroInternals();
  fillTemplates();
}


void AnitaTemplate::zeroInternals(){

  
  memset(waisResults,0,sizeof(TemplateResults)); 
  memset(impulseResponseResults,0,sizeof(TemplateResults)); 
  memset(measuredCRResults,0,sizeof(TemplateResults)); 
  memset(generatedCRResults,0,sizeof(TemplateResults)); 


  return;
}


void AnitaTemplate::fillTemplates(){
  waisTemplate = getWaisTemplate();
  impulseResponseTemplate = getImpulseResponseTemplate();

  return;
}


void AnitaTemplate::fillTemplateResults(TemplateResults *theResults, FFTWComplex *theTemplateFFT, UCorrelator::Analyzer *analyzer) {
  for (int poli=0; poli<2; poli++) {
    const TGraphAligned *coherentAligned = analyzer->getCoherent((AnitaPol::AnitaPol_t)poli,0)->even();
    TGraph *coherent = new TGraph(coherentAligned->GetN(),coherentAligned->GetX(),coherentAligned->GetY());
    //make sure it is the same length as the template
    TGraph *coherent2 = FFTtools::padWaveToLength(coherent,AnitaTemplate::templateLength);
    delete coherent;
    //normalize it
    TGraph *normCoherent = normalizeWaveform(coherent2);
    delete coherent2;
    FFTWComplex *coherentFFT=FFTtools::doFFT(AnitaTemplate::templateLength,normCoherent->GetY());
    delete normCoherent;
    
    double *dCorr = getCorrelationFromFFT(AnitaTemplate::templateLength,theTemplateFFT,coherentFFT);
    double max = TMath::MaxElement(AnitaTemplate::templateLength,dCorr);
    double min = TMath::Abs(TMath::MinElement(AnitaTemplate::templateLength,dCorr));
    
    if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kHorizontal ) {
      theResults->templateValueH = TMath::Max(max,min);
    }
    if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kVertical ) {
      theResults->templateValueV = TMath::Max(max,min);
    }
    
    delete[] coherentFFT;
    delete[] dCorr;
  }
  return;
}



FFTWComplex* AnitaTemplate::getImpulseResponseTemplate() {
  
  //and get the "averaged" impulse response as the template"
  char* templateDir = getenv("ANITA_UTIL_INSTALL_DIR");
  std::stringstream name;
  name.str("");
  name << templateDir << "/share/AnitaAnalysisFramework/responses/SingleBRotter/all.imp";
  TGraph *grTemplateRaw = new TGraph(name.str().c_str());
  //waveforms are normally a little over 1024 so lets pad to 2048 (defined in length above)
  TGraph *grTemplatePadded = FFTtools::padWaveToLength(grTemplateRaw,AnitaTemplate::templateLength);
  delete grTemplateRaw;
  //and normalize it
  TGraph *grTemplate = normalizeWaveform(grTemplatePadded);
  delete grTemplatePadded;
  
  //and get the FFT of it as well, since we don't want to do this every single event
  FFTWComplex *theTemplateFFT=FFTtools::doFFT(AnitaTemplate::templateLength,grTemplate->GetY());
  delete grTemplate;
  
  return theTemplateFFT;
}



FFTWComplex* AnitaTemplate::getWaisTemplate() {
    
  //and get the "averaged" impulse response as the template"
  TFile *inFile = TFile::Open("waisTemplate.root");
  TGraph *grTemplateRaw = (TGraph*)inFile->Get("wais01TH");
  //the wais waveform is like N=2832, but most of it is dumb, so cut off the beginning
  TGraph *grTemplateCut = new TGraph();
  for (int pt=0; pt<grTemplateRaw->GetN(); pt++) {
    if (grTemplateRaw->GetX()[pt] > 0) grTemplateCut->SetPoint(grTemplateCut->GetN(),grTemplateRaw->GetX()[pt],grTemplateRaw->GetY()[pt]);
  }
  inFile->Close();
  //of course this way of doing it probably makes it too short :P
  TGraph *grTemplatePadded = FFTtools::padWaveToLength(grTemplateCut,AnitaTemplate::templateLength);
  delete grTemplateCut;
  //and then normalize it
  TGraph *grTemplate = normalizeWaveform(grTemplatePadded);
  delete grTemplatePadded;
  
  //and get the FFT of it as well, since we don't want to do this every single event
  FFTWComplex *theTemplateFFT=FFTtools::doFFT(AnitaTemplate::templateLength,grTemplate->GetY());
  delete grTemplate;
  
  return theTemplateFFT;
}





TGraph *AnitaTemplate::normalizeWaveform(TGraph *inGraph) {
  
  TGraph *outGraph = (TGraph*)inGraph->Clone();
  
  //normalize it ( as seen in macros/testTemplate.C )
  double waveSum = 0;
  for (int pt=0; pt<outGraph->GetN(); pt++) waveSum += pow(outGraph->GetY()[pt],2);
  for (int pt=0; pt<outGraph->GetN(); pt++) outGraph->GetY()[pt] /= TMath::Sqrt(waveSum / (outGraph->GetN()/4));

  return outGraph;

}




double *AnitaTemplate::getCorrelationFromFFT(int length,const FFTWComplex *theFFT1, const FFTWComplex *theFFT2) 
{


    int newLength=(length/2)+1;
//     cout << "newLength " << newLength << endl;
    FFTWComplex *tempStep = new FFTWComplex [newLength];
    int no2=length>>1;
    for(int i=0;i<newLength;i++) {
	double reFFT1=theFFT1[i].re;
	double imFFT1=theFFT1[i].im;
	double reFFT2=theFFT2[i].re;
	double imFFT2=theFFT2[i].im;

	//Real part of output 
	tempStep[i].re=(reFFT1*reFFT2+imFFT1*imFFT2)/double(no2/2);
	//Imaginary part of output 
	tempStep[i].im=(imFFT1*reFFT2-reFFT1*imFFT2)/double(no2/2);
    }
//    cout << "finished messing around" << endl;
    double *theOutput=FFTtools::doInvFFT(length,tempStep);
//    cout << "got inverse" << endl;
    delete [] tempStep;
    return theOutput;

}
