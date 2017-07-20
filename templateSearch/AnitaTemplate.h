#ifndef _ANITA_TEMPLATE_H_
#define _ANITA_TEMPLATE_H_

#include "TObject.h" 
#include "AnitaConventions.h" 
#include "Analyzer.h"
#include "FFTtools.h"
#include "TFile.h"

/** Common analysis output format for template correlations
 *  
 *  Needless to say, there's no guarantee that everything will be filled, so be wary if something is 0 (it may not have been filled).   
 *
 *  Shamelessly copied format from AnitaEventSummary
 * */ 


class AnitaTemplate 
{
 public:

  static const int templateLength = 2048;


  class TemplateResults
  {
  public:
    TemplateResults () { ; }
    double templateValue[NUM_POLS];
    double templatePeakLoc[NUM_POLS];

    ClassDefNV(TemplateResults,1); 
  };


  FFTWComplex *waisTemplate;
  FFTWComplex *impulseResponseTemplate;
  FFTWComplex *measuredCRTemplate;
  FFTWComplex *generatedCRTemplate;

  TemplateResults *waisResults;
  TemplateResults *impulseResponseResults;
  TemplateResults *measuredCRResults;
  TemplateResults *generatedCRResults;

  void fillTemplates();
  void getImpulseResponseTemplate();  
  void getWaisTemplate();  

  

  void fillTemplateResults(UCorrelator::Analyzer *analyzer);

  TGraph *normalizeWaveform(TGraph *inGraph);

  double calcTemplateValue(FFTWComplex* FFT1, FFTWComplex* FFT2);
  double* getCorrelationFromFFT(int length, const FFTWComplex *FFT1, const FFTWComplex *FFT2);

  void zeroInternals();

  AnitaTemplate();


 private:
  ClassDefNV(AnitaTemplate, 0); 

};





#endif
  
