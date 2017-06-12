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
    double templateValueH;
    double templatePeakLocH;
    double templateValueV;
    double templatePeakLocV;

    ClassDefNV(TemplateResults,12); 
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
  FFTWComplex* getImpulseResponseTemplate();  
  FFTWComplex* getWaisTemplate();  

  
  void fillTemplateResults(TemplateResults *theResults, FFTWComplex *theTemplateFFT, UCorrelator::Analyzer *analyzer);
  TGraph *normalizeWaveform(TGraph *inGraph);


  double* getCorrelationFromFFT(int length, const FFTWComplex *FFT1, const FFTWComplex *FFT2);

  void zeroInternals();

  AnitaTemplate();


 private:
  ClassDefNV(AnitaTemplate, 15); 

};





#endif
  
