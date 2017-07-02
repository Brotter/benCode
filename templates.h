#ifndef ANITA_TEMPLATES
#define ANITA_TEMPLATES

#include "AnitaEventSummary.h"
#include "AnitaConventions.h"

TGraph *normalizeWaveform(TGraph *inGraph) {
  
  TGraph *outGraph = (TGraph*)inGraph->Clone();
  
  //normalize it ( as seen in macros/testTemplate.C )
  double waveSum = 0;
  for (int pt=0; pt<outGraph->GetN(); pt++) waveSum += pow(outGraph->GetY()[pt],2);
  for (int pt=0; pt<outGraph->GetN(); pt++) outGraph->GetY()[pt] /= TMath::Sqrt(waveSum / (outGraph->GetN()/4));

  return outGraph;

}




double *getCorrelationFromFFT(int length,const FFTWComplex *theFFT1, const FFTWComplex *theFFT2) 
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



void getImpulseResponseTemplate(int length, FFTWComplex *theTemplateFFT, TGraph *theTemplate) {
  
  //and get the "averaged" impulse response as the template"
  char* templateDir = getenv("ANITA_UTIL_INSTALL_DIR");
  stringstream name;
  name.str("");
  name << templateDir << "/share/AnitaAnalysisFramework/responses/SingleBRotter/all.imp";
  TGraph *grTemplateRaw = new TGraph(name.str().c_str());
  //waveforms are normally a little over 1024 so lets pad to 2048 (defined in length above)
  TGraph *grTemplatePadded = FFTtools::padWaveToLength(grTemplateRaw,length);
  delete grTemplateRaw;
  //then cut it back down with a window function
  //  int peakHilb = -1;
  //  TGraph *grTemplateCut = windowDispersed(grTemplatePadded,peakHilb);
  //  delete grTemplatePadded;
  //and finally normalize it (last step!)
  TGraph *grTemplate = normalizeWaveform(grTemplatePadded);
  delete grTemplatePadded;

  //give it a name
  grTemplate->SetName("templateImp");

  theTemplate = grTemplate;

  cout << "Impulse Response Template Length: " << grTemplate->GetN() << endl;

  //and get the FFT of it as well, since we don't want to do this every single event
  theTemplateFFT=FFTtools::doFFT(grTemplate->GetN(),grTemplate->GetY());
  
  return;
}

void getCRTemplates(int length, const int numTemplates,FFTWComplex** theTemplateFFTs, TGraph** theTemplates) {

  stringstream name;


  //  TFile *inFile = TFile::Open("/Users/brotter/benCode/ZHAiresReader/convolveCRWithSigChain.root");
  TFile *inFile = TFile::Open("convolveCRWithSigChain.root");
  
  for (int i=0; i<numTemplates; i++) {
    //want to get graphs 13 through 24 (like in makeTemplate.C)
    int wave = i+13; //peak seems to be at around the 13th one, then by 23 it is basically zero
    name.str("");
    name << "dispersedCR" << wave;
    TGraph *grTemplateRaw = (TGraph*)inFile->Get(name.str().c_str());
    
    //waveforms are super long so we can just cut it to the window dimentions
    int peakHilb = -1;
    TGraph *grTemplateCut = windowCut(grTemplateRaw,length);
    delete grTemplateRaw;
    
    
    //and finally normalize it (last step!)
    TGraph *grTemplate = normalizeWaveform(grTemplateCut);
    delete grTemplateCut;
    
    //give it a name
    grTemplate->SetName(name.str().c_str());
    
    cout << "CR Template " << i << " Length: " << grTemplate->GetN() << endl;
    
    //and get the FFT of it as well, since we don't want to do this every single event
    FFTWComplex *theTemplateFFT=FFTtools::doFFT(grTemplate->GetN(),grTemplate->GetY());
    delete grTemplate;
    theTemplates[i] = grTemplate;
      theTemplateFFTs[i] = theTemplateFFT;
  }
  
  inFile->Close();

  return;
}

    
void getWaisTemplate(int length, FFTWComplex *theTemplateFFT, TGraph* theTemplate) {
  
  //and get the "averaged" impulse response as the template"
  TFile *inFile = TFile::Open("waisTemplate.root");
  TGraph *grTemplateRaw = (TGraph*)inFile->Get("wais01TH");
  //the wais waveform is like N=2832, but most of it is dumb, so cut off the beginning
  //actually just window it!
  int peakHilb = -1;
  TGraph *grTemplateCut = windowCut(grTemplateRaw,length);
  delete grTemplateRaw;

  //and then normalize it
  TGraph *grTemplate = normalizeWaveform(grTemplateCut);
  delete grTemplateCut;

  //give it a name
  grTemplate->SetName("templateWais");

  cout << "Wais Template Length: " << grTemplate->GetN() << endl;

  //and get the FFT of it as well, since we don't want to do this every single event
  theTemplateFFT=FFTtools::doFFT(grTemplate->GetN(),grTemplate->GetY());
  theTemplate = grTemplate;
  
  return;
}


class AnitaTemplates
{
 public:

  /** Constructor.  Loads all the templates **/
  AnitaTemplates();

  /** Destructor. **/
  virtual ~AnitaTemplates();

  /* Length of the templates */
  static const int length = 2048;

  /* Impulse Response Template */
  TGraph *theImpTemplate;
  FFTWComplex *theImpTemplateFFT;

  /* Impulse Response Template Loader */
  void getImpulseResponseTemplate(FFTWComplex *theTemplateFFT, TGraph *grTemplate);


  /* Wais Template */
  TGraph *theWaisTemplate;
  FFTWComplex *theWaisTemplateFFT;

  /* Wais Template Loader */
  void getWaisTemplate(FFTWComplex *theTemplateFFT, TGraph *grTemplate);


  /* Cosmic Ray Templates */
  /* Number of points on the cone */
  static const int numCRTemplates = 10;

  /* Templates */
  TGraph *theCRTemplates[numCRTemplates];
  FFTWComplex *theCRTemplateFFTss[numCRTemplates];
  
  /* CR Template Loader */
  void getCRTemplates(const int numTemplates, FFTWComplex** theTemplateFFTs, TGraph** grTemplates);

 private:

  ClassDefNV(AnitaTemplates,1);

};



class AnitaTemplateResults
{
 public:

  
  /** The maximum number of hypotheses storable per polarization */ 
  //
  static const Int_t maxDirectionsPerPol = AnitaEventSummary::maxDirectionsPerPol; 

  /*The template correlatin values for a single coherent waveform comparison*/
  class SingleTemplateResult {
  public:
    SingleTemplateResult() {; }
    //impulse response
    Double_t templateImp;
    
    //one for the WAIS template too
    Double_t templateWais;
    
    //and for the bigger multi-coherence-angle one
    Double_t templateCRay[AnitaTemplates::numCRTemplates];
    
    ClassDefNV(SingleTemplateResult,1);
  };

  
  SingleTemplateResult coherentV[maxDirectionsPerPol];
  SingleTemplateResult coherentH[maxDirectionsPerPol];
  
  SingleTemplateResult deconvolvedV[maxDirectionsPerPol];
  SingleTemplateResult deconvolvedH[maxDirectionsPerPol];

    
 private:
  ClassDefNV(AnitaTemplateResults, 1); 
};
 
 
 

#endif
 
