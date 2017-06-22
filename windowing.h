#include <iostream>
#include <string>
#include <sstream>
//root
#include "TGraph.h"
//anita
#include "FFTtools.h"


using namespace std;



/*====================
windowing functions
*/

TGraph *windowWave(TGraph*, int&, const int, const int, const int, const int);

TGraph *windowDispersed(TGraph *inGraph, int &peakHilbertLoc) {  
  return windowWave(inGraph,peakHilbertLoc,
		    50,250,20,20);
}

TGraph *windowEField(TGraph *inGraph, int &peakHilbertLoc) {
  return windowWave(inGraph,peakHilbertLoc,
		    30,30,2,2);
}


TGraph *windowWave(TGraph *inGraph, int &peakHilbertLoc,
		   const int upRampLoc = 50,   const int downRampLoc = 600,
		   const int upRampLen = 100, const int downRampLen = 400) {
  //defaults are for the impulse response

  bool debug = false;

  /*

    The noise after the waveform part is useless.  I need to window it to increase the correlation value

    Find the peak of the hilbert envelope, then go 5ns before it (50pts) and then do a hamming maybe like 60 after it?

   */


  //the following are window config params, in POINTS (not nanoseconds)
  // downRampLoc - how far after peak hilbert env to start tail hamming
  // upRampLoc - how far before peak hilbert to start hamming (well close)
  
  // upRampLen: how "long" the hamming should be (half period)
  // downRampLen - how "long" the hamming should be at the end (well close)

  
  //If I don't tell it where to window, it should figure it out
  if (peakHilbertLoc == -1) {
    TGraph *hilbert = FFTtools::getHilbertEnvelope(inGraph);
    peakHilbertLoc = TMath::LocMax(hilbert->GetN(),hilbert->GetY());
    delete hilbert; //no memory leaks!
  }

  int startLoc = peakHilbertLoc - upRampLoc;
  int stopLoc  = peakHilbertLoc + downRampLoc;

  if (stopLoc+downRampLen > inGraph->GetN()) {
    if (debug) cout << "****";
    int overrun = (stopLoc+downRampLen) - inGraph->GetN() + 1;
    startLoc -= overrun;
    stopLoc -= overrun;
  }

  if (debug) cout << "inGraph->GetN()=" << inGraph->GetN() << " startLoc=" << startLoc << " stopLoc=" << stopLoc;



  TGraph *outGraph = new TGraph();
  for (int pt=0; pt<inGraph->GetN(); pt++) {
    if (pt <= (startLoc-upRampLen)) {
      outGraph->SetPoint(outGraph->GetN(),inGraph->GetX()[pt],0);
      if (debug) cout << pt << " - 1" << endl;
    }
    else if (pt > (startLoc-upRampLen) && pt <= startLoc ) {
      int ptMod = pt - (startLoc-upRampLen);
      double modValue = 0.5-(TMath::Cos(ptMod * ( TMath::Pi()/upRampLen ))/2.);
      double value =  modValue * inGraph->GetY()[pt];
      outGraph->SetPoint(outGraph->GetN(),inGraph->GetX()[pt],value);
      if (debug) cout << pt << " - 2 - " << ptMod << " - " << modValue << endl;
    }
    else if (pt > startLoc && pt <= stopLoc) {
      outGraph->SetPoint(outGraph->GetN(),inGraph->GetX()[pt],inGraph->GetY()[pt]);
      if (debug) cout << pt << " - 3" << endl;
    }
    else if (pt > stopLoc && pt <= (stopLoc+downRampLen)) {
      double ptMod = pt - stopLoc;
      double modValue = (1+TMath::Cos(ptMod*( TMath::Pi()/downRampLen ))/2.) - 0.5;
      double value = modValue * inGraph->GetY()[pt];
      outGraph->SetPoint(outGraph->GetN(),inGraph->GetX()[pt],value);
      if (debug) cout << pt << " - 4 - " << ptMod << " - " << modValue << endl;
    }
    else if (pt > stopLoc+downRampLen) {
      outGraph->SetPoint(outGraph->GetN(),inGraph->GetX()[pt],0);
      if (debug) cout << pt << " - 5" << endl;
    }
  }

  if (debug) cout << " outGraph->GetN()=" << outGraph->GetN() << endl;
  return outGraph;

}
/*-----------
 */



