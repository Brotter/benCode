
//I need to make this a library or something
/*

  And here it is I guess.  I think you can just do #include "whatever.C"

*/



TH1* makeNormCumulative(TH1* inHist,Bool_t forward=true) {
  /*
    Makes a "normalized cumulative cut fraction" graph I think

    Should let me set cuts on a specific amount of "reduction"

    forward specifies the direction of the cumulative sum
   */


  TH1* copyHist = (TH1*)inHist->Clone();
  
  double integral = 0;
  for (int i=0; i<copyHist->GetNbinsX(); i++) {
    double value = copyHist->GetBinContent(i);
    integral += value;
  }

  copyHist->Scale(1./integral);
  TH1* outHist = (TH1*)copyHist->GetCumulative(forward);

  for (int i=0; i<copyHist->GetNbinsX()+1; i++) {
    double value = outHist->GetBinContent(i);
    outHist->SetBinContent(i,1.-value);

  }
  delete copyHist;

  if (forward) outHist->GetYaxis()->SetTitle("Surviving Fraction");
  else         outHist->GetYaxis()->SetTitle("Fraction of Events Cut");
  return outHist;
}



TH1* makeNormHist(TH1* inHist) {
  /*

    Make a normalized histogram

  */

  TH1* copyHist = (TH1*)inHist->Clone();
  copyHist->GetYaxis()->SetTitle("Occupancy");

  
  double integral = 0;
  for (int i=0; i<copyHist->GetNbinsX(); i++) {
    double value = copyHist->GetBinContent(i);
    integral += value;
  }

  copyHist->Scale(1./integral);

  return copyHist;
}

