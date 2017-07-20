{

  //Ben Rotter August 2016

  //This macro easily makes the "heatmap" featured in Abby's dissertation.
  // From her thesis:
  // "A two dimensional histogram over Antarctica with one entry per event weighted by the peak value
  // of the interferometric map, normalized by the total number of events within that bin"
  //
  //  This discribes a TProfile2D, which can be plotted over Antarctica with BenS's really nice aMap object



  // define the TChains where abby stores the data (most of these are actually TNtuples but whatever)
  TChain *ndata = new TChain("ndata","ndata");
  TChain *ndata2 = new TChain("ndata2","ndata2");
  TChain *ndata3 = new TChain("ndata3","ndata3");
  TChain *ndata4 = new TChain("ndata4","ndata4");
  TChain *tdataEvent = new TChain("tdataEvent","tdataEvent");
  TChain *tdataTracing = new TChain("tdataTracing","tdataTracing");



  //link all the processed runs to those chains
  stringstream name;

  //should start at 130 (when the flight starts), but only have Abby processed code starting at 150
  for (int i=150; i<440; i++) {
    //these runs are missing
    if ((i>=257 && i<=263) || (i==368) || (i==393)) continue;

    name.str("");
    name << "/Volumes/ANITA3Data/rootOutputs/output" << i << "_1.root";
    ndata->Add(name.str().c_str());
    ndata2->Add(name.str().c_str());
    ndata3->Add(name.str().c_str());
    ndata4->Add(name.str().c_str());
    tdataEvent->Add(name.str().c_str());
    tdataTracing->Add(name.str().c_str());
  }

  //Link all those trees and ntuples together
  tdataTracing->AddFriend(ndata,"ndata");
  tdataTracing->AddFriend(ndata2,"ndata2");
  tdataTracing->AddFriend(ndata3,"ndata3");
  tdataTracing->AddFriend(ndata4,"ndata4");

  //Figure out how many entries we are dealing with
  int numEntries = tdataTracing->GetEntries();

  
  //Make the AntarcticaMapPlotter object and populate it with a TProfile2D
  AntarcticaMapPlotter *aMap = new AntarcticaMapPlotter("mapPlotter","mapPlotter",1000,1000);
  aMap->addProfile("heatMap","heatMap",1000,1000);
  TProfile2D* heatMap = (TProfile2D*)aMap->getCurrentHistogram();


  //pull out the variables that we want to plot
  Double_t sourceLat,sourceLon;
  Float_t peakVal;
  tdataTracing->SetBranchAddress("sourceLat",&sourceLat);
  tdataTracing->SetBranchAddress("sourceLon",&sourceLon);
  ndata->SetBranchAddress("peakVal",&peakVal);

  //Fill up the histogram
  for (int entry=0; entry<numEntries; entry++) {
    if (entry%1000 == 0) cout << entry << " / " << numEntries << "\r";
    tdataTracing->GetEntry(entry);
    aMap->Fill(sourceLat,sourceLon,peakVal);
  }
   
  //Then plot it!
  aMap->DrawHist("colz");

  //All done!
  return;
}
