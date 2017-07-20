{

  TChain *tdataPointed = new TChain("tdataPointed","tdataPointed");
  TChain *tdataTracing = new TChain("tdataTracing","tdataTracing");
  int startRun,endRun;
    
  stringstream name;
  for (int i=150; i<440; i++) {
    name.str("");
    name << "/Volumes/NO NAME/ANITA3/rootOutputs/output" << i << "_1.root";
    tdataPointed->Add(name.str().c_str());
    tdataTracing->Add(name.str().c_str());
  }
  
  
  
  //  tdataPointed->Draw("sourceLat:sourceLon:eventNumber","","colz");
  gStyle->SetOptStat(0);
  tdataTracing->SetMarkerStyle(0);
  tdataPointed->SetMarkerStyle(0);
  tdataPointed->SetMarkerColor(kRed);
  
  TH1I *hEventNumber = new TH1I("hEventNumber","hEventNumber",1000,0,90e6);
  tdataTracing->Draw("eventNumber>>hEventNumber");
  
  hEventNumber->Draw();
  
  return;
}
