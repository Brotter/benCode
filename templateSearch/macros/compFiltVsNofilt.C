/*

  I processed decimated data sets, one with adaptive sin subtraction on, and one with no filtering at all.

  The sin subtraction one takes like four times longer, so I'm going to look to see what it gains us.

 */




void compFiltVsNofilt() {

  string baseDir = "/Volumes/ANITA3data/bigAnalysisFiles/templateSearch/";
  string noFiltDir = baseDir + "06.06.17_13h/";
  string filtDir = baseDir + "06.07.17_17h/";

  TChain *noFiltTree = new TChain("summaryTree");
  TChain *filtTree = new TChain("summaryTree");

  


  stringstream name;
  for (int run=130; run<440; run++) {
    //except "bad" runs...
    if (run >= 214 && run <= 221) continue;
    if (run >= 257 && run <= 263) continue;

    name.str("");
    name << noFiltDir << run << ".root";
    noFiltTree->Add(name.str().c_str());
    name.str("");
    name << filtDir << run << ".root";
    filtTree->Add(name.str().c_str());
  }

  cout << "Imported trees" << endl;


  AnitaEventSummary *filtSummary = NULL;
  AnitaEventSummary *noFiltSummary = NULL;
  double filtTempH,filtTempV,noFiltTempH,noFiltTempV;

  noFiltTree->SetBranchAddress("eventSummary",&noFiltSummary);
  noFiltTree->SetBranchAddress("templateValueH",&noFiltTempH);
  noFiltTree->SetBranchAddress("templateValueV",&noFiltTempV);
  filtTree->SetBranchAddress("eventSummary",filtSummary);
  filtTree->SetBranchAddress("templateValueH",&filtTempH);
  filtTree->SetBranchAddress("templateValueV",&filtTempV);


  cout << "Checking filtered tree..." << endl;
  int lenFiltTree = filtTree->GetEntries();
  cout << "Checking non-filtered tree..." << endl;
  int lenNoFiltTree = noFiltTree->GetEntries();


  cout << "filtTree->GetEntries() = " << lenFiltTree << endl;
  cout << "noFiltTree->GetEntries() = " << lenNoFiltTree << endl;

  

  //make storage histograms
  TH2D *nFiltH = new TH2D("nFiltH","Sin Subtract Filtered - HPol Noise;Interferometric Map Peak; Template Correlation",
			   100,0,50,100,0,1);

  TH2D *nFiltV = new TH2D("nFiltV","Sin Subtract Filtered - VPol Noise;Interferometric Map Peak; Template Correlation",
			   100,0,50,100,0,1);

  TH2D *nNoFiltH = new TH2D("nFiltH","No Filter - HPol Noise;Interferometric Map Peak; Template Correlation",
			   100,0,50,100,0,1);

  TH2D *nNoFiltV = new TH2D("nFiltV","No Filter - VPol Noise;Interferometric Map Peak; Template Correlation",
			   100,0,50,100,0,1);

  
  TH2D *pFiltH = new TH2D("pFiltH","Sin Subtract Filtered - HPol Pulse;Interferometric Map Peak; Template Correlation",
			   100,0,50,100,0,1);

  TH2D *pFiltV = new TH2D("pFiltV","Sin Subtract Filtered - VPol Pulse;Interferometric Map Peak; Template Correlation",
			   100,0,50,100,0,1);

  TH2D *pNoFiltH = new TH2D("pFiltH","No Filter - HPol Pulse;Interferometric Map Peak; Template Correlation",
			   100,0,50,100,0,1);

  TH2D *pNoFiltV = new TH2D("pFiltV","No Filter - VPol Pulse;Interferometric Map Peak; Template Correlation",
			   100,0,50,100,0,1);



  //go through filtered tree
  for (int entry=0; entry<lenFiltTree; entry++) {
    if (entry%1000 == 0) cout << "FiltTree: " << entry << "/" << lenFiltTree << "\r";
    filtTree->GetEntry();
    if (filtSummary->flags.isVPolTrigger) {
      if (filtSummary->
    





  return;

}
