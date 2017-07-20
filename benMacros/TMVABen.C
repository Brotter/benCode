#include "AnitaTMVA.h" 

/*  Macro to demonstrate usage of TMVA with AnitaEventSummary's
 *  
 *  Note that the TMVA interface changed completely between ROOT 6.06 and ROOT 6.08. This macro
 *  defaults to the ROOT 6.08 interface, but it's easy to modify for ROOT 6.06. 
 *
 *  Please don't use this for real analysis, it's just meant to chose how to do things. In this case, cut_signal selects for a WAIS pulser and cut_bg selects for not a wais pulser  
 *
 *  Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 */ 


/* These are the names of our branches */ 
 
const char * vars[] = { "mapPeak","hilbertPeak","linPolFraction" }; 

/* These are the expressions for our branches */ 
const char * expr[] = {"peak.value[0][0]", "coherent.peakHilbert[0][0]",
		       "TMath::Sqrt(pow(coherent[0][0].Q,2) + pow(coherent[0][0].U,2) ) / coherent[0][0].I"}; 

/* These is our selection cuts.
 *
 * This is just an example, so I want to show how to both train TMVA and use it to apply cuts, so I'll use odd events for one and even for the other 
 *
 */
stringstream cut_signal,cut_bg,cut_eval;


void setupCuts() {
  string evenCut = "(Entry$ %2) == 1";
  string oddCut = "(Entry$ %2) == 0";
  string globalCuts1 = "peak.value[0][0] > 0 && flags.isRF==1 && peak[0][0].theta < 60";
  string globalCuts2 = "peak[0][0].theta > -60 && flags.isPayloadBlast == 0";
  string hPolCut = "flags.isHPolTrigger == 1";
  string waisPulser = "flags.pulser == 1";
  string noPulser = "flags.pulser == 0";    
  string cut450M = "!(coherent[0][0].peakFrequency[0] > 0.45 && coherent[0][0].peakFrequency[0] < 0.47  && coherent[0][0].peakPower[0] > 16.2 && coherent[0][0].peakPower[0] < 17.4)";
  string plus = " && ";
  
  cut_signal.str("");
  cut_signal << "(" << evenCut << plus << globalCuts1 << plus << globalCuts2 << plus << waisPulser;
  cut_signal << plus << cut450M << plus << hPolCut << ")";
  
  cut_bg.str("");
  cut_bg << "(" << evenCut << plus << globalCuts1 << plus << globalCuts2 << plus << noPulser;
  cut_bg  << plus << cut450M << plus << hPolCut << ")";
  
  
  cut_eval.str("");
  cut_eval << "(" << oddCut << plus << globalCuts1 << plus << globalCuts2  << plus << cut450M;
  cut_eval << plus << hPolCut << ")";
 
  return;
}

/* filename should be a ROOT file containg tree with AnitaEventSummary treename */ 
void TMVABen()
{

  string baseDir = "/Volumes/ANITA3Data/bigAnalysisFiles/findCosmicRays/01.19.17_22h/";


  /* This is just a lazy way of loading the tree in this case, but you could also add multiple files trivially */ 

  TChain signalChain("headTree"); 
  TChain backChain("headTree");
  TChain evalChain("headTree");


  stringstream name;
  for (int run=340; run<350; run++) {
    name.str("");
    name << baseDir << run << ".root";
    signalChain.Add(name.str().c_str());
    evalChain.Add(name.str().c_str());
  }
  for (int run=287; run<300; run++) {
    name.str("");
    name << baseDir << run << ".root";
    backChain.Add(name.str().c_str()); 
    evalChain.Add(name.str().c_str());
  }



  cout << "Got : " << evalChain.GetEntries() << " entries " << endl;

  /* Create our variable set */ 
  AnitaTMVA::MVAVarSet varset(3,vars,expr); 


  //set up the cuts
  setupCuts();

  /* We'll do everything in the same file here. */ 
  TFile *treeFile = TFile::Open("tmvaTrees.root","RECREATE"); 
  TTree *sigtree = AnitaTMVA::makeTMVATree(&signalChain, treeFile, "signal_in", varset, cut_signal.str().c_str()); 
  TTree *bgtree = AnitaTMVA::makeTMVATree(&backChain, treeFile, "bg_in", varset, cut_bg.str().c_str()); 
  TTree *evaltree = AnitaTMVA::makeTMVATree(&backChain, treeFile, "eval_in", varset, cut_eval.str().c_str()); 

  cout << "made Trees" << endl;

  /** Create the TMVA factory */ 
  TFile *tvmaFile = TFile::Open("tmva.root","RECREATE");
  TMVA::Factory factory("tmva_example", tvmaFile); 

  cout << "made factory" << endl;


  /** If you're using ROOT < 6.08, comment out everything with the data_loader and uncomment the following lines out */ 
  varset.setUpData(&factory); 
  factory.AddSignalTree(sigtree); 
  factory.AddBackgroundTree(bgtree); 
  factory.BookMethod(TMVA::Types::kFisher,"Fisher"); // book a Fisher discriminant 
  //  factory.BookMethod(TMVA::Types::kBDT,"BDT"); //also a BDT because cosmin had luck with it


  cout << "added stuff to factory" << endl;


  //TMVA::DataLoader loader; 
  //varset.setUpData(&loader); 
  //loader.AddSignalTree(sigtree); 
  //loader.AddBackgroundTree(bgtree); 
  //factory.BookMethod(&loader, TMVA::Types::kFisher,"Fisher"); // book a Fisher discriminant 


  //Tell TMVA to train, test and evaluate
  factory.TrainAllMethods(); 
  cout << "trained" << endl;
  factory.TestAllMethods(); 
  cout << "tested" << endl;
  factory.EvaluateAllMethods(); 
  cout << "evaluated" << endl;
  


  // now let's go ahead and tag our eval tree with the output. 
  // not sure why this is the default path on my system, maybe it's different on yours? 
  evaluateTMVA( evaltree, varset, "Fisher",  "weights/tmva_example_Fisher.weights.xml"); 

  cout << "done except closing files" << endl;


  cout << "done with everything" << endl;

  treeFile->Close();
  tvmaFile->Close();

}

