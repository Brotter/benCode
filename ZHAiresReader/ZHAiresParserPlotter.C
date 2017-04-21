void ZHAiresParserPlotter(int energy=195) {


  stringstream name;
  name.str("");
  name << "ZHAiresParser_En" << energy << ".root";

  TFile *inFile = TFile::Open(name.str().c_str());
  cout << inFile << endl;
  
  int  off = 0;
  for (int sub=0; sub<4; sub++) {
    for (int ant=0; ant<8; ant++) {
      name.str("");
      name << "off" << off << "sub" << sub << "ant" << ant+1 << "_X";
      TGraph *currGraphX = (TGraph*)inFile->Get(name.str().c_str());
      name.str("");
      name << "off" << off << "sub" << sub << "ant" << ant+1 << "_Y";
      TGraph *currGraphY = (TGraph*)inFile->Get(name.str().c_str());
      name.str("");
      name << "off" << off << "sub" << sub << "ant" << ant+1 << "_Z";
      TGraph *currGraphZ = (TGraph*)inFile->Get(name.str().c_str());
      if (currGraphX == NULL || currGraphY == NULL || currGraphZ == NULL) {
	cout << name.str() << endl;
	continue;
      }
      
      currGraphY->SetLineColor(kRed);
      currGraphY->SetMarkerColor(kRed);
      
      currGraphZ->SetLineColor(kBlue);
      currGraphZ->SetMarkerColor(kBlue);
      
      if (off==0 && sub==0 && ant==0) {
	currGraphX->Draw("alp"); 
	currGraphX->GetYaxis()->SetRangeUser(-0.06,0.1);
	currGraphX->GetXaxis()->SetRangeUser(0,2000);

	TLegend *leg = new TLegend(0.7,0.9,0.9,0.6);
	leg->AddEntry(currGraphX,"Ex","lp");
	leg->AddEntry(currGraphY,"Ey","lp");
	leg->AddEntry(currGraphZ,"Ez","lp");
	leg->Draw();

      }
	else {
	  currGraphX->Draw("lpSame");
	  currGraphY->Draw("lpSame");
	  currGraphZ->Draw("lpSame");
	}


    }
  }
  


  return;

}
