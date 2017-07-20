{

  TGraph *lines = new TGraph();
  lines->SetName("gLines");

  for (int i=0; i<16; i++) {
    lines->SetPoint(lines->GetN(),1.418938412e9,i*22.5);
    lines->SetPoint(lines->GetN(),1.42088e9,i*22.5);
    lines->SetPoint(lines->GetN(),1.418938412e9,i*22.5);
  }


  lines->Draw("lSame");
}
