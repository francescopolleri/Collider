{
  TFile f("FlavTag.root");
  auto FlavTag = f.Get<TTree>("FlavTag");

  TH1D hb("hb","",40,-1.2,1.2);
  TH1D hc("hc","",40,-1.2,1.2);
  TH1D hl("hl","",40,-1.2,1.2);

  FlavTag->Draw("sv1_l3d>>hb","label==5");
  FlavTag->Draw("sv1_l3d>>hc","label==4");
  FlavTag->Draw("sv1_l3d>>hl","label==0");

  gStyle->SetOptStat(0);
  hb.Draw("SAME");
  hb.SetLineColor(2);
  hc.Draw("SAME");
  hc.SetLineColor(4);
  hl.Draw("SAME");
  hl.SetLineColor(1);

  TLegend l(0.7, 0.8, 0.9, 0.9);
  l.AddEntry(&hb,"b-jets");
  l.AddEntry(&hc,"c-jets");
  l.AddEntry(&hl,"l-jets");
  l.Draw();
}
