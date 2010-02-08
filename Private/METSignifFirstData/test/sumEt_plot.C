{
  TChain *chdata = new TChain("t","t");
  chdata->Add("outputhistos_data.root");
  TChain *chmc = new TChain("t","t");
  chmc->Add("outputhistos_mc.root");
  
  Double_t scaler = chdata->GetEntries()/(Double_t) chmc->GetEntries();

  TH1D *hist1 = new TH1D("hist1","CMS at #sqrt{s}=900 GeV",100,0,100);
  hist1->SetXTitle("#Sigma E_{T} (GeV)");
  hist1->SetYTitle("events / 1 GeV");
  TH1D *hist2 = (TH1D*) hist1->Clone("hist2");
  hist1->SetLineWidth(2);
  hist1->SetMarkerStyle(20);
  hist1->SetMarkerSize(1.1*hist1->GetMarkerSize());
  hist1->Sumw2();
  hist1->SetLineColor(1);
  hist2->SetLineColor(2);
  hist2->SetFillColor(hist2->GetLineColor());
  hist2->SetFillStyle(1001);
  
  chdata->Draw("met.sumEt>>hist1","met.sumEt<100","goff");
  chmc->Draw("met.sumEt>>hist2","met.sumEt<100","goff");
  //  scaler = hist1->GetEntries()/hist2->GetEntries();
  scaler = hist1->GetMaximum()/hist2->GetMaximum();
  hist2->Scale(scaler);
  TCanvas *c1 = new TCanvas();
  c1->cd();
  c1->SetLogy();
  hist1->Draw("pe");
  hist2->Draw("histsame");
  hist1->Draw("psame");

  TLegend *leg = new TLegend(0.6,0.6,1.0,0.9,hist1->GetTitle());
  hist1->SetTitle("");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hist1,"data","pl");
  leg->AddEntry(hist2,"MC","l");
  leg->Draw("same");
  c1->Modified();
  c1->Print("sumEt_CMSSW341_900GeV.pdf");
  c1->Print("sumEt_CMSSW341_900GeV.jpg");

  TCanvas *canv = new TCanvas();
  canv->cd();
  TH1D *ratio = (TH1D*) hist1->Clone("ratio");
  ratio->Divide(hist2);
  ratio->SetYTitle("data/MC");
  ratio->Draw();
  
}
