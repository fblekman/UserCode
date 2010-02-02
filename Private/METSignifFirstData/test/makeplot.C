{
  TChain *ch = new TChain("metsignifntup2/tree");
  ch->Add("/tmp/fblekman/*.root");
  TChain *ch2 = new TChain("metsignifntup1/tree");
  ch2->Add("/tmp/fblekman/*.root");
  ch2->Draw("met.met");
  ch->Draw("met.met","","same");
//   new TCanvas();
//   ch->Draw("met.phi");
//   new TCanvas();
//   ch->Draw("met.signif");
//   new TCanvas();
//   ch->Draw("met.sumEt");

  


}
