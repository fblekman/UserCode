{

  gSystem->Load("histolib_cxx.so");
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  TFile theFile("PATLayer0_Output.fromAOD_full.root");  
  std::cout << "opening file" << std::endl;
  TTree *events = (TTree*)theFile.Get("Events");
  events->GetEntry();
 
// connection of products and branches
 std::vector<pat::Electron> patElectrons;
 std::vector<pat::Muon> patMuons;
 std::vector<pat::Jet> patJets;
 std::vector<pat::MET> patMet;
 
 // TBranch *eleBranch = events->GetBranch(events->GetAlias("selectedLayer1Electrons"));
 // eleBranch->SetAddress(&patElectrons);

 TBranch *muBranch = events->GetBranch(events->GetAlias("selectedLayer1Muons"));
 muBranch->SetAddress(&patMuons);

 TBranch *jetBranch = events->GetBranch(events->GetAlias("selectedLayer1Jets"));
 jetBranch->SetAddress(&patJets);

 TBranch *metBranch = events->GetBranch(events->GetAlias("selectedLayer1METs"));
 metBranch->SetAddress(&patMet);

  gStyle->SetOptStat(0);
  using namespace std;

  double muPtcut=20;
  double jetPtcut=30;
  double MZ=91.1876;

  TFile * myfile = new TFile("myfile.root","recreate");
  myfile->cd();
  TH1D * hist_muPt_precut = new TH1D("hist_muPt_precut","PAT muons",40,0,400);
  settitles(hist_muPt_precut,"#mu p_{T} (GeV/c)","entries / 10 GeV/c");
  TH1D *hist_muPt_postcut = (TH1D*) hist_muPt_precut->Clone("hist_muPt_postcut");
  TH1D * hist_muEta_precut = new TH1D("hist_muEta_precut","PAT muons",20,-2.5,2.5);
  change_color(hist_muPt_postcut,103);
  settitles(hist_muEta_precut,"#mu #eta","entries / 0.25");
  TH1D *hist_muEta_postcut = (TH1D*) hist_muEta_precut->Clone("hist_muEta_postcut");
  change_color(hist_muEta_postcut,103);

  TH1D * hist_jetPt_precut = new TH1D("hist_jetPt_precut","PAT jets",40,0,400);
  settitles(hist_jetPt_precut,"jet p_{T} (GeV/c)","entries / 10 GeV/c");
  TH1D *hist_jetPt_postcut = (TH1D*) hist_jetPt_precut->Clone("hist_jetPt_postcut");
  TH1D * hist_jetEta_precut = new TH1D("hist_jetEta_precut","PAT jets",20,-2.5,2.5);
  change_color(hist_jetPt_postcut,103);
  settitles(hist_jetEta_precut,"jet #eta","entries / 0.25");
  TH1D *hist_jetEta_postcut = (TH1D*) hist_jetEta_precut->Clone("hist_jetEta_postcut");
  change_color(hist_jetEta_postcut,103);
  TH1D *hist_dPhi_jetZ = new TH1D("hist_dPhi_jetZ","Z+jet #rightarrow #mu #mu jet",20, -TMath::Pi(), TMath::Pi());
  settitles(hist_dPhi_jetZ,"#Delta #phi (Z boson, jet) (rad)","entries/ 0.1 #pi");
  

  TH1D * hist_muZPt = new TH1D("hist_muZPt", "Z #rightarrow #mu #mu candidates", 20, 0, 200 );
  settitles(hist_muZPt,"Z p_{T} (GeV/c)","entries/10 GeV/c");
  TH1D * hist_muZPhi = new TH1D("hist_muZPhi", "Z #rightarrow #mu #mu candidates", 20, -TMath::Pi(), TMath::Pi());
  settitles(hist_muZPhi,"Z #phi (rad)","entries/ 0.1 #pi");
  TH1D * hist_muZM = new TH1D("hist_muZM","Z #rightarrow #mu #mu candidates", 20, 0, 200 );
  settitles(hist_muZM,"Z boson mass (GeV/c^{2})","entries / 10 GeV/c^{2}");


  std::cout << "number of events is : " << events->GetEntries() << std::endl;
for( unsigned int index = 0; index < events->GetEntries(); index++)
//for( unsigned int index = 0; index < 2; index++)
 {
    jetBranch->GetEntry(index);
    //    eleBranch->GetEntry(index);
    muBranch->GetEntry(index);
    metBranch->GetEntry(index);

    events->GetEntry(index,0);
    if(index<10 || index%100==0)
      std::cout << "event " << index <<"/" <<events->GetEntries() <<   std::endl;

    std::cout << "number of muons = " << patMuons.size() << std::endl;
    TLorentzVector  bestzcand, jetcand, bestjetcand;
    int nzcand=0;
    for(size_t i=0; i< patMuons.size();++i){
      hist_muPt_precut->Fill(patMuons[i].pt());
      hist_muEta_precut->Fill(patMuons[i].eta());
      if(patMuons[i].pt()<muPtcut)
	continue;
      hist_muPt_postcut->Fill(patMuons[i].pt());
      hist_muEta_postcut->Fill(patMuons[i].eta());
      for(size_t j=i+1; j< patMuons.size(); ++j){
	math::XYZTLorentzVector diMus;
	diMus = patMuons[i].p4()+patMuons[j].p4();
	std::cout << patMuons[i].pt() <<" " << patMuons[i].charge() << " " << patMuons[i].eta() << std::endl;
	std::cout << patMuons[j].pt() <<" " << patMuons[j].charge() << " " << patMuons[j].eta() << std::endl;
	std::cout << "z mass " << diMus.mass() << std::endl;
	if(nzcand==0)
	  bestzcand.SetPtEtaPhiE(diMus.pt(),diMus.eta(),diMus.phi(),diMus.e());
	else{
	  if(fabs(diMus.mass()-MZ)<fabs(bestzcand.M()-MZ)){
	    bestzcand.SetPtEtaPhiE(diMus.pt(),diMus.eta(),diMus.phi(),diMus.e());
	  }
	}
	nzcand++;
      }      
    }
    if(nzcand==0)
      continue;
    hist_muZPt->Fill(bestzcand.Pt());
    hist_muZPhi->Fill(bestzcand.Phi());
    hist_muZM->Fill(bestzcand.M());
    std::cout << "number of jets = " << patJets.size() << std::endl;
    int ngoodjets=0;
    for(size_t ij=0; ij<patJets.size();++ij){
      hist_jetPt_precut->Fill(patJets[ij].pt());
      hist_jetEta_precut->Fill(patJets[ij].eta());
      if(patJets[i].pt()<jetPtcut)
	continue;
      hist_jetPt_postcut->Fill(patJets[ij].pt());
      hist_jetEta_postcut->Fill(patJets[ij].eta());
      jetcand.SetPtEtaPhiE(patJets[ij].pt(),patJets[ij].eta(),patJets[ij].phi(),patJets[ij].e());
      if(ngoodjets==0)
	bestjetcand=jetcand;
      else{
	if(fabs(jetcand.DeltaPhi(bestzcand)-TMath::Pi())<fabs(bestjetcand.DeltaPhi(bestzcand)-TMath::Pi()))
	  bestjetcand=jetcand;
      }
    }

 } // end loop on events

  TCanvas *mus = new TCanvas();
  mus->Divide(2,1);
  mus->cd(1);
  hist_muPt_precut->Draw();
  hist_muPt_postcut->Draw("same");
  mus->cd(2);
  hist_muEta_precut->Draw();
  hist_muEta_postcut->Draw("same");
  
  TCanvas *jets = new TCanvas();
  jets->Divide(2,1);
  jets->cd(1);
  hist_jetPt_precut->Draw();
  hist_jetPt_postcut->Draw("same");
  jets->cd(2);
  hist_jetEta_precut->Draw();
  hist_jetEta_postcut->Draw("same");
  
  
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);
  hist_muZPt->Draw();
  c1->cd(2);
  hist_muZM->Draw();
  c1->cd(3);
  hist_muZPhi->Draw();
  c1->cd(4);
  hist_dPhi_jetZ->Draw();

  myfile->Write();

}

