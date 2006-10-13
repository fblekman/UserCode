// -*- C++ -*-
//
// Package:    IcTauFakeRateAnalyzer
// Class:      IcTauFakeRateAnalyzer
// 
/**\class IcTauFakeRateAnalyzer IcTauFakeRateAnalyzer.cc MyTauAndHLTAnalyzer/IcTauFakeRateAnalyzer/src/IcTauFakeRateAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Dr Freya Blekman
//         Created:  Fri Oct  6 13:10:23 BST 2006
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// necessary objects:
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BTauReco/interface/IsolatedTauTagInfo.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include <DataFormats/Common/interface/RefVector.h>
#include <DataFormats/Common/interface/RefVectorIterator.h>
#include <DataFormats/Common/interface/Ref.h>
// basic cpp:

// basic cpp:
#include <iostream>
#include <string>
// some root includes
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
//
// class decleration
//

class IcTauFakeRateAnalyzer : public edm::EDAnalyzer {
   public:
      explicit IcTauFakeRateAnalyzer(const edm::ParameterSet&);
      ~IcTauFakeRateAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  std::string muonname;
  std::string jetname;
  std::string elename;
  std::string tauname;
  std::string rootfilename;
  
  int nEvt;// used to count the number of events

  
  double m_Riso;
  double m_Rsig;  
  double m_Rm;
  double m_Ptmin_lt;
  double m_Pt_tr;
  double m_z_ltr;

  // to be used for root output tree
  TFile *thefile;
  TH1D *jetptspectrum;
  TH1D *jetetaspectrum;
  TH1D *jetphispectrum;
  TH1D *jetptspectrumLt10Eta3;;
  TH1D *jetetaspectrumLt10Eta3;;
  TH1D *jetphispectrumLt10Eta3;;
  TH1D *njetsperevent;
  TH1D *njetspereventLt10Eta3;
  
  TH1D *ptZ_ele;
  TH1D *ptZ_mu;
  TH1D *mZ_ele;
  TH1D *mZ_mu;
  TH1D *eleptspectrum;
  TH1D *eleetaspectrum;
  TH1D *elephispectrum;
  TH1D *eleptspectrumLt5Eta25;
  TH1D *eleetaspectrumLt5Eta25;
  TH1D *elephispectrumLt5Eta25;
  TH1D *neleperevent;
  TH1D *nelepereventLt5Eta25;
  TH1D *muptspectrum;
  TH1D *muetaspectrum;
  TH1D *muphispectrum;
  TH1D *muptspectrumLt5Eta25;
  TH1D *muetaspectrumLt5Eta25;
  TH1D *muphispectrumLt5Eta25;
  TH1D *nmuperevent;
  TH1D *nmupereventLt5Eta25;
  // tree
  TTree *thetree;
  // with branches
  // four-momenta of Z candidate
  double z_pt;
  double z_phi;
  double z_eta;
  double z_e;
  double z_m;
  int z_leptonid;
  // some jet stuff
  int njets;
  double jet_pt[10];
  double jet_phi[10];
  double jet_eta[10];
  double jet_e[10];
  int ntaus;
  double tau_pt[10];
  double tau_phi[10];
  double tau_eta[10];
  double tau_e[10];
  double tau_leadingtrkpt[10];
  double tau_leadingtrkdR[10];
  int tau_tag[10];
  int tau_ntrk[10];
  int ntracks;
  double trk_pt[50];
  double trk_phi[50];
  double trk_eta[50];
  double trk_e[50];
  double trk_z[50];
  double trk_dR[50];
  int trk_taumatch[50];
  
 
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
IcTauFakeRateAnalyzer::IcTauFakeRateAnalyzer(const edm::ParameterSet& iConfig):
  m_Riso(iConfig.getUntrackedParameter<double>("Riso",0.4)),
  m_Rsig(iConfig.getUntrackedParameter<double>("Rsig",0.07)),
  m_Rm(iConfig.getUntrackedParameter<double>("Rm",0.1)),
  m_Ptmin_lt(iConfig.getUntrackedParameter<double>("Ptmin_lt",10.0)),
  m_Pt_tr(iConfig.getUntrackedParameter<double>("Pt_tr",1.0)),         
  m_z_ltr(iConfig.getUntrackedParameter<double>("z_ltr",20)),
  muonname(iConfig.getUntrackedParameter<std::string>("MuonHandle","globalMuons")),
  jetname(iConfig.getUntrackedParameter<std::string>("JetHandle","iterativeCone5CaloJets")),
  elename(iConfig.getUntrackedParameter<std::string>("ElectronHandle","pixelMatchElectrons")),
  tauname(iConfig.getUntrackedParameter<std::string>("TauHandle","coneIsolation")),
  rootfilename(iConfig.getUntrackedParameter<std::string>("RootFileName","ExampleSmallAnalysisOutput.root")),
  nEvt(0)
{
   //now do what ever initialization is needed
  thefile = new TFile(rootfilename.c_str(),"recreate");
  thefile->cd();
  TDirectory *jetdir = new TDirectory("jets","jets");
  jetdir->cd();
  jetptspectrum = new TH1D("jetptspectrum","jet p_{T}",100,0,100);
  jetetaspectrum = new TH1D("jetetaspectrum","jet #eta",100,-6,6);
  jetphispectrum = new TH1D("jetphispectrum","jet #phi",100,-TMath::Pi(),TMath::Pi());
  jetptspectrumLt10Eta3 = new TH1D("jetptspectrumLt10Eta3","jet p_{T} (E_{T} > 10 GeV, |#eta| < 3 )",100,0,100);
  jetetaspectrumLt10Eta3 = new TH1D("jetetaspectrumLt10Eta3","jet #eta (E_{T} > 10 GeV, |#eta| < 3 )",100,-6,6);
  jetphispectrumLt10Eta3 = new TH1D("jetphispectrumLt10Eta3","jet #phi (E_{T} > 10 GeV, |#eta| < 3 )",100,-TMath::Pi(),TMath::Pi());
  njetsperevent = new TH1D("njetsperevent","jet multiplicity (no cuts)", 50,0,50);
  njetspereventLt10Eta3 = new TH1D("njetspereventLt10Eta3","jet multiplicity (E_{T} > 10 GeV, |#eta| < 3 )",50,0,50);
  thefile->cd();
  TDirectory *zs = new TDirectory("Z_candidates","Z_candidates");
  zs->cd();
  
  eleptspectrum = new TH1D("eleptspectrum","ele p_{T}",100,0,100);
  eleetaspectrum = new TH1D("eleetaspectrum","ele #eta",100,-6,6);
  elephispectrum = new TH1D("elephispectrum","ele #phi",100,-TMath::Pi(),TMath::Pi());
  eleptspectrumLt5Eta25 = new TH1D("eleptspectrumLt5Eta25","ele p_{T} (E_{T} > 5 GeV, |#eta| < 2.5 )",100,0,100);
  eleetaspectrumLt5Eta25 = new TH1D("eleetaspectrumLt5Eta25","ele #eta (E_{T} > 5 GeV, |#eta| < 2.5 )",100,-6,6);
  elephispectrumLt5Eta25 = new TH1D("elephispectrumLt5Eta25","ele #phi (E_{T} > 5 GeV, |#eta| < 2.5 )",100,-TMath::Pi(),TMath::Pi());
  muptspectrumLt5Eta25 = new TH1D("muptspectrumLt5Eta25","mu p_{T} (E_{T} > 5 GeV, |#eta| < 2.5 )",100,0,100);
  muetaspectrumLt5Eta25 = new TH1D("muetaspectrumLt5Eta25","mu #eta (E_{T} > 5 GeV, |#eta| < 2.5 )",100,-6,6);
  muphispectrumLt5Eta25 = new TH1D("muphispectrumLt5Eta25","mu #phi (E_{T} > 5 GeV, |#eta| < 2.5 )",100,-TMath::Pi(),TMath::Pi());
  muptspectrum = new TH1D("muptspectrum","mu p_{T}",100,0,100);
  muetaspectrum = new TH1D("muetaspectrum","mu #eta",100,-6,6);
  muphispectrum = new TH1D("muphispectrum","mu #phi",100,-TMath::Pi(),TMath::Pi());
  ptZ_ele = new TH1D("ptZ_ele","p_{T} of Z candidate (e^{+}e^{-})",200,0,200);
  TH1D *ptZ_ele_chargeinc = (TH1D*) ptZ_ele->Clone("ptZ_ele_chargeinc");
  ptZ_mu = new TH1D("ptZ_mu","p_{T} of Z candidate (#mu^{+}#mu^{-})",200,0,200);
  TH1D *ptZ_mu_chargeinc = (TH1D*) ptZ_mu->Clone("ptZ_mu_chargeinc");
  mZ_ele = new TH1D("mZ_ele","mass of Z candidate (e^{+}e^{-})",200,0,200);
  TH1D *mZ_ele_chargeinc = (TH1D*) mZ_ele->Clone("mZ_ele_chargeinc");
  mZ_mu = new TH1D("mZ_mu","mass of Z candidate (#mu^{+}#mu^{-})",200,0,200);
  TH1D *mZ_mu_chargeinc = (TH1D*) mZ_mu->Clone("mZ_mu_chargeinc");
  neleperevent = new TH1D("neleperevent","electron multiplicity (no cuts)", 50,0,50);
  nelepereventLt5Eta25 = new TH1D("nelepereventLt5Eta2.5","electron multiplicity (E_{T} > 5 GeV, |#eta| < 2.5 )",50,0,50);
  nmuperevent = new TH1D("nmuperevent","#mu multiplicity (no cuts)", 50,0,50);
  nmupereventLt5Eta25 = new TH1D("nmupereventLt5Eta2.5","#mu multiplicity (E_{T} > 5 GeV, |#eta| < 2.5 )",50,0,50);
 
  thetree = new TTree("thetree","thetree");
  thetree->Branch("z_pt",&z_pt,"z_pt/D");
  thetree->Branch("z_phi",&z_phi,"z_phi/D");
  thetree->Branch("z_eta",&z_eta,"z_eta/D");
  thetree->Branch("z_e",&z_e,"z_e/D");
  thetree->Branch("z_m",&z_m,"z_m/D");
  thetree->Branch("z_leptonid",&z_leptonid,"z_leptonid/I");
  thetree->Branch("njets",&njets,"njets/I");
  thetree->Branch("jet_pt",jet_pt,"jet_pt[njets]/D");
  thetree->Branch("jet_phi",jet_phi,"jet_phi[njets]/D");
  thetree->Branch("jet_eta",jet_eta,"jet_eta[njets]/D");
  thetree->Branch("jet_e",jet_e,"jet_e[njets]/D");
  thetree->Branch("ntaus",&ntaus,"ntaus/I");
  thetree->Branch("tau_pt",tau_pt,"tau_pt[ntaus]/D");
  thetree->Branch("tau_phi",tau_phi,"tau_phi[ntaus]/D");
  thetree->Branch("tau_eta",tau_eta,"tau_eta[ntaus]/D");
  thetree->Branch("tau_e",tau_e,"tau_e[ntaus]/D");
  thetree->Branch("tau_leadingtrkpt",tau_leadingtrkpt,"tau_leadingtrkpt[ntaus]/D");
  thetree->Branch("tau_leadingtrkdR",tau_leadingtrkdR,"tau_leadingtrkdR[ntaus]/D");
  thetree->Branch("tau_tag",tau_tag,"tau_tag[ntaus]/I");
  thetree->Branch("tau_ntrk",tau_ntrk,"tau_ntrk[ntaus]/I");
  thetree->Branch("ntracks",&ntracks,"ntracks/I");
  thetree->Branch("trk_pt",trk_pt,"trk_pt[ntracks]/D");
  thetree->Branch("trk_phi",trk_phi,"trk_phi[ntracks]/D");
  thetree->Branch("trk_eta",trk_eta,"trk_eta[ntracks]/D");
  thetree->Branch("trk_e",trk_e,"trk_e[ntracks]/D");
  thetree->Branch("trk_z",trk_z,"trk_z[ntracks]/D");
  thetree->Branch("trk_dR",trk_dR,"trk_dR[ntracks]/D");
  thetree->Branch("trk_taumatch",trk_taumatch,"trk_taumatch[ntracks]/I");
  }


IcTauFakeRateAnalyzer::~IcTauFakeRateAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
IcTauFakeRateAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
 using namespace edm;

   // this analyzer produces a small root file with basic candidates and some MC information
   // some additional print statements
 nEvt++;
 if((nEvt%10==0 && nEvt<=100)||(nEvt%100==0 && nEvt>100))
   std::cout << "reading event " << nEvt << std::endl;
   


 Handle<HepMCProduct> evt;
 iEvent.getByLabel("VtxSmeared", evt);
 
 TLorentzVector z_momentum,l_momplus,l_mommin;
//  HepMC::GenEvent * myGenEvent = new HepMC::GenEvent(*(evt->GetEvent()));
//  for ( HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin();
//        p != myGenEvent->particles_end() ; ++p ) {
//    if(abs((*p)->pdg_id())==23)
//      std::cout << "found Z!" << std::endl;
//    if(nEvt==1)
//      std::cout << (*p)->pdg_id() << std::endl;
//  }

 Handle<reco::ElectronCollection> electrons;
 iEvent.getByLabel(elename,electrons);
 bool foundzele=false;
 bool foundzmu=false;
 int elecounter=0;
 if(electrons->size()>=0){
   reco::ElectronCollection::const_iterator ele;
   reco::ElectronCollection::const_iterator ele2;
   for( ele = electrons->begin(); ele != electrons->end(); ++ ele ) {
     
     eleptspectrum->Fill(ele->et());
     eleetaspectrum->Fill(ele->eta());
     elephispectrum->Fill(ele->phi());
     if(ele->et()>5. && fabs(ele->eta())<2.5){
       elecounter++;
       eleptspectrumLt5Eta25->Fill(ele->et());
       eleetaspectrumLt5Eta25->Fill(ele->eta());
       elephispectrumLt5Eta25->Fill(ele->phi());
     }
     for(ele2 = electrons->begin() ; ele2 != electrons->end(); ++ ele2 ){
       if(ele2!=ele){
	 if(ele2->et()>5. && fabs(ele2->eta())<2.5){
	   std::cout <<" electrons: " << ele->energy() << " "  << ele2->energy() << std::endl;
	   l_momplus.SetXYZT(ele->px(),ele->py(),ele->pz(),ele->energy());
	   l_mommin.SetXYZT(ele2->px(),ele2->py(),ele2->pz(),ele2->energy());
	   z_momentum = l_momplus + l_mommin;
	   mZ_ele->Fill(z_momentum.M());
	   foundzele=true;
	 }  
       }
     }     
   }
 }
 
 neleperevent->Fill(electrons->size());
 nelepereventLt5Eta25->Fill(elecounter);
 Handle<reco::TrackCollection> muons;
 
 //   iEvent.getByLabel("globalMuons",muons);
 iEvent.getByLabel(muonname,muons);
 int mucounter=0;
 if(muons->size()>=0){
   reco::TrackCollection::const_iterator muon;
   reco::TrackCollection::const_iterator muon2;
   
   for( muon = muons->begin(); muon != muons->end() ; ++ muon ) {
     muptspectrum->Fill(muon->pt());
     muetaspectrum->Fill(muon->eta());
     muphispectrum->Fill(muon->phi());
     if(muon->pt()>5. && fabs(muon->eta())<2.5){
       mucounter++;
       muptspectrumLt5Eta25->Fill(muon->pt());
       muetaspectrumLt5Eta25->Fill(muon->eta());
       muphispectrumLt5Eta25->Fill(muon->phi());
     }
     for( muon2 = muons->begin(); muon2 != muons->end() ; ++ muon2 ) {
       if(muon != muon2){
	 if(muon2->pt()>5. && fabs(muon2->eta())<2.5){
	   std::cout <<" muons: " << muon->p() << " "  << muon2->p() << std::endl;
	   l_momplus.SetXYZT(muon->px(),muon->py(),muon->pz(),muon->p());
	   l_mommin.SetXYZT(muon2->px(),muon2->py(),muon2->pz(),muon2->p());
	   z_momentum = l_momplus + l_mommin;
	   mZ_mu->Fill(z_momentum.M());
	   foundzmu=true;
	 }
       }
     }
   }
 } 
 nmuperevent->Fill(muons->size());
 nmupereventLt5Eta25->Fill(mucounter);
 // fill ntuple info:
 if(foundzmu|| foundzele){
   z_pt=z_momentum.Pt();
   z_phi=z_momentum.Phi();
   z_eta=z_momentum.Eta();
   z_e=z_momentum.E();
   z_m=z_momentum.M();
   if(foundzmu)
     z_leptonid=13;
   else
     z_leptonid=11;
 }
 Handle<reco::CaloJetCollection> jets;
 iEvent.getByLabel(jetname,jets);
 reco::CaloJetCollection::const_iterator jet;
 // std::cout << "number of jets: " << jets->size() << std::endl;
 njetsperevent->Fill(jets->size());
 int jetcounter =0;
 for( jet = jets->begin(); jet != jets->end(); ++ jet ) {
   //   std::cout << "jet pt: " << jet->et() << ", eta " << jet->eta() << ", phi " <<jet->phi() << ", energy " << jet->energy() << std::endl;
   jetptspectrum->Fill(jet->et());
   jetetaspectrum->Fill(jet->eta());
   jetphispectrum->Fill(jet->phi());
   if(jet->et()>10. && fabs(jet->eta())<3){
     if(jetcounter<10){
       jet_pt[jetcounter]=jet->et();
       jet_eta[jetcounter]=jet->eta();
       jet_phi[jetcounter]=jet->phi();
       jet_e[jetcounter]=jet->energy();
       njets=jetcounter;
     }
     jetcounter++;
     jetptspectrumLt10Eta3->Fill(jet->et());
     jetetaspectrumLt10Eta3->Fill(jet->eta());
     jetphispectrumLt10Eta3->Fill(jet->phi());
   }
 }
 njetspereventLt10Eta3->Fill(jetcounter);


 // do tau stuff

 Handle<reco::IsolatedTauTagInfoCollection> tauTagInfoHandle;
 iEvent.getByLabel("coneIsolation",tauTagInfoHandle);
 reco::IsolatedTauTagInfoCollection::const_iterator tau;
 ntaus=0;
 ntracks=0;
 for(tau=tauTagInfoHandle->begin();tau!=tauTagInfoHandle->end();++tau){
   
   if(tau->discriminator()){
     // access the tau jet info with tau->jet()
     
       tau_pt[ntaus]=tau->jet().et();
       tau_eta[ntaus]=tau->jet().eta();
       tau_phi[ntaus]=tau->jet().phi();
       TVector3 jetvector(tau->jet().px(),tau->jet().py(),tau->jet().pz());
       TVector3 leadingtrackvec(tau->leadingSignalTrack(0.45,6)->momentum().X(),tau->leadingSignalTrack(0.45,6)->momentum().Y(),tau->leadingSignalTrack(0.45,6)->momentum().Z());
       tau_e[ntaus]=tau->jet().energy();
       tau_leadingtrkpt[ntaus]=tau->leadingSignalTrack(0.45,6)->pt();
       tau_leadingtrkdR[ntaus]=jetvector.DeltaR(leadingtrackvec);
       // and examine some tau properties
       reco::TrackRefVector RecoTauJetTracks=tau->selectedTracks();//Returns All tracks in jet
       reco::TrackRefVector::iterator rvi=RecoTauJetTracks.begin();
       tau_ntrk[ntaus]=0;
       tau_tag[ntaus]=0;
       if(tau->discriminator(m_Rm,m_Rsig,m_Riso,m_Ptmin_lt,m_Pt_tr))
	 tau_tag[ntaus]=1;
       
       for(;rvi!=RecoTauJetTracks.end();rvi++)
	 {
	   trk_pt[ntracks]=(*rvi)->pt();
	   trk_eta[ntracks]=(*rvi)->eta();
	   trk_phi[ntracks]=(*rvi)->phi();
	   trk_e[ntracks]=(*rvi)->p();
	   trk_taumatch[ntracks]=ntaus;
	   trk_z[ntracks]=(*rvi)->dz(); 
	   TVector3 trackmom((*rvi)->momentum().X(),(*rvi)->momentum().Y(),(*rvi)->momentum().Z());

	   trk_dR[ntracks]=jetvector.DeltaR(trackmom);
	   tau_ntrk[ntaus]++;
	   ntracks++;
	   //	   cout<<"Track Pt::"<<(*rvi)->pt()<<endl;
	 }
       ntaus++;
   }
 } 
 // and last but not least:
 // fill the tree
 // maybe have some criterion
 if((foundzele==1||foundzmu==1)&&njets>0)
   thetree->Fill();
}



// ------------ method called once each job just before starting event loop  ------------
void 
IcTauFakeRateAnalyzer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
IcTauFakeRateAnalyzer::endJob() { 
  std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "analyzed " << nEvt << " events " << std::endl;
  std::cout << "writing information into file: " << thefile->GetName() << std::endl;
  thetree->Write();
  thefile->Write();
  thefile->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(IcTauFakeRateAnalyzer)
