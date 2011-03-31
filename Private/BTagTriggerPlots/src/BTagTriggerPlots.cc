// -*- C++ -*-
//
// Package:    BTagTriggerPlots
// Class:      BTagTriggerPlots
// 
/**\class BTagTriggerPlots BTagTriggerPlots.cc Private/BTagTriggerPlots/src/BTagTriggerPlots.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  local user
//         Created:  Wed Mar 30 11:14:35 CEST 2011
// $Id: BTagTriggerPlots.cc,v 1.1 2011/03/30 13:55:24 fblekman Exp $
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"



#include <TString.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <map>
#include <string>
//
// class declaration
//

class BTagTriggerPlots : public edm::EDAnalyzer {
   public:
      explicit BTagTriggerPlots(const edm::ParameterSet&);
      ~BTagTriggerPlots();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  void makeNewHistos(std::string bla);

      // ----------member data ---------------------------
  // input:
  //  jet info:
  std::string jetSrc_;
  std::string muSrc_;
  std::string trigSrc_;
  std::string trigname_;
  std::string montrigname_;

  std::map<std::string,TH1D*> histos;


  edm::Service<TFileService> fs;

  TTree *tree;

  Int_t nJet;

  Float_t jettrig_pt[20];
  Float_t jettrig_eta[20];
  Float_t jettrig_phi[20];
  Float_t jettrig_e[20];
  Float_t jet_pt[20];
  Float_t jet_eta[20];
  Float_t jet_phi[20];
  Float_t jet_e[20];
  Float_t jetmu_pt[20];
  Float_t jetmu_eta[20];
  Float_t jetmu_phi[20];
  Float_t jetmu_e[20];
  Float_t jet_mupt[20];
  Float_t jet_muptrel[20];
  Float_t jet_ssvhe[20];
  Float_t jet_tche[20];
  Float_t jet_ssvhp[20];
  Float_t jet_sltDisc[20];
  Float_t jet_tchp[20];
  Int_t runnum;
  Int_t npv;
  

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
BTagTriggerPlots::BTagTriggerPlots(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  jetSrc_=iConfig.getParameter<std::string>("jetSrc");
  muSrc_=iConfig.getParameter<std::string>("muSrc");
  trigSrc_=iConfig.getParameter<std::string>("trigSrc");
  trigname_=iConfig.getParameter<std::string>("trigName");
  montrigname_=iConfig.getParameter<std::string>("monTrigName");
  std::cout << "examining " << trigname_ << " vs " << montrigname_ <<  std::endl;

  makeNewHistos("any");
  
  histos["nevents"]= fs->make<TH1D>("nevents","number of events processed",1,0,1);

  tree = new TTree("tree","tree");
  tree->Branch("nJet",&nJet,"nJet/I");
  tree->Branch("runnum",&runnum,"runnum/I");
  tree->Branch("npv",&npv,"npv/I");
  tree->Branch("jetmu_pt",jetmu_pt,"jetmu_pt[nJet]/F");
  tree->Branch("jetmu_eta",jetmu_eta,"jetmu_eta[nJet]/F");
  tree->Branch("jetmu_phi",jetmu_phi,"jetmu_phi[nJet]/F");
  tree->Branch("jetmu_e",jetmu_e,"jetmu_e[nJet]/F");
  tree->Branch("jettrig_pt",jettrig_pt,"jettrig_pt[nJet]/F");
  tree->Branch("jettrig_eta",jettrig_eta,"jettrig_eta[nJet]/F");
  tree->Branch("jettrig_phi",jettrig_phi,"jettrig_phi[nJet]/F");
  tree->Branch("jettrig_e",jettrig_e,"jettrig_e[nJet]/F");
  tree->Branch("jet_pt",jet_pt,"jet_pt[nJet]/F");
  tree->Branch("jet_eta",jet_eta,"jet_eta[nJet]/F");
  tree->Branch("jet_phi",jet_phi,"jet_phi[nJet]/F");
  tree->Branch("jet_e",jet_e,"jet_e[nJet]/F");
  tree->Branch("jet_mupt",jet_mupt,"jet_mupt[nJet]/F");
  tree->Branch("jet_muptrel",jet_muptrel,"jet_muptrel[nJet]/F");
  tree->Branch("jet_ssvhe",jet_ssvhe,"jet_ssvhe[nJet]/F");
  tree->Branch("jet_tche",jet_tche,"jet_tche[nJet]/F");
  tree->Branch("jet_ssvhp",jet_ssvhp,"jet_ssvhp[nJet]/F");
  tree->Branch("jet_sltDisc",jet_sltDisc,"jet_sltDisc[nJet]/F");
  tree->Branch("jet_tchp",jet_tchp,"jet_tchp[nJet]/F");


}


BTagTriggerPlots::~BTagTriggerPlots()
{
 
  // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
BTagTriggerPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  std::cout << "analyze loop! " << std::endl;
  using namespace edm;
  using namespace trigger;

  edm::Handle<TriggerEvent> triggerObj;
  iEvent.getByLabel(trigSrc_,triggerObj); 
  if(!triggerObj.isValid()) { 
    edm::LogInfo("Status") << "Summary HLT object (TriggerEvent) not found, "
      "skipping event"; 
    return;
  }
  
  const trigger::TriggerObjectCollection & toc(triggerObj->getObjects());

  size_t itrig1=-1, itrig2=-1;
  bool foundsomething=false;
  bool foundmontrig=false;
  for ( size_t ia = 0; ia < triggerObj->sizeFilters(); ++ ia) {
    std::string fullname = triggerObj->filterTag(ia).encode();

      std::string name;
      size_t p = fullname.find_first_of(':');
      if ( p != std::string::npos) {
	name = fullname.substr(0, p);
      }
      else {
	name = fullname;
      }


      TString workerstring = name;
      if(workerstring.Contains(trigname_.c_str())){
	foundsomething=true;
	//	std::cout << "found trigger yaay!" << std::endl;
	//	k = triggerObj->filterKeys(ia);
	itrig1=ia;
      }
      if(workerstring.Contains(montrigname_.c_str())){
	foundmontrig=true;
	//	kb=triggerObj->filterKeys(ia);
	itrig2=ia;
      }
      //      std::cout << "filter " << ia << ", full name = " << fullname
      //		<< ", p = " << p 
      //		<< ", abbreviated = " << name  << std::endl;
  }
  //  LogDebug("Status") << "filling ... " ;
   
  if(!foundsomething)
    return;
  //  std::cout << triggerObj->precale() << std::endl;
  const trigger::Keys & k = triggerObj->filterKeys(itrig1);
  
  
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muSrc_,muons);
  
//   edm::Handle<edm::View<reco::Vertex> > vertices;
//   iEvent.getByLabel("offlinePrimaryVertex",vertices);
//   npv=vertices->size();

  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);
      
  nJet=0;

  std::vector<std::pair<std::string,float> > discrimcombo;
  
  TLorentzVector jetworker,muworker;

  for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
    for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
      //	  std::cout << reco::deltaR(toc[*ki].eta(),toc[*ki].phi(),jet->eta(),jet->phi()) << std::endl;
      //	  std::cout << toc[*ki].pt() << " " << toc[*ki].eta() << " " << jet->pt() << std::endl;
      if(reco::deltaR(toc[*ki].eta(),toc[*ki].phi(),jet->eta(),jet->phi())>0.3)
	continue;
      nJet++;
      //      discrimcombo= jet->getPairDiscri();
      //      for(size_t jj=0; jj<discrimcombo.size(); jj++)
      //	std::cout << discrimcombo[jj].first << std::endl;
      histos["jetPt_any"]->Fill(jet->pt());
      jettrig_pt[nJet-1]=-1;
      jettrig_eta[nJet-1]=-10;
      jettrig_phi[nJet-1]=-10;
      jettrig_e[nJet-1]=-1;
      jetmu_pt[nJet-1]=-1;
      jetmu_eta[nJet-1]=-10;
      jetmu_phi[nJet-1]=-10;
      jetmu_e[nJet-1]=-1;
      jet_pt[nJet-1]=jet->pt();
      jet_eta[nJet-1]=jet->eta();
      jet_phi[nJet-1]=jet->phi();
      jet_e[nJet-1]=jet->energy();
      jet_ssvhp[nJet-1]=jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
      jet_ssvhe[nJet-1]=jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      jet_tchp[nJet-1]=jet->bDiscriminator("trackCountingHighPurBJetTags");
      jet_tche[nJet-1]=jet->bDiscriminator("trackCountingHighEffBJetTags");
      jet_sltDisc[nJet-1]=jet->bDiscriminator("softMuonBJetTags");
      jetworker.SetPtEtaPhiE(jet->pt(),jet->eta(),jet->phi(),jet->energy());
      muworker.SetPxPyPzE(0,0,0,0);
      jet_mupt[nJet-1]=	jet_muptrel[nJet-1]=-1;
      if(jet-> tagInfoSoftLepton("softMuonBJetTags")){
	const reco::SoftLeptonTagInfo *infotag = jet->tagInfoSoftLepton("softMuonBJetTags");
	muworker.SetPtEtaPhiE(infotag->lepton(0)->pt(),infotag->lepton(0)->eta(),infotag->lepton(0)->phi(),infotag->lepton(0)->p());
	jet_mupt[nJet-1]=infotag->lepton(0)->pt();
	jet_muptrel[nJet-1]=infotag->properties(0).ptRel;
	jetmu_pt[nJet-1]=(jetworker+muworker).Pt();
	jetmu_eta[nJet-1]=(jetworker+muworker).Eta();
	jetmu_phi[nJet-1]=(jetworker+muworker).Phi();
	jetmu_e[nJet-1]=(jetworker+muworker).E();
      }
      
      if(!foundmontrig)
	continue;

      const trigger::Keys & kb = triggerObj->filterKeys(itrig2);
      for (trigger::Keys::const_iterator kib = kb.begin(); kib !=kb.end(); ++kib ) {
	if(reco::deltaR(toc[*kib].eta(),toc[*kib].phi(),jet->eta(),jet->phi())<0.3){
	  histos["jetPtTrig_any"]->Fill(jet->pt());
	  jettrig_pt[nJet-1]=toc[*kib].pt();
	  jettrig_eta[nJet-1]=toc[*kib].eta();
	  jettrig_phi[nJet-1]=toc[*kib].phi();
	  jettrig_e[nJet-1]=toc[*kib].energy();
	}
      }  
    }
  }
  runnum=iEvent.id().run();
  if(nJet>0)
    tree->Fill();
  histos["nevents"]->Fill(0.5);

}


// ------------ method called once each job just before starting event loop  ------------
void BTagTriggerPlots::makeNewHistos(std::string name){
  
  std::string name1="jetPt_"+name;
  std::string name2="jetPtTrig_"+name;
  std::string name3="jetPtEff_"+name;
  
  histos[name1] = fs->make<TH1D>(name1.c_str(),"jet p_{T}",100,0,400);
  histos[name2] = fs->make<TH1D>(name2.c_str(),"jet p_{T}",100,0,400);
  histos[name3] = fs->make<TH1D>(name3.c_str(),"jet p_{T}",100,0,400);
  histos[name1]->Sumw2();
  histos[name2]->Sumw2();
  histos[name3]->Sumw2();
}

void 
BTagTriggerPlots::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
BTagTriggerPlots::endJob() {
  //  tree->Write();

  histos["jetPtEff_any"]->Divide(histos["jetPtTrig_any"],histos["jetPt_any"],1,1,"B");
}

//define this as a plug-in
DEFINE_FWK_MODULE(BTagTriggerPlots);
