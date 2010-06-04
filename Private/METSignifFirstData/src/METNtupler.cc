// -*- C++ -*-
//
// Package:    METSignifFirstDataTest
// Class:      METSignifFirstDataTest
// 
/**\class METSignifFirstDataTest METSignifFirstDataTest.cc Private/METSignifFirstDataTest/src/METSignifFirstDataTest.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Freya Blekman,6 R-029,+41227678914,
//         Created:  Sun Apr 18 22:45:10 CEST 2010
// $Id: METNtupler.cc,v 1.3 2010/06/03 09:22:03 fblekman Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/InputTag.h"


#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"


#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoJets/JetAlgorithms/interface/JetIDHelper.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include <map>
#include <vector>
#include <iostream>

#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"

//
// class declaration
//

class METNtupler : public edm::EDAnalyzer {
public:
  explicit METNtupler(const edm::ParameterSet&);
  ~METNtupler();
  
private:
  virtual void beginJob() ;  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag src_;
  edm::InputTag jetsrc_;
  bool  usePF_;  
  bool usePAT_;
  reco::helper::JetIDHelper *jetID;


  std::map<std::string,TH1D*> histocontainer_; // simple map to contain all histograms. Histograms are booked in the beginJob() method
  

  TTree *tree_;
  TTree *pftree_;
  struct stuff{
    Float_t met;
    Float_t sig;
    Float_t phi;
    Float_t sumEt;
    Float_t prob;
  };
  stuff container_;

  struct eventstuff{
    Int_t runnum;
    Int_t evnum;
    Int_t lumiblock;
    Int_t eventclass;
  };
  eventstuff event_;
  struct pfstuff{
    Float_t p;
    Float_t pt;
    Float_t dpt;
    Float_t dphi;
    Float_t e;
    Float_t et;
    Float_t phi;
    Float_t eta;
    Float_t ehad;
    Float_t eem;
    Int_t pid;    
  };
  pfstuff pfccontainer_;
  Int_t numJet;
  Float_t jetPt[41];
  Float_t jetPhi[41];
  Float_t jetEta[41];
  Float_t jetE[41];
  Float_t jetEMF[41];
  Float_t jetFHPD[41];
  Float_t jetN90Hits[41];

  Int_t numEle;
  Float_t elePt[11];
  Float_t eleSCEta[11];
  Float_t eleEta[11];
  Float_t elePhi[11];
  Float_t eleE[11];
  Float_t eleD0[11];
  Float_t eleD0Sig[11];
  Float_t eleHcalIso[11];
  Float_t eleTrkIso[11];
  Float_t eleEcalIso[11];
      
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
METNtupler::METNtupler(const edm::ParameterSet& iConfig):
  src_(iConfig.getParameter<edm::InputTag>("src")), 
  jetsrc_(iConfig.getParameter<edm::InputTag>("jetsrc")),
  usePF_(iConfig.getUntrackedParameter<bool>("usePF",false)),
  usePAT_(iConfig.getUntrackedParameter<bool>("usePAT",false))
{
  jetID = new reco::helper::JetIDHelper(iConfig.getParameter<edm::ParameterSet>("JetIDParams"));

  edm::Service<TFileService> fs;
  histocontainer_["eventcount"]=fs->make<TH1D>("eventcount","events processed",1,-0.5,+0.5);
  
  tree_= new TTree("t","t");
  tree_->Branch("ev",&event_,"run/I:ev/I:lumi/I:eventclass/I");
  tree_->Branch("met",&container_,"met/F:sig/F:phi/F:sumEt/F");

  tree_->Branch("numJet",&numJet,"numJet/I");
  tree_->Branch("jetPt",jetPt,"jetPt[numJet]/F");
  tree_->Branch("jetPhi",jetPhi,"jetPhi[numJet]/F");
  tree_->Branch("jetEta",jetEta,"jetEta[numJet]/F");
  tree_->Branch("jetE",jetE,"jetE[numJet]/F");
  tree_->Branch("jetEMF",jetEMF,"jetEMF[numJet]/F");
  tree_->Branch("jetFHPD",jetFHPD,"jetFHPD[numJet]/F");
  tree_->Branch("jetN90Hits",jetN90Hits,"jetN90Hits[numJet]/F");

  if(usePAT_){
    tree_->Branch("numEle",&numEle,"numEle/I");
    tree_->Branch("eleEta",eleEta,"eleEta[numEle]/F");
    tree_->Branch("elePt",elePt,"elePt[numEle]/F");
    tree_->Branch("eleE",eleE,"eleE[numEle]/F");
    tree_->Branch("eleD0",eleD0,"eleD0[numEle]/F");
    tree_->Branch("eleD0Sig",eleD0Sig,"eleD0Sig[numEle]/F");
    tree_->Branch("eleHcalIso",eleHcalIso,"eleHcalIso[numEle]/F");
    tree_->Branch("eleEcalIso",eleEcalIso,"eleEcalIso[numEle]/F");
    tree_->Branch("eleTrkIso",eleTrkIso,"eleTrkIso[numEle]/F");
  }
  if(usePF_){
    pftree_= new TTree("pft","pft");
    pftree_->Branch("pfc",&pfccontainer_,"p/F:pt/F:dpt/F:dphi/F:e/F:et/F:phi/F:eta/F:ehad/F:eem/F:pid/I");
    pftree_->Branch("ev",&event_,"run/I:ev/I:lumi/I:eventclass/I");
  }
  
}


METNtupler::~METNtupler()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
METNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm; 

  histocontainer_["eventcount"]->Fill(0.0);
  

  event_.runnum = iEvent.id().run();
  event_.evnum = iEvent.id().event();
  event_.lumiblock = iEvent.luminosityBlock();
  event_.eventclass=0;

  numEle=0;
  numJet=0;
  if(usePAT_){
    edm::Handle<std::vector<pat::Electron> > eleHandle;
    iEvent.getByLabel("cleanPatElectrons",eleHandle);
    edm::Handle<std::vector<pat::Jet> > jetHandle;
    iEvent.getByLabel("cleanPatJets",jetHandle);
    edm::Handle<std::vector<pat::MET> > metHandle;
    iEvent.getByLabel("patMETs",metHandle);
    numEle=0;
    numJet=0;
    container_.met=metHandle->front().et();
    container_.phi=metHandle->front().phi();
    container_.sumEt=metHandle->front().sumEt();
    container_.sig=metHandle->front().significance();
    container_.prob = TMath::Prob(container_.sig,2);
    numEle=0;
    for(std::vector<pat::Electron>::const_iterator ele_iter = eleHandle->begin(); ele_iter!=eleHandle->end() && numEle<10; ++ele_iter){
      elePt[numEle]=ele_iter->pt();
      elePhi[numEle]=ele_iter->phi();
      eleE[numEle]=ele_iter->energy();
      eleEta[numEle]=ele_iter->eta();
      eleSCEta[numEle]=ele_iter->superCluster()->eta();
      eleHcalIso[numEle]=ele_iter->dr04HcalTowerSumEt();
      eleEcalIso[numEle]=ele_iter->dr04EcalRecHitSumEt();
      eleTrkIso[numEle]=ele_iter->dr04TkSumPt();
      eleD0[numEle]=ele_iter->dB();
      eleD0Sig[numEle]=0.0;//
      numEle++;
    }
    numJet=0;
    for(std::vector<pat::Jet>::const_iterator jet_iter = jetHandle->begin(); jet_iter!=jetHandle->end() && numJet<40; ++jet_iter){
      if(jet_iter->pt()<20)
	continue;
      if(fabs(jet_iter->eta())>3)
	continue;
      if(jet_iter->emEnergyFraction()<0.01)
	continue;
      if(jet_iter->jetID().fHPD>=0.98)
	continue;
      if(jet_iter->jetID().n90Hits<=1)
	continue;
      jetPt[numJet]=jet_iter->pt();
      jetEta[numJet]=jet_iter->eta();
      jetPhi[numJet]=jet_iter->phi();
      jetE[numJet]=jet_iter->energy();
      jetEMF[numJet]=jet_iter->emEnergyFraction();
      jetFHPD[numJet]=jet_iter->jetID().fHPD;
      jetN90Hits[numJet]=jet_iter->jetID().n90Hits;
      numJet++;
    }
  }
  else{
    edm::Handle<std::vector<reco::CaloJet> > jetHandle;
    iEvent.getByLabel(jetsrc_,jetHandle);
    numJet=0;
    for(std::vector<reco::CaloJet>::const_iterator jet_iter = jetHandle->begin(); jet_iter!=jetHandle->end() && numJet<40; ++jet_iter){
      const reco::CaloJet jet = *jet_iter;
      jetID->calculate(iEvent,jet);
      if(jet_iter->emEnergyFraction()<0.01)
	continue;
      if(jetID->n90Hits()<=1)
	continue;
      if(jetID->fHPD()>=0.98)
	continue;
      if(fabs(jet_iter->eta())>2.4)
	continue;
      jetE[numJet]=jet_iter->energy();
      jetPt[numJet]=jet_iter->pt();
      jetPhi[numJet]=jet_iter->phi();
      jetEta[numJet]=jet_iter->eta();
      jetEMF[numJet]=jet_iter->emEnergyFraction();
      jetFHPD[numJet]=jetID->fHPD();
      jetN90Hits[numJet]=jetID->n90Hits();
      numJet++;
    }// jet loop
  }
  if(!usePF_){
    Handle<std::vector<reco::CaloMET> > mets;
    iEvent.getByLabel(src_,mets);
        
    std::vector<reco::CaloMET>::const_iterator met=mets->begin();
    container_.met=met->et();
    container_.phi=met->phi();
    container_.sig=met->significance();
    container_.sumEt=met->sumEt();
    container_.prob = TMath::Prob(container_.sig,2);
    tree_->Fill();
  }//no PF
  if(usePF_){
    Handle<std::vector<reco::PFMET> > mets;
    iEvent.getByLabel(src_,mets);
        
    std::vector<reco::PFMET>::const_iterator met=mets->begin();
    container_.met=met->et();
    container_.phi=met->phi();
    container_.sig=met->significance();
    container_.sumEt=met->sumEt();
    container_.prob = TMath::Prob(container_.sig,2);
    
    Handle<PFCandidateCollection> pfcands;
    iEvent.getByLabel("particleFlow",pfcands);
    reco::PFCandidateCollection::const_iterator pfc;
    for(pfc=pfcands->begin(); pfc!=pfcands->end(); ++pfc){
      pfccontainer_.p=pfc->p();
      pfccontainer_.pt=pfc->pt();
      pfccontainer_.e=pfc->energy();
      pfccontainer_.et=pfc->et();
      pfccontainer_.phi=pfc->phi();
      pfccontainer_.eta=pfc->eta();
      pfccontainer_.ehad=pfc->hcalEnergy();
      pfccontainer_.eem=pfc->ecalEnergy();
      pfccontainer_.pid=pfc->particleId();
      pfccontainer_.dpt=-1;
      pfccontainer_.dphi=-1;
      if(pfc->trackRef().isAvailable()){
	pfccontainer_.dpt=pfc->trackRef()->ptError();
	pfccontainer_.dphi=pfc->trackRef()->phiError();
	//	    std::cout << "TRK : " << pfc->trackRef()->ptError() << " " << pfc->trackRef()->phiError() << std::endl;
      }
      else if(pfc->gsfTrackRef().isAvailable()){
	pfccontainer_.dpt= pfc->gsfTrackRef()->ptError();
	pfccontainer_.dphi=pfc->gsfTrackRef()->phiError();
      }
    }
    pftree_->Fill();
    tree_->Fill();
  }// usePF

}


// ------------ method called once each job just before starting event loop  ------------
void 
METNtupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METNtupler::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(METNtupler);
