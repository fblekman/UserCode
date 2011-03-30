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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TH1D.h>
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

      // ----------member data ---------------------------
  // input:
  //  jet info:
  std::string jetSrc_;
  std::string muSrc_;
  std::string trigSrc_;
  std::string trigname_;

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
  Float_t jet_tchp[20];
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
  jetSrc_=iConfig.getParameter<string>("jetSrc");
  muSrc_=iConfig.getParameter<string>("muSrc");
  trigSrc_=iConfig.getParameter<string>("trigSrc");
  trigname_=iConfig.getParameter<string>("trigName");
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
  //   using namespace edm;
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_,muons);


  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);

  nJet=0;
  for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
    if(jet->pt()>50)
      ++nJet;
  }



}


// ------------ method called once each job just before starting event loop  ------------
void 
BTagTriggerPlots::beginJob()
{
  histos["jetPt"] = fs->make<TH1D>("jetPt","jet p_{T}",100,0,400);
  histos["jetPtTrig"] = fs->make<TH1D>("jetPtTrig","jet p_{T}",100,0,400);

  tree = fs->make<TTree>("tree","tree");
  tree->Branch("nJet",&nJet,"nJet/I");
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
  tree->Branch("jet_tchp",jet_tchp,"jet_tchp[nJet]/F");


}

// ------------ method called once each job just after ending the event loop  ------------
void 
BTagTriggerPlots::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(BTagTriggerPlots);
