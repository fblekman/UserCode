// -*- C++ -*-
//
// Package:    PATMissingHTSignificanceProducer
// Class:      PATMissingHTSignificanceProducer
// 
/**\class PATMissingHTSignificanceProducer PATMissingHTSignificanceProducer.cc myana/PATMissingHTSignificanceProducer/src/PATMissingHTSignificanceProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Freya Blekman
//         Created:  Tue Apr 22 20:14:24 CEST 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"

#include "RecoMET/METAlgorithms/interface/SigInputObj.h"
#include "RecoMET/METAlgorithms/interface/significanceAlgo.h"

#include <sstream>
#include <iostream>
//
// class decleration
//

class PATMissingHTSignificanceProducer : public edm::EDProducer {
   public:
      explicit PATMissingHTSignificanceProducer(const edm::ParameterSet&);
      ~PATMissingHTSignificanceProducer();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data --------------------------- 
  edm::InputTag eleLabel_;
  edm::InputTag muoLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag tauLabel_;
  edm::InputTag metLabel_;
  edm::InputTag phoLabel_;
  std::string significanceLabel_;
  bool includeElectrons_;
  bool includeMuons_;
  bool includeJets_;
  bool includePhotons_;
  bool includeTaus_;
  bool guessEmptyResolutions_;
  bool includeJES_;
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
PATMissingHTSignificanceProducer::PATMissingHTSignificanceProducer(const edm::ParameterSet& iConfig) : 
  significanceLabel_("PATMissingHTSignificance"),  
  eleLabel_(iConfig.getParameter<edm::InputTag>("electronTag")),
  muoLabel_(iConfig.getParameter<edm::InputTag>("muonTag")),
  jetLabel_(iConfig.getParameter<edm::InputTag>("jetTag")),
  tauLabel_(iConfig.getParameter<edm::InputTag>("tauTag")),
  metLabel_(iConfig.getParameter<edm::InputTag>("metTag")),
  phoLabel_(iConfig.getParameter<edm::InputTag>("photonTag")),
  includeElectrons_(iConfig.getParameter<bool>("useElectrons")),
  includeMuons_(iConfig.getParameter<bool>("useMuons")),
  includeJets_(iConfig.getParameter<bool>("useJets")),
  includePhotons_(iConfig.getParameter<bool>("usePhotons")),
  includeTaus_(iConfig.getParameter<bool>("useTaus")),
  guessEmptyResolutions_(iConfig.getParameter<bool>("fillInResolutionWhenEmpty")),
  includeJES_(iConfig.getParameter<bool>("includeJESuncertaintyForJets"))
{
   //register your products

  produces<double>() .setBranchAlias(significanceLabel_.c_str());
  produces<reco::METCollection>() .setBranchAlias(significanceLabel_.c_str());
   //now do what ever other initialization is needed

}


PATMissingHTSignificanceProducer::~PATMissingHTSignificanceProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PATMissingHTSignificanceProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
   edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByLabel(muoLabel_,muonHandle);
   std::vector<pat::Muon> muons = *muonHandle;
   
   edm::Handle<std::vector<pat::Jet> > jetHandle;
   iEvent.getByLabel(jetLabel_,jetHandle);
   std::vector<pat::Jet> jets = *jetHandle;

   edm::Handle<std::vector<pat::Electron> > electronHandle;
   iEvent.getByLabel(eleLabel_,electronHandle);
   std::vector<pat::Electron> electrons = *electronHandle;

   edm::Handle<std::vector<pat::MET> > metHandle;
   iEvent.getByLabel(metLabel_,metHandle);
   std::vector<pat::MET> mets = *metHandle;

   edm::Handle<std::vector<pat::Photon> > phoHandle;
   iEvent.getByLabel(phoLabel_,phoHandle);
   std::vector<pat::Photon> photons = *phoHandle;
   
   edm::Handle<std::vector<pat::Tau> > tauHandle;
   iEvent.getByLabel(tauLabel_,tauHandle);
   std::vector<pat::Tau> taus = *tauHandle;

   // this is the worker class for the MET significance.
   std::vector<metsig::SigInputObj> signInputVec;
   double phi_val,phi_err,et_val,et_err;
   std::string type="unknown";
   std::ostringstream infomsg;
   // loop over jets.
   for(std::vector<pat::Jet>::const_iterator jet_iter = jets.begin(); jet_iter!=jets.end() && includeJets_; ++jet_iter){
     phi_err = jet_iter->resolutionPhi();
     et_err=jet_iter->resolutionEt();
     phi_val=jet_iter->phi();
     et_val=jet_iter->et();
     type="jet";
     double jes_var=0;
     if(includeJES_){ // check JES uncertainty:
       jes_var = fabs(jet_iter->jetCorrFactors().scaleDefault() - jet_iter->jetCorrFactors().scaleUds()) ;
       jes_var += fabs(fabs(jet_iter->jetCorrFactors().scaleDefault() - jet_iter->jetCorrFactors().scaleGlu()));
       jes_var /= 2.;
       jes_var *=et_val;
       //       std::cout << jes_var << " " << jes_var* et_val<<  " " << et_val << " " << et_err << " " ; 
       //       std::cout << et_err << std::endl;
     }
     et_err = sqrt(pow(et_err,2)+pow(jes_var,2));
     infomsg  << "MHTSig input > \t"<< type <<",\t pT:" <<  et_val << "+/-" << et_err << " (includes dJES " <<jes_var<< "), phi:" << phi_val << "+/-" << phi_err << " eta:" << jet_iter->eta()  << std::endl;
     metsig::SigInputObj tempjet(type,et_val,phi_val,et_err,phi_err);
     signInputVec.push_back(tempjet);
     
   }
   // loop over electrons.
   for(std::vector<pat::Electron>::const_iterator electron_iter = electrons.begin(); electron_iter!=electrons.end() && includeElectrons_; ++electron_iter){
     
     phi_err = electron_iter->resolutionPhi();
     et_err=electron_iter->resolutionEt();
     phi_val=electron_iter->phi();
     et_val=electron_iter->et();
     type="electron";
     infomsg  << "MHTSig input > \t"<< type <<",\t pT:" <<  et_val << "+/-" << et_err << ", phi:" << phi_val << "+/-" << phi_err << " eta:" << electron_iter->eta()<< std::endl;
     metsig::SigInputObj tempelectron(type,et_val,phi_val,et_err,phi_err);
     signInputVec.push_back(tempelectron);
   }
   // loop over photons.
   for(std::vector<pat::Photon>::const_iterator photon_iter = photons.begin(); photon_iter!=photons.end() && includePhotons_; ++photon_iter){
     phi_err = photon_iter->resolutionPhi();
     et_err=photon_iter->resolutionEt();
     phi_val=photon_iter->phi();
     et_val=photon_iter->et();
     type="photon";
     if(phi_err == 0 && guessEmptyResolutions_)
       phi_err = 0.01*et_val;
     if(et_err ==0 && guessEmptyResolutions_)
       et_err = 0.2 * et_val;
     infomsg  << "MHTSig input > \t"<< type <<",\t pT:" <<  et_val << "+/-" << et_err << ", phi:" << phi_val << "+/-" << phi_err << " eta:" << photon_iter->eta()<< std::endl;
     metsig::SigInputObj tempphoton(type,et_val,phi_val,et_err,phi_err);
     signInputVec.push_back(tempphoton);
   }
   // loop over muons.
   for(std::vector<pat::Muon>::const_iterator muon_iter = muons.begin(); muon_iter!=muons.end() && includeMuons_; ++muon_iter){
     phi_err = muon_iter->resolutionPhi();
     et_err=muon_iter->resolutionEt();
     phi_val=muon_iter->phi();
     et_val=muon_iter->et();
     type="muon";
     infomsg  << "MHTSig input > \t"<< type <<",\t pT:" <<  et_val << "+/-" << et_err << ", phi:" << phi_val << "+/-" << phi_err << " eta:" << muon_iter->eta()<< std::endl;
     metsig::SigInputObj tempmuon(type,et_val,phi_val,et_err,phi_err); 
     signInputVec.push_back(tempmuon);
   }  
   // loop over taus.
   for(std::vector<pat::Tau>::const_iterator tau_iter = taus.begin(); tau_iter!=taus.end() && includeTaus_; ++tau_iter){
     phi_err = tau_iter->resolutionPhi();
     et_err=tau_iter->resolutionEt();
     phi_val=tau_iter->phi();
     et_val=tau_iter->et();
     type="tau";
     if(phi_err == 0 && guessEmptyResolutions_)
       phi_err = 0.01*et_val;
     if(et_err ==0 && guessEmptyResolutions_)
       et_err = 0.2 * et_val;
     infomsg  << "MHTSig input > \t"<< type <<",\t pT:" <<  et_val << "+/-" << et_err << ", phi:" << phi_val << "+/-" << phi_err << " eta:" << tau_iter->eta()<< std::endl;
     metsig::SigInputObj temptau(type,et_val,phi_val,et_err,phi_err); 
     signInputVec.push_back(temptau);
   }
   
   // now calculate the object:
   double sign_scalar_et=0;
   double sign_met_phi=0;
   double sign_met=0;
   double significance = metsig::ASignificance(signInputVec,sign_met,sign_met_phi,sign_scalar_et);
   infomsg << "MHTSig output > \t" << "now done with HT Calculation: mHT:" << sign_met << ",HT:" << sign_scalar_et << ", phi:" << sign_met_phi << ", mHTSignificance:" << significance << std::endl;

   edm::LogInfo("PATMissingHTSignificanceProducer::produce()") << infomsg.str() ;

   math::XYZTLorentzVector p4( sign_met * cos(sign_met_phi), sign_met * sin(sign_met_phi), 0.0, sign_met);
   math::XYZPoint vtx(0,0,0);
   reco::MET met( sign_scalar_et, p4, vtx );
   std::auto_ptr<METCollection> metcoll;
   metcoll.reset(new METCollection);
   metcoll->push_back( met );
   iEvent.put( metcoll );
   std::auto_ptr<double> themetsig(new double(significance));
   iEvent.put( themetsig );

   // make sure worker vector is empty.
   signInputVec.clear();
}

// ------------ method called once each job just before starting event loop  ------------
void 
PATMissingHTSignificanceProducer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATMissingHTSignificanceProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATMissingHTSignificanceProducer);
