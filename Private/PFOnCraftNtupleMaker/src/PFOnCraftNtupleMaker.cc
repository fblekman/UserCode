// -*- C++ -*-
//
// Package:    PFOnCraftNtupleMaker
// Class:      PFOnCraftNtupleMaker
// 
/**\class PFOnCraftNtupleMaker PFOnCraftNtupleMaker.cc Private/PFOnCraftNtupleMaker/src/PFOnCraftNtupleMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Freya Blekman
//         Created:  Tue Sep 22 13:26:47 CEST 2009
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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoParticleFlow/PFAnalyses/interface/EventDelegate.h"

#include "TTree.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TObject.h"
#include "TBranchElement.h"
#include <iostream>
//
// class decleration
//

class PFOnCraftNtupleMaker : public edm::EDAnalyzer {
   public:
      explicit PFOnCraftNtupleMaker(const edm::ParameterSet&);
      ~PFOnCraftNtupleMaker();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  edm::Service<TFileService> fs_;
  TTree *tree_;
  //  pftools::EventDelegate* ed_;
  
  
  class cosmictracks{
  public:
    int n;
    float p[50];
    float pt[50];
    float eta[50];
    float phi[50];
    float dz[50];
    float d0[50];
    float mutime[50];
    int   nhitsstrip[50];
    int   nhitspixel[50];
    int   nhitsmissed[50];
    float chi2[50];
    float chi2ndof[50];
  };
  cosmictracks tracks_;
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
PFOnCraftNtupleMaker::PFOnCraftNtupleMaker(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  tree_ = fs_->make<TTree>("pfinfo","");
  //  ed_->init(tree_, iConfig);
  TBranchElement* br = tree_->Branch("cos.",&tracks_);
}


PFOnCraftNtupleMaker::~PFOnCraftNtupleMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PFOnCraftNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   Handle<std::vector<reco::Track> > pTRK;
   iEvent.getByLabel("cosmictrackfinderP5",pTRK);
   const std::vector<reco::Track> & tracks = *pTRK;



//    ed_->startEvent(iEvent, iSetup);
//    ed_->processEvent(iEvent, iSetup);
   //and add the track info:
   // clear arrays:
   for(int i=0; i<50; ++i){
     tracks_.n=0;
     tracks_.p[i]=tracks_.pt[i]=tracks_.eta[i]=0;
     tracks_.phi[i]=tracks_.dz[i]=tracks_.d0[i]=tracks_.mutime[i]= tracks_.chi2[i]=-999;
     tracks_.chi2ndof[i]=tracks_.nhitsstrip[i]=tracks_.nhitspixel[i]=-1;
    
   }
   for(size_t itr=0; itr<tracks.size(); itr++){
     tracks_.p[tracks_.n]=tracks[itr].p();
     tracks_.pt[tracks_.n]=tracks[itr].pt();
     tracks_.eta[tracks_.n]=tracks[itr].eta();
     tracks_.d0[tracks_.n]=tracks[itr].d0();
     tracks_.dz[tracks_.n]=tracks[itr].dz();
     tracks_.chi2[tracks_.n]=tracks[itr].chi2();
     tracks_.chi2ndof[tracks_.n]=tracks[itr].ndof();
     // for pixel and strip hits one needs to loop over the hits:
     const reco::HitPattern& p = tracks[itr].hitPattern();
      
     tracks_.nhitspixel[tracks_.n]=p.numberOfValidPixelHits();
     tracks_.nhitsstrip[tracks_.n]=p.numberOfValidStripHits();
     tracks_.nhitsmissed[tracks_.n]=p.numberOfLostTrackerHits();
       
     tracks_.n++;
   }
   //   ed_->endEvent();
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
PFOnCraftNtupleMaker::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
PFOnCraftNtupleMaker::endJob() {
  //  ed_->finish();

}

//define this as a plug-in
DEFINE_FWK_MODULE(PFOnCraftNtupleMaker);
