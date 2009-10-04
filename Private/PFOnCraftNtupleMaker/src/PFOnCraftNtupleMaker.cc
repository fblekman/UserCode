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
// $Id: PFOnCraftNtupleMaker.cc,v 1.1 2009/09/22 13:48:15 fblekman Exp $
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
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoParticleFlow/PFAnalyses/interface/EventDelegate.h"
#include "RecoParticleFlow/PFAnalyses/interface/TestbeamDelegate.h"
#include "RecoParticleFlow/PFAnalyses/interface/DipionDelegate.h"
#include "Private/PFOnCraftNtupleMaker/interface/CRAFTDelegate.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"

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
  pftools::EventDelegate* ed_;
  TrackAssociatorParameters trackasspars_;
  const MagneticField* magfield_;

  
  int n_ct;
  float p_ct[50];
  float pt_ct[50];
  float eta_ct[50];
  float phi_ct[50];
  float etatop_ct[50];
  float phitop_ct[50];
  float etadown_ct[50];
  float phidown_ct[50];
  float dz_ct[50];
  float d0_ct[50];
  float mutime_ct[50];
  int   nhitsstrip_ct[50];
  int   nhitspixel_ct[50];
  int   nhitsmissed_ct[50];
  float chi2_ct[50];
  float chi2ndof_ct[50];
  
  int n_pfc;
  float e_pfc[200];
  float p_pfc[200];
  float momGSF_pfc[200];
  float eta_pfc[200];
  float phi_pfc[200];
  float pt_pfc[200];
  int type_pfc[200];
  float energyEcal_pfc[200];
  float energyHcal_pfc[200];
  float dRclosest_pfc[200];
  float dRclosestCal_pfc[200];

    
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
 
  

  edm::ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  trackasspars_.loadParameters( parameters );

  tree_ = fs_->make<TTree>("pfinfo","");
  //  ed_ = new pftools::CRAFTDelegate();
  //  ed_->init(tree_, iConfig);
  tree_->Branch("n_ct",&n_ct,"n_ct/I");
  tree_->Branch("p_ct",p_ct,"p_ct[n_ct]/F");
  tree_->Branch("pt_ct",pt_ct,"pt_ct[n_ct]/F");
  tree_->Branch("eta_ct",eta_ct,"eta_ct[n_ct]/F");
  tree_->Branch("phi_ct",phi_ct,"phi_ct[n_ct]/F");
  tree_->Branch("etatop_ct",etatop_ct,"etatop_ct[n_ct]/F");
  tree_->Branch("phitop_ct",phitop_ct,"phitop_ct[n_ct]/F");
  tree_->Branch("etadown_ct",etadown_ct,"etadown_ct[n_ct]/F");
  tree_->Branch("phidown_ct",phidown_ct,"phidown_ct[n_ct]/F");
  tree_->Branch("dz_ct",dz_ct,"dz_ct[n_ct]/F");
  tree_->Branch("d0_ct",d0_ct,"d0_ct[n_ct]/F");
  tree_->Branch("mutime_ct",mutime_ct,"mutime_ct[n_ct]/F");
  tree_->Branch("chi2_ct",chi2_ct,"chi2_ct[n_ct]/F");
  tree_->Branch("chi2ndof_ct",chi2ndof_ct,"chi2ndof_ct[n_ct]/F");
  tree_->Branch("nhitsstrip_ct",nhitsstrip_ct,"nhitsstrip_ct[n_ct]/I");
  tree_->Branch("nhitspixel_ct",nhitspixel_ct,"nhitspixel_ct[n_ct]/I");
  tree_->Branch("nhitsmissed_ct",nhitsmissed_ct,"nhitsmissed_ct[n_ct]/I");

  tree_->Branch("n_pfc",&n_pfc,"n_pfc/I");
  tree_->Branch("e_pfc",e_pfc,"e_pfc[n_pfc]/F");
  tree_->Branch("p_pfc",p_pfc,"p_pfc[n_pfc]/F");
  tree_->Branch("pt_pfc",pt_pfc,"pt_pfc[n_pfc]/F");
  tree_->Branch("eta_pfc",eta_pfc,"eta_pfc[n_pfc]/F");
  tree_->Branch("phi_pfc",phi_pfc,"phi_pfc[n_pfc]/F");
  tree_->Branch("type_pfc",type_pfc,"type_pfc[n_pfc]/I");
  tree_->Branch("energyEcal_pfc",energyEcal_pfc,"energyEcal_pfc[n_pfc]/F");
  tree_->Branch("energyHcal_pfc",energyHcal_pfc,"energyHcal_pfc[n_pfc]/F");
  tree_->Branch("momGSF_pfc",momGSF_pfc,"momGSF_pfc[n_pfc]/F");
  tree_->Branch("dRclosest_pfc",dRclosest_pfc,"dRclosest_pfc[n_pfc]/F");
  tree_->Branch("dRclosestCal_pfc",dRclosestCal_pfc,"dRclosestCal_pfc[n_pfc]/F");
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
   // get tracking infos:
   edm::ESHandle<MagneticField> esmagfield;
   iSetup.get<IdealMagneticFieldRecord>().get(esmagfield);
   magfield_=&(*esmagfield);

   // create propagators for track propagation to calorimeters, as this is a cosmic track this needs to be done for up and down from center of detector:
  // these are deleted at the end of the analysis loop because they're auto-pointers

   Propagator *opposite = new SteppingHelixPropagator(magfield_,oppositeToMomentum);
   Propagator *along = new SteppingHelixPropagator(magfield_,alongMomentum);
   TrackDetectorAssociator *tda_up = new TrackDetectorAssociator();
   TrackDetectorAssociator *tda_down = new TrackDetectorAssociator();
   tda_up->useDefaultPropagator();
   tda_down->useDefaultPropagator();
   //   tda_up->setPropagator(opposite);
   //   tda_down->setPropagator(along);
 
   Handle<std::vector<reco::Track> > pTRK;
   iEvent.getByLabel("cosmictrackfinderP5",pTRK);
   const std::vector<reco::Track> & tracks = *pTRK;

   bool save_event=false;
   //   ed_->startEvent(iEvent, iSetup);
   //   ed_->processEvent(iEvent, iSetup);
   //and add the track info:
   // clear arrays:
   for(size_t i=0; i<50; ++i){
     n_ct=0;
     p_ct[i]=pt_ct[i]=eta_ct[i]=0;
     phi_ct[i]=dz_ct[i]=d0_ct[i]=mutime_ct[i]= chi2_ct[i]=-999;
     chi2ndof_ct[i]=nhitsstrip_ct[i]=nhitspixel_ct[i]=-1;
    
   }
   for(size_t i=0; i<200; ++i){
     n_pfc=0;
     type_pfc[i]=-1;
     phi_pfc[i]=eta_pfc[i]=-999;
     energyEcal_pfc[i]= energyEcal_pfc[i]= e_pfc[i]=pt_pfc[i]=-1;
   }
   for(size_t itr=0; itr<tracks.size(); itr++){
     if(n_ct>=50)
       continue;
     p_ct[n_ct]=tracks[itr].p();
     pt_ct[n_ct]=tracks[itr].pt();
     eta_ct[n_ct]=tracks[itr].eta();
     phi_ct[n_ct]=tracks[itr].phi();
     d0_ct[n_ct]=tracks[itr].d0();
     dz_ct[n_ct]=tracks[itr].dz();
     chi2_ct[n_ct]=tracks[itr].chi2();
     chi2ndof_ct[n_ct]=tracks[itr].ndof();
     // for pixel and strip hits one needs to loop over the hits:
     const reco::HitPattern& p = tracks[itr].hitPattern();
      
     nhitspixel_ct[n_ct]=p.numberOfValidPixelHits();
     nhitsstrip_ct[n_ct]=p.numberOfValidStripHits();
     nhitsmissed_ct[n_ct]=p.numberOfLostTrackerHits();
       
     // do track - calorimeter association. Important: these tracks have two entries as the completely cross the tracking volume. This makes our life more difficult.
     // first: get a free trajectory ('the track'):
     
     FreeTrajectoryState startingstate_up = tda_up->getFreeTrajectoryState(iSetup,tracks[itr]);
     FreeTrajectoryState startingstate_down = tda_down->getFreeTrajectoryState(iSetup,tracks[itr]);
     // associate this state two ways, 'up' and 'down':
     TrackDetMatchInfo upstate=tda_up->associate(iEvent,iSetup,tracks[itr],trackasspars_,TrackDetectorAssociator::OutsideIn);
     TrackDetMatchInfo downstate=tda_down->associate(iEvent,iSetup,tracks[itr],trackasspars_,TrackDetectorAssociator::InsideOut);
     etatop_ct[n_ct]=upstate.trkGlobPosAtEcal.eta();
     phitop_ct[n_ct]=upstate.trkGlobPosAtEcal.phi();
     etadown_ct[n_ct]=downstate.trkGlobPosAtEcal.eta();
     phidown_ct[n_ct]=downstate.trkGlobPosAtEcal.phi();


     n_ct++;
   }
   // loop over PF Candidates:
   if(n_ct>0)
     save_event=true;
   
   Handle<std::vector<reco::PFCandidate> > pPFC;
   iEvent.getByLabel("particleFlow",pPFC);
   const std::vector<reco::PFCandidate> & pfcs = *pPFC;

   for(size_t ipfc=0; ipfc<pfcs.size(); ++ipfc){  
     if(n_pfc>=200)
       continue;
     e_pfc[n_pfc]=pfcs[ipfc].energy();
     phi_pfc[n_pfc]=pfcs[ipfc].phi();
     eta_pfc[n_pfc]=pfcs[ipfc].eta();
     type_pfc[n_pfc]=pfcs[ipfc].particleId();
     pt_pfc[n_pfc]=pfcs[ipfc].pt();
     energyEcal_pfc[n_pfc]=pfcs[ipfc].ecalEnergy();
     energyEcal_pfc[n_pfc]=pfcs[ipfc].hcalEnergy();
     
     // calculate deltaR with respect to entry and exit points.
     dRclosest_pfc[n_pfc]=dRclosestCal_pfc[n_pfc]=1000;
     
     for(int itr=0;itr<n_ct; ++itr){
       if(dRclosest_pfc[n_pfc]>deltaR(eta_pfc[n_pfc],phi_pfc[n_pfc],eta_ct[itr],phi_ct[itr]))
	 dRclosest_pfc[n_pfc]=deltaR(eta_pfc[n_pfc],phi_pfc[n_pfc],eta_ct[itr],phi_ct[itr]);
       if(dRclosestCal_pfc[n_pfc]>deltaR(eta_pfc[n_pfc],phi_pfc[n_pfc],etatop_ct[itr],phitop_ct[itr]))
	 dRclosestCal_pfc[n_pfc]=deltaR(eta_pfc[n_pfc],phi_pfc[n_pfc],etatop_ct[itr],phitop_ct[itr]);
       if(dRclosestCal_pfc[n_pfc]>deltaR(eta_pfc[n_pfc],phi_pfc[n_pfc],etadown_ct[itr],phidown_ct[itr]))
	 dRclosestCal_pfc[n_pfc]=deltaR(eta_pfc[n_pfc],phi_pfc[n_pfc],etadown_ct[itr],phidown_ct[itr]);
     }
     n_pfc++;
   }
   if(save_event)
     tree_->Fill();
   //   ed_->endEvent();
  
   // delete pointers:
   delete opposite;
   delete along;
   delete tda_up;
   delete tda_down;
 
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
  //  tree_->Write(); // not necessary, TFileService does that for you
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFOnCraftNtupleMaker);
