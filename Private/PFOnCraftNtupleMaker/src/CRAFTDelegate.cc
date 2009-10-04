#include "Private/PFOnCraftNtupleMaker/interface/CRAFTDelegate.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowReco/interface/PFSimParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFSimParticleFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFNuclearInteraction.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"

#include "RecoParticleFlow/Benchmark/interface/PFBenchmarkAlgo.h"

#include <iostream>

#include <cmath>

using namespace pftools;
using namespace reco;
using namespace edm;

CRAFTDelegate::CRAFTDelegate() {
  LogDebug("CRAFTDelegate") << __PRETTY_FUNCTION__ << std::endl;
}

void CRAFTDelegate::initCore(const edm::ParameterSet& parameters) {

  inputTagCandidates_ = parameters.getParameter<InputTag> ("PFCandidates");
  inputTagRecTracks_ = parameters.getParameter<InputTag> ("PFRecTracks");
  inputTagClustersEcal_
    = parameters.getParameter<InputTag> ("PFClustersEcal");
  inputTagClustersHcal_
    = parameters.getParameter<InputTag> ("PFClustersHcal");
  inputTagRecHitsEcal_ = parameters.getParameter<InputTag> ("PFRecHitsEcal");
  inputTagRecHitsHcal_ = parameters.getParameter<InputTag> ("PFRecHitsHcal");

  inputTagRawHitsEcalEB_ = parameters.getParameter<InputTag> (
							      "RawRecHitsEcalEB");
  inputTagRawHitsEcalEE_ = parameters.getParameter<InputTag> (
							      "RawRecHitsEcalEE");
  inputTagRawHitsHcal_ = parameters.getParameter<InputTag> ("RawRecHitsHcal");

}

bool CRAFTDelegate::processEvent(const edm::Event& event,
				 const edm::EventSetup& setup) {
  LogDebug("CRAFTDelegate") << __PRETTY_FUNCTION__ << "\n";
  //start with connecting sim particles with pf candidates, if simulation

  if (debug_ > 1)
    LogDebug("CRAFTDelegate") << "\tProceeding in MC mode..."
			      << std::endl;
	
  std::vector<Track> cosmictracks = **cosmicTracks_;
  PFCandidateCollection candidates = **pfCandidates_;
  PFRecTrackCollection tracks = **recTracks_;
  PFClusterCollection clustersEcal = **clustersEcal_;
  PFClusterCollection clustersHcal = **clustersHcal_;

  thisEventPasses_=true;
  return thisEventPasses_;
}

void CRAFTDelegate::startEventCore(const edm::Event& event,
				   const edm::EventSetup& setup) {

  pfCandidates_ = new Handle<PFCandidateCollection> ;
  cosmicTracks_ = new Handle<TrackCollection>;
  clustersEcal_ = new Handle<PFClusterCollection> ;
  clustersHcal_ = new Handle<PFClusterCollection> ;
  recHitsEcal_ = new Handle<PFRecHitCollection> ;
  recHitsHcal_ = new Handle<PFRecHitCollection> ;
  recTracks_ = new Handle<PFRecTrackCollection> ;

  rawRecHitsEcalEB_ = new Handle<EcalRecHitCollection> ;
  rawRecHitsEcalEE_ = new Handle<EcalRecHitCollection> ;
  rawRecHitsHcal_ = new Handle<HBHERecHitCollection> ;

  getCollection(*pfCandidates_, inputTagCandidates_, event);
  getCollection(*clustersEcal_, inputTagClustersEcal_, event);
  getCollection(*clustersHcal_, inputTagClustersHcal_, event);
  getCollection(*recHitsEcal_, inputTagRecHitsEcal_, event);
  getCollection(*recHitsHcal_, inputTagRecHitsHcal_, event);
  getCollection(*recTracks_, inputTagRecTracks_, event);


}

void CRAFTDelegate::startParticleCore() {
}

void CRAFTDelegate::endParticleCore() {

}

//Checks vetos
//Return true if you want the particle written to the tree
bool CRAFTDelegate::endEventCore() {
  if (thisEventPasses_) {
    ++nWrites_;
  } else {
    ++nFails_;
  }
  delete pfCandidates_;
  delete clustersEcal_;
  delete clustersHcal_;
  delete recHitsEcal_;
  delete recHitsHcal_;
  delete recTracks_;


  return true;
}

