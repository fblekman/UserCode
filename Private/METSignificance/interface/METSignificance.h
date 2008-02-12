// -*- C++ -*-
//
// Package:    METSignificance
// Class:      METSignificance
// 
/**\class METSignificance METSignificance.cc Wenu/METSignificance/src/METSignificance.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Freya Blekman
//         Created:  Mon Feb 11 10:31:24 CET 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TF1.h"
#include "Wenu/METSignificance/interface/significanceAlgo.h"
#include "Wenu/METSignificance/interface/SigInputObj.h"
//
// class decleration
//

class METSignificance : public edm::EDProducer {
   public:
      explicit METSignificance(const edm::ParameterSet&);
      ~METSignificance();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // not necessary  double calc_sigma_et(double par_n, double par_s, double par_c, double et);
      // ----------member data ---------------------------

  // method that does all checking 
  void fillVector(edm::Event&, const edm::EventSetup&);
  
  class uncertaintyFunctions{
  public:
    TF1 *etUncertainty;
    TF1 *phiUncertainty;
  };

  void setUncertaintyParameters();// fills the following uncertaintyFunctions objects:
  uncertaintyFunctions jetUncertainty;
  uncertaintyFunctions eleUncertainty;
  uncertaintyFunctions muonUncertainty;
  
  bool useJets_;
  bool useMuons_;
  bool useMuonDepositCor_;
  bool useElectrons_;
  
  float jetPtThreshold_;
  float jetEtaThreshold_;
  float muonPtThreshold_;
  float muonEtaThreshold_;
  float muonTrackD0Max_;
  float muonTrackDzMax_;
  int muonNHitsMin_;
  float muonDPtMax_;
  float muonChiSqMax_;

  
  
  edm::InputTag CaloJetAlgorithmTag_; 
  edm::InputTag CorJetAlgorithmTag_;
  std::string   JetCorrectionService_;
  edm::InputTag MuonTag_;
  edm::InputTag ElectronTag_;

  //  TrackDetectorAssociator   trackAssociator_; // for muon calo deposits
  //  TrackAssociatorParameters trackAssociatorParameters_;
  
  std::string significanceLabel_;
  std::vector<metsig::SigInputObj*> physobjvector_;
};
