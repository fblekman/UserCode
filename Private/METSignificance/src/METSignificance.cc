#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "CondFormats/JetMETObjects/interface/StandaloneJetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/METReco/interface/METCollection.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"

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
#include "FWCore/Framework/interface/ESHandle.h"

#include "Wenu/METSignificance/interface/METSignificance.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
METSignificance::METSignificance(const edm::ParameterSet& iConfig):
  // aLL THESE SHOULD BE SET VIA CFG FILE!!!
  useJets_(true),
  useMuons_(true),
  useMuonDepositCor_(false),
  useElectrons_(true),
  jetPtThreshold_(20.),
  jetEtaThreshold_(3.0),
  elePtThreshold_(10.),
  eleEtaThreshold_(3.0),
  muonPtThreshold_(10.),
  muonEtaThreshold_(2.5),
  muonTrackD0Max_(999.0),
  muonTrackDzMax_(999.0),
  muonNHitsMin_(5),
  muonDPtMax_(0.5),
  muonChiSqMax_(1000.0),
  significanceLabel_("ourMetSignificance"),
  CaloJetAlgorithmTag_( iConfig.getParameter<edm::InputTag>( "CaloJetAlgorithm" ) ),
  CorJetAlgorithmTag_( iConfig.getParameter<edm::InputTag>( "CorJetAlgorithm" ) ), 
  JetCorrectionService_( iConfig.getParameter<std::string>( "JetCorrectionService" ) ), 
  MuonTag_( iConfig.getParameter<edm::InputTag>("MuonTag") ),
  ElectronTag_( iConfig.getParameter<edm::InputTag>("ElectronTag"))
  
{

   //now do what ever other initialization is needed

  //edm::ParameterSet trackAssociatorParams = 
  //    iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  //  trackAssociatorParameters_.loadParameters(trackAssociatorParams);  
  //  trackAssociator_.useDefaultPropagator();   
  produces<reco::METCollection>() .setBranchAlias(significanceLabel_.c_str());
  produces<double>() .setBranchAlias(significanceLabel_.c_str());
}


METSignificance::~METSignificance()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  
}


void
METSignificance::fillVector(edm::Event& iEvent, const edm::EventSetup & iSetup)
{
  double tmp_sigma_e;
  double tmp_sigma_phi;
  // FB: Personally, I would split this up even more... but it's a start

  edm::Handle<reco::CaloJetCollection> caloJets;
  iEvent.getByLabel( CaloJetAlgorithmTag_, caloJets );
  //-- Set Jet Energy Corrections --//
  const JetCorrector* corrector = JetCorrector::getJetCorrector (JetCorrectionService_, iSetup);


  //Get the Muon collection
  edm::Handle<reco::MuonCollection> muCollection;
  iEvent.getByLabel(MuonTag_,muCollection);
  
  //--- Records Required for Muons -------------------------------//
  edm::ESHandle<MagneticField> magneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);

  edm::ESHandle<CaloGeometry> caloGeometryHandle;
  iSetup.get<IdealGeometryRecord>().get(caloGeometryHandle);

  // Muon Geometry
  edm::ESHandle<MuonDetLayerGeometry> muonDetLayerGeometryHandle;
  iSetup.get<MuonRecoGeometryRecord>().get(muonDetLayerGeometryHandle);
  
  edm::ESHandle<CSCGeometry> cscGeometryHandle;
  iSetup.get<MuonGeometryRecord>().get(cscGeometryHandle);

  edm::ESHandle<GlobalTrackingGeometry> globalTrackingGeometryHandle;
  iSetup.get<GlobalTrackingGeometryRecord>().get(globalTrackingGeometryHandle);

  edm::ESHandle<TrackerGeometry> trackerGeometryHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(trackerGeometryHandle);
  //--------------------------------------------------------------//

  edm::Handle<PixelMatchGsfElectronCollection> gsfElectrons;
 iEvent.getByLabel(ElectronTag_,gsfElectrons);
 
 // electron loop
 for (reco::PixelMatchGsfElectronCollection::const_iterator ele=gsfElectrons->begin(); ele!=gsfElectrons->end() && useElectrons_; ++ele){
  //-- Use min Pt cut of 20 Gev from MET note AN-2007/041 --//
    if (ele->pt() < elePtThreshold_) {
      continue; 
    }
    if(fabs(ele->eta())>eleEtaThreshold_)
      continue; 

    tmp_sigma_e = eleUncertainty.etUncertainty->Eval(ele->et());
    tmp_sigma_phi = eleUncertainty.phiUncertainty->Eval(ele->et());

    //-- Set up input for Significance Calculation --//
    metsig::SigInputObj* tmp_electron = new metsig::SigInputObj("electron",ele->et(),ele->phi(),tmp_sigma_e,tmp_sigma_phi);
    physobjvector_.push_back(tmp_electron);
  }
  
  //--- Calo Jets loop ------------------------------------------------//
  for( reco::CaloJetCollection::const_iterator cal = caloJets->begin(); 
       cal != caloJets->end() && useJets_; 
       ++cal ) {
    //-- Use min Pt cut of 20 Gev from MET note AN-2007/041 --//
    if (cal->pt() < jetPtThreshold_) {
      continue; 
    }
    if(fabs(cal->eta())>jetEtaThreshold_)
      continue;

    //-- Treat all calo jets with emEnergyFraction>0.9 as electrons. Cut value is from MET note AN-2007/041 --//
    if(cal->emEnergyFraction()>=0.9)
      continue;
   
    //-- Set Jet Energy Corrections --//
    double scale = corrector->correction(*cal);  // Jet energy correction
    
    //-- Apply Jet Energy Corrections --//
    double cor_et = scale*cal->et();  
    double tmp_phi = cal->phi();
    
    tmp_sigma_e = jetUncertainty.etUncertainty->Eval(cor_et);
    tmp_sigma_phi = jetUncertainty.phiUncertainty->Eval(cor_et);
    
    //-- Set up input for Significance Calculation --//
    metsig::SigInputObj* tmp_jet = new metsig::SigInputObj("jet",cor_et,tmp_phi,tmp_sigma_e,tmp_sigma_phi);
    physobjvector_.push_back(tmp_jet);
  }

  // muons:
  for( reco::MuonCollection::const_iterator muon = muCollection->begin();
       muon != muCollection->end() && useMuons_;
       ++muon)
    {
      //-- Impose cuts on muons --//
      if( muon->combinedMuon()->pt() > muonPtThreshold_ &&
	  fabs(muon->combinedMuon()->eta()) < muonEtaThreshold_ &&
	  fabs(muon->combinedMuon()->d0()) < muonTrackD0Max_ &&
	  fabs(muon->combinedMuon()->dz()) < muonTrackDzMax_ &&
	  muon->combinedMuon()->numberOfValidHits() > muonNHitsMin_ ) {
	
	float dpt_track = muon->combinedMuon()->error(0)/(muon->combinedMuon()->qoverp());
	float chisq = muon->combinedMuon()->normalizedChi2();

	if (dpt_track < muonDPtMax_ &&
	    chisq < muonChiSqMax_) {
	  
	  //-- Muon Energy --//
	  float muonEt = muon->et();
	  double muon_sigma_e = muonUncertainty.etUncertainty->Eval( muonEt);
	  double muon_sigma_phi = muonUncertainty.phiUncertainty->Eval( muonEt );

	  metsig::SigInputObj* tmp_muon = new metsig::SigInputObj("muon",muonEt,muon->phi(),muon_sigma_e,muon_sigma_phi);
	  physobjvector_.push_back(tmp_muon);
	  
	  //-- Calculate muon energy deposition in the calorimeters --//

	  // FB: I could not get this to work. I think it needs newer code than CMSSW167.
	  
// 	  if (useMuonDepositCor_) {
// 	    reco::TrackRef mu_track = muon->combinedMuon();
// 	    TrackDetMatchInfo info = 
// 	      trackAssociator_.associate(iEvent, iSetup,
// 					 trackAssociator_.getFreeTrajectoryState(iSetup, *mu_track),
// 					 trackAssociatorParameters_);
// 	    double ene = info.crossedEnergy(TrackDetMatchInfo::TowerTotal);//total calotower energy
// 	    std::cout << "muon energy = " << ene << std::endl;  //DEBUGGING
// 	    float muonCaloEtDepX = ene*sin((*mu_track).theta())*cos((*mu_track).phi());
// 	    float muonCaloEtDepY = ene*sin((*mu_track).theta())*sin((*mu_track).phi());
// 	    float muonCaloEtDepSum = ene*sin((*mu_track).theta());		  
// 	    
	    
// 	  }
	  
	}
      }
    }
  //--- End Muons -----------------------------------------------------//

}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
METSignificance::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // call the fillVector method. This is where all the various objects are selected:
  fillVector(iEvent,iSetup);

  //-- Calculate the MET Significance --//
  double calo_met_total=0;
  double calo_met_phi=0;
  double calo_met_set=0;
  double significance = metsig::ASignificance(physobjvector_, calo_met_total, calo_met_phi, calo_met_set);

  double met_x = calo_met_total * cos(calo_met_phi);
  double met_y = calo_met_total * sin(calo_met_phi);

  // now book a missing et object:
  if(significance>0){
    std::cout << "strange significance " << significance << " " << calo_met_total << " " <<  calo_met_phi << " " <<  calo_met_set << " " << met_x << " " << met_y << std::endl;
  }  
  math::XYZTLorentzVector p4( met_x, met_y, 0.0, calo_met_total);
  math::XYZPoint vtx(0,0,0);
  reco::MET met( calo_met_set, p4, vtx );
  std::auto_ptr<METCollection> metcoll;
  metcoll.reset(new METCollection);
  metcoll->push_back( met );
  iEvent.put( metcoll );
  std::auto_ptr<double> themetsig(new double(significance));
  iEvent.put( themetsig );

  physobjvector_.clear();

}
//----------------------
void
METSignificance::setUncertaintyParameters(){
 // set the various functions here:

  jetUncertainty.etUncertainty = new TF1("jetEtFunc","x*sqrt(([2]*[2])+([1]*[1]/x)+([0]*[0]/(x*x)))",3);
  // [0]= par_n, [1]=par_s, [2]= par_c
  // values from PTDR 1, ch 11.4
  jetUncertainty.etUncertainty->SetParameter(0,1.25);
  jetUncertainty.etUncertainty->SetParameter(1,5.6);
  jetUncertainty.etUncertainty->SetParameter(2,0.033);
  // phi value from our own fits:
  jetUncertainty.phiUncertainty = new TF1("jetPhiFunc","[0]*x",1);
  jetUncertainty.phiUncertainty->SetParameter(0,2.635*(3.14159/180.));
  
  // completely ambiguious values for electron-like jets...
  // the egamma group keeps track of these here:
  // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCMSSWVal
  eleUncertainty.etUncertainty = new TF1("eleEtFunc","[0] * x",1);
  eleUncertainty.etUncertainty->SetParameter(0,0.034); // electron resolution in energy is around 3.4%, measured for 10 < pT < 50 at realistic events with pile-up.
  eleUncertainty.phiUncertainty = new TF1("elePhiFunc","[0] * x",1);
  eleUncertainty.phiUncertainty->SetParameter(0,1*(3.14159/180.));

  // and ambiguious values for the muons...
  muonUncertainty.etUncertainty = new TF1("muonEtFunc","[0] * x",1);
  muonUncertainty.etUncertainty->SetParameter(0,0.01);
  muonUncertainty.phiUncertainty = new TF1("muonPhiFunc","[0] * x",1);
  muonUncertainty.phiUncertainty->SetParameter(0,1*(3.14159/180.));
 
}
// ------------ method called once each job just before starting event loop  ------------
void 
METSignificance::beginJob(const edm::EventSetup&)
{
  setUncertaintyParameters();
 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METSignificance::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(METSignificance);
