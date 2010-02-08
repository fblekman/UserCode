// -*- C++ -*-
//
// Package:    METSignifFirstData
// Class:      METSignifFirstData
// 
/**\class METSignifFirstData METSignifFirstData.cc Private/METSignifFirstData/src/METSignifFirstData.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Freya Blekman,6 R-029,+41227678914,
//         Created:  Tue Jan 26 17:00:23 CET 2010
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

#include "FWCore/ParameterSet/interface/ParameterSjÄe "FWCore/ParameterSet/interfaceöÇÇà⁄[ò€YHëï–€‹ôK‘Ÿ\ùöXŸistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TTree.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include <vector>
//
// class declaration
//

class METSignifFirstData : public edm::EDAnalyzer {
   public:
      explicit METSignifFirstData(const edm::ParameterSet&);
      ~METSignifFirstData();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

  std::vector<uint64_t> goodruns_;
  std::vector<uint32_t> techTrigs_;

  std::map<std::string,TH1D*> histos_;
  TTree *tree;

  struct evtstr{
    Int_t run;
    Int_t evt;
  };
  evtstr tree_evt;

  struct metstr{
    Float_t met;
    ™vt;

  struct metstr{
∆ˆE˜B6ñvÊñc∞¢f∆ˆ†t;

  struct mt_t sumEt;
    Float_t px;
    Float_t py;
    Float_t phi;
    
  };
  metstr tree_metstr;
  
  Int_t tree_nct;
  Float_t tree_ee[5000];
  Float_t tree_R4[5000];
  Float_t tree_eb[5000];
  Float_t tree_hb[5000];
  Float_t tree_he[5000];
  Float_t tree_ho[5000];
  Float_t tree_hf[5000];
  Float_t tree_alpha[5000];
  
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
METSignifFirstData::METSignifFirstData(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  // dummy run for testing:
  goodruns_.push_back(124115);
  // got these from Artur
  goodruns_.push_back(123815);
  goodruns_.push_back(123818);
  goodruns_.push_back( 123906);
  goodruns_.push_back( 123908);
  goodruns_.push_back( 123909);
  goodruns_.push_back( 123970);
  goodruns_.push_back( 123976);
  goodruns_.push_back( 123977);
  goodruns_.push_back( 123978);
  goodruns_.push_back(123985);
  goodruns_.push_back( 123987);
  goodruns_.push_back( 124006);
  goodruns_.push_back( 124008);
  goodruns_.push_back( 124009);
  goodruns_.push_back( 124017);
  goodruns_.push_back( 124020);
  goodruns_.push_back( 124022);
  goodruns_.push_back( 124023);
  goodruns_.push_back(   124024);
  goodruns_.push_back( 124025);
  goodruns_.push_back( 124027);
  goodruns_.push_back( 124030);

  edm::Service<TFileService> fs;
  TH1D *temp = fs->make<TH1D>("empty","empty",1,0,1);
  //tree = new TTree("tree","tree");
  //tree->Branch("met",&tree_met,"met/F");
  //  tree->Branch("metsignif",&tree_signif,"metsignif/F");
  tree=fs->make<TTree>("tree","tree");
  tree->Branch("met",&tree_metstr,"met/F:signif/F:sumEt/F:px/F:py/F:phi/F");
  tree->Branch("evt",&tree_evt,"run/I:evt/I");

  //  tree->Branch("nct",&tree_nct,"nct/I");
  
}


METSignifFirstData::~METSignifFirstData()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
METSignifFirstData::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   if( std::find(goodruns_.begin(),goodruns_.end(),iEvent.id().run())==goodruns_.end() && iEvent.id().run()>100){
     std::cout << "run " << iEvent.id().run() << " is not in good run list" << std::endl;
     return;
   }

   //   std::cout << "now looking at run: " << iEvent.id().run() << ", event " << iEvent.id().event() << std::endl;


   edm::Handle<L1GlobalTriggerReadoutRecord> gtReadoutRecord;
   iEvent.getByLabel( edm::InputTag("gtDigis"), gtReadoutRecord);

   //   const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = gtReadoutRecord->technicalTriggerWord();
      
   //   for (int ttr = 0; ttr< 50; ttr++) {
   //     std::cout << technicalTriggerWordBeforeMask.at(ttr) << " ";
   //   }
   //   std::cout <<std::endl;

   // require bit 0 to have fired.
   if(technicalTriggerWordBeforeMask.at(0)==0)
     return;

   // print out the other relevant bits:
   //   std::cout << "bit 40: " << technicalTriggerWordBeforeMask.at(40) << std::endl;
   //   std::cout << "bit 41: " << technicalTriggerWordBeforeMask.at(41) << std::endl;
   //   edm::Handle<edm::View<PFCandidate> > pfCandidates;
   //   iEvent.getByLabel("ParticleFlow" pfCandidates);
   


   tree_evt.run=iEvent.id().run();
   tree_evt.evt= iEvent.id().event();

   edm::Handle<reco::CaloMETCollection > metHandle;
   iEvent.getByLabel("met",metHandle);
   const reco::CaloMETCollection *metColl = metHandle.product();

   if(metColl->size()<=0)
     return;
   const reco::CaloMET met = metColl->front();
  //Fill the variables

      
   //   std::cout << met.pt() << " " << met.phi() << " " << met.sumEt() << std::endl;
   tree_metstr.met=met.pt();
   tree_metstr.signif=met.metSignificance();
   tree_metstr.phi = met.phi();
   tree_metstr.sumEt= met.sumEt();
   tree_metstr.px=met.px();
   tree_metstr.py=met.py();
   tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
METSignifFirstData::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METSignifFirstData::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(METSignifFirstData);
