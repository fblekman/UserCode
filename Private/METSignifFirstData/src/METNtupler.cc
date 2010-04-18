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

#include "DataFormats/METReco/interface/MET.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class declaration
//

class METNtupler : public edm::EDProducer {
   public:
      explicit METNtupler(const edm::ParameterSet&);
      ~METNtupler();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag src_;
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
  src_(iConfig.getParameter<edm::InputTag>("src"))
{

  produces<std::vector<float> >("et").setBranchAlias("et");
  produces<std::vector<float> >("sumEt").setBranchAlias("sumEt");
  produces<std::vector<float> >("phi").setBranchAlias("phi");
  produces<std::vector<float> >("significance").setBranchAlias("significance");
  
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
METNtupler::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm; using namespace std;
   Handle<std::vector<reco::MET> > mets;
   iEvent.getByLabel(src_,mets);
   auto_ptr<vector<float> > metPhi( new vector<float>);
   auto_ptr<vector<float> > metEt( new vector<float>);
   auto_ptr<vector<float> > metSumEt(new vector<float>);
   auto_ptr<vector<float> > metSignif(new vector<float>);
   const reco::MET &met = (*mets)[0];
   metPhi->push_back(met.phi());
   metEt->push_back(met.et());
   metSumEt->push_back(met.sumEt());
   metSignif->push_back(met.significance());
   
   iEvent.put(metPhi,"phi");
   iEvent.put(metEt,"et");
   iEvent.put(metSumEt,"sumEt");
   iEvent.put(metSignif,"significance");
   
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
