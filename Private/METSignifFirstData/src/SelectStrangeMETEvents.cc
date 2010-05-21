// -*- C++ -*-
//
// Package:    SelectStrangeMETEvents
// Class:      SelectStrangeMETEvents
// 
/**\class SelectStrangeMETEvents SelectStrangeMETEvents.cc Private/SelectStrangeMETEvents/src/SelectStrangeMETEvents.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Freya Blekman,40 4-B11,+41227671552,
//         Created:  Mon Apr 19 13:13:20 CEST 2010
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/METReco/interface/CaloMET.h"


//
// class declaration
//

class SelectStrangeMETEvents : public edm::EDFilter {
   public:
      explicit SelectStrangeMETEvents(const edm::ParameterSet&);
      ~SelectStrangeMETEvents();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  
      // ----------member data ---------------------------
  edm::InputTag src_;
  double metmin_;
  double signifmin_;

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
SelectStrangeMETEvents::SelectStrangeMETEvents(const edm::ParameterSet& iConfig):
  src_(iConfig.getParameter<edm::InputTag>("src")), 
  metmin_(iConfig.getUntrackedParameter<double>("MinMetToSave",30)),
  signifmin_(iConfig.getUntrackedParameter<double>("MinSignifToSave",10))
{
   //now do what ever initialization is needed

}


SelectStrangeMETEvents::~SelectStrangeMETEvents()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SelectStrangeMETEvents::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  Handle<std::vector<reco::CaloMET> > mets;
  iEvent.getByLabel(src_,mets);

  std::vector<reco::CaloMET>::const_iterator met=mets->begin();
  if(met->et()>metmin_)
    return true;
  if(met->significance()>signifmin_)
    return true;
  // don't do anything with other events...
  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
SelectStrangeMETEvents::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SelectStrangeMETEvents::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SelectStrangeMETEvents);
