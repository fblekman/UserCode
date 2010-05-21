// -*- C++ -*-
//
// Package:    SelectGoodJetEvents
// Class:      SelectGoodJetEvents
// 
/**\class SelectGoodJetEvents SelectGoodJetEvents.cc Private/SelectGoodJetEvents/src/SelectGoodJetEvents.cc

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
#include "DataFormats/PatCandidates/interface/Jet.h"


//
// class declaration
//

class SelectGoodJetEvents : public edm::EDFilter {
   public:
      explicit SelectGoodJetEvents(const edm::ParameterSet&);
      ~SelectGoodJetEvents();

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
SelectGoodJetEvents::SelectGoodJetEvents(const edm::ParameterSet& iConfig):
  src_(iConfig.getParameter<edm::InputTag>("src"))
{
   //now do what ever initialization is needed

}


SelectGoodJetEvents::~SelectGoodJetEvents()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SelectGoodJetEvents::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  Handle<std::vector<pat::Jet> > handle_jets;
  iEvent.getByLabel(src_,handle_jets);
  size_t ngood=0;
  if(handle_jets->size()<2)
     return false;
  

  for(std::vector<pat::Jet>::const_iterator jet=handle_jets->begin();jet!=handle_jets->end(); ++jet){

     if(jet->pt()<20)
        continue;
     if(fabs(jet->eta())>3)
        continue;
     if(jet->jetID().n90Hits<=1)
        continue;
     if(jet->emEnergyFraction()<0.01)
        continue;
     if(jet->jetID().fHPD>0.98)
        continue;
//      std::cout << jet->eta() << " " << jet->pt() << " " << jet->corrStep() <<std::endl;
     //    for(jj=0; jj<jet->corrFactorSetLabels().size(); jj++)
// 	std::cout << jet->corrFactorSetLabels()[jj] << std::endl;
     ngood++;
     if(ngood>=2)
        return true;
  }
  // don't do anything with other events...
  if(ngood>=2)
     return true;
  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
SelectGoodJetEvents::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SelectGoodJetEvents::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SelectGoodJetEvents);
