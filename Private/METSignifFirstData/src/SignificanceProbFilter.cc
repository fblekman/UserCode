// -*- C++ -*-
//
// Package:    SignificanceProbFilter
// Class:      SignificanceProbFilter
// 
/**\class SignificanceProbFilter SignificanceProbFilter.cc Private/SignificanceProbFilter/src/SignificanceProbFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Freya Blekman,40 4-B11,+41227671552,
//         Created:  Mon Apr 19 13:13:20 CEST 2010
// $Id: SignificanceProbFilter.cc,v 1.1 2010/05/21 16:58:31 fblekman Exp $
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
#include "DataFormats/METReco/interface/PFMET.h"
#include "TMath.h"
//
// class declaration
//

class SignificanceProbFilter : public edm::EDFilter {
   public:
      explicit SignificanceProbFilter(const edm::ParameterSet&);
      ~SignificanceProbFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  
      // ----------member data ---------------------------
      edm::InputTag src_;
      bool useCalo_;
      double metmin_;
      double signifmin_;
      double probmin_;

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
SignificanceProbFilter::SignificanceProbFilter(const edm::ParameterSet& iConfig):
  src_(iConfig.getParameter<edm::InputTag>("src")), 
  useCalo_(iConfig.getParameter<bool>("isCaloMET")),
  metmin_(iConfig.getParameter<double>("MinMetToSave")),
  signifmin_(iConfig.getParameter<double>("MinSignifToSave")),
  probmin_(iConfig.getParameter<double>("MaxProbabilityToSave"))
{
   //now do what ever initialization is needed

}


SignificanceProbFilter::~SignificanceProbFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SignificanceProbFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  if(!useCalo_){
     Handle<std::vector<reco::PFMET> > mets;
     iEvent.getByLabel(src_,mets);
     
     std::vector<reco::PFMET>::const_iterator met=mets->begin();
     if(met->et()>metmin_)
        return true;
     if(met->significance()>signifmin_)
        return true;
     if(TMath::Prob(met->significance(),2)<probmin_)
        return true;
     // don't do anything with other events...
  }
  else{
     Handle<std::vector<reco::CaloMET> > mets;
     iEvent.getByLabel(src_,mets);
     
     std::vector<reco::CaloMET>::const_iterator met=mets->begin();
     if(met->et()>metmin_)
        return true;
     if(met->significance()>signifmin_)
        return true;
     if(TMath::Prob(met->significance(),2)<probmin_)
        return true;
  }
     return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
SignificanceProbFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SignificanceProbFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SignificanceProbFilter);
