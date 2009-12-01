// -*- C++ -*-
//
// Package:    PixelPlots
// Class:      PixelPlots
// 
/**\class PixelPlots PixelPlots.cc Private/PixelPlots/src/PixelPlots.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Freya Blekman
//         Created:  Tue Dec  1 14:06:25 CET 2009
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
//
// class decleration
//

class PixelPlots : public edm::EDAnalyzer {
   public:
      explicit PixelPlots(const edm::ParameterSet&);
      ~PixelPlots();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
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
PixelPlots::PixelPlots(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


PixelPlots::~PixelPlots()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PixelPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
PixelPlots::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixelPlots::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelPlots);
