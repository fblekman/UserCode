// -*- C++ -*-
//
// Package:    PixelEventFilter
// Class:      PixelEventFilter
// 
/**\class PixelEventFilter PixelEventFilter.cc Private/PixelEventFilter/src/PixelEventFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Freya Blekman
//         Created:  Tue Dec  1 15:16:14 CET 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class PixelEventFilter : public edm::EDFilter {
   public:
      explicit PixelEventFilter(const edm::ParameterSet&);
      ~PixelEventFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  uint32_t pixelcounter;
  uint32_t maxevents;
  uint32_t digithreshold;
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
PixelEventFilter::PixelEventFilter(const edm::ParameterSet& iConfig):
  pixelcounter(0),
  maxevents(iConfig.getUntrackedParameter<uint32_t>("saveEventTail",10)),
  digithreshold(iConfig.getUntrackedParameter<uint32_t>("minDigisToReset",20))
{
   //now do what ever initialization is needed

}


PixelEventFilter::~PixelEventFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
PixelEventFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< edm::DetSetVector<PixelDigi> > pixelDigis;
  iEvent.getByLabel( "siPixelDigis", pixelDigis );
  
    // loop over the data and reset counter if digis are found
  pixelcounter++;
  uint32_t ndigis=0;
  edm::DetSetVector<PixelDigi>::const_iterator digiIter;
  for(digiIter=pixelDigis->begin(); digiIter!=pixelDigis->end() && ndigis<=digithreshold; ++digiIter){// ITERATOR OVER DET IDs
    uint32_t detid = digiIter->id;
    edm::DetSet<PixelDigi>::const_iterator ipix; // ITERATOR OVER DIGI DATA  

    for(ipix = digiIter->data.begin(); ipix!=digiIter->end() && ndigis<=digithreshold; ++ipix){
      ndigis++;
      if(ndigis>=digithreshold)
	pixelcounter=0;
    }
  }
  

  if(pixelcounter<=maxevents)
    return true;
  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
PixelEventFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixelEventFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelEventFilter);
