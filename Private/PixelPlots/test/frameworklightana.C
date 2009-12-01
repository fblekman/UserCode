#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <vector>
#include <string>
#include <TMath.h>
#include "DataFormats/Math/interface/deltaR.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
// INCLUDE ALL CMSSW FORMATS HERE!

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#endif

void frameworklightana(void){

  std::vector<std::string> fileNames;
  fileNames.push_back("reco_RAW2DIGI_RECO.root");
  fwlite::ChainEvent ev(fileNames);

  fwlite::Handle<edm::DetSetVector<PixelDigi> > h_digis;
  edm::DetSetVector<PixelDigi> pixelDigis;

  TH2D *pixelhistory = new TH2D("pixelhistory","",100,0,100,10,0,10);
  //  pixelhistory->SetTitleX("number of hit pixels per detector");
  //  pixelhistory->SetTitleY("number of events next hit in same detid");
  
  std::map<uint64_t,uint64_t> hitmap;
  size_t iev=0;
  size_t ii,jj,kk;

  for( ev.toBegin();
       ! ev.atEnd();
       ++ev) {
    iev++;
    if(iev%10==0)
      std::cout << iev << std::endl;
    bool firedtrigger=false; // triggers are tough in fwlite!
    
    h_digis.getByLabel(ev,"siPixelDigis");
    pixelDigis=*h_digis;
    edm::DetSetVector<PixelDigi>::const_iterator digiIter;
    //    std::cout << pixelDigis.size() << std::endl;
    for(digiIter=pixelDigis.begin(); digiIter!=pixelDigis.end(); ++digiIter){// ITERATOR OVER DET IDs
      uint32_t detid = digiIter->id;
      edm::DetSet<PixelDigi>::const_iterator ipix; // ITERATOR OVER DIGI DATA  
      
      for(ipix = digiIter->data.begin(); ipix!=digiIter->end(); ++ipix){
	std::cout << detid << " " << ipix->adc() << " " << ipix->row() << " " << ipix->column() << std::endl;
	    
      }
    }
    
  }// end of event loop

  TCanvas *canv  = new TCanvas();
 
}// end of method
