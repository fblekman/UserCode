// -*- C++ -*-
//
// Package:    ExampleAnalysis
// Class:      ExampleAnalysis
// 
/**\class ExampleAnalysis ExampleAnalysis.cc ExampleAnalysis/ExampleAnalysis/src/ExampleAnalysis.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Dr Freya Blekman
//         Created:  Wed Sep 20 10:08:38 BST 2006
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

// necessary objects:
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"
#include "DataFormats/TrackReco/interface/Track.h"

// basic cpp:
#include <iostream>
// some root includes
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
//
// class decleration
//

class ExampleAnalysis : public edm::EDAnalyzer {
   public:
      explicit ExampleAnalysis(const edm::ParameterSet&);
      ~ExampleAnalysis();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  int nEvt;// used to count the number of events

  
  // to be used for root output tree
  TFile *thefile;
  TTree *smalltree;
  int nMCPar;
  int MCPARMAX;// used to set maximum of arrays
  double MCPar_px[20];
  double MCPar_py[20];
  double MCPar_pz[20];
  double MCPar_e[20];
  int MCPar_pdgid[20];
  int nEleCand;
  int ELEMAX;// used to set maximum of arrays
  double EleCand_px[10];
  double EleCand_py[10];
  double EleCand_pz[10];
  double EleCand_e[10];
  int nMuonCand;
  int MUONMAX;// used to set maximum of arrays
  double MuonCand_px[10];
  double MuonCand_py[10];
  double MuonCand_pz[10];
  double MuonCand_p[10];
  int MuonCand_charge[10];
  int nJetCand;
  int JETMAX;// used to set maximum of arrays
  double JetCand_px[30];
  double JetCand_py[30];
  double JetCand_pz[30];
  double JetCand_e[30];
  // a histogram
  TH1D *zmassmuons;
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
ExampleAnalysis::ExampleAnalysis(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

  // set the maximal values of the arrays. This is a clutch but needs to be done to make sure your arrays don't run out of bounds 
  nEvt=0;
  MCPARMAX=20;
  ELEMAX=10;
  MUONMAX=10;
  JETMAX=30;
  
  thefile = new TFile("ExampleSmallAnalysisOutput.root","recreate");
  thefile->cd();
  zmassmuons = new TH1D("zmassmuons","invariant mass of opposite sign #mu",100,0,200);
  smalltree= new TTree("ExampleTree","ExampleTree");
  smalltree->Branch("nMCPar",&nMCPar,"nMCPar/I");
  smalltree->Branch("MCPar_px",MCPar_px,"MCPar_px[nMCPar]/D");
  smalltree->Branch("MCPar_py",MCPar_py,"MCPar_py[nMCPar]/D");
  smalltree->Branch("MCPar_pz",MCPar_pz,"MCPar_pz[nMCPar]/D");
  smalltree->Branch("MCPar_e",MCPar_e,"MCPar_e[nMCPar]/D");
  smalltree->Branch("MCPar_pdgid",MCPar_pdgid,"MCPar_pdgid[nMCPar]/I");
  smalltree->Branch("nEleCand",&nEleCand,"nEleCand/I");
  smalltree->Branch("EleCand_px",EleCand_px,"EleCand_px[nEleCand]/D");
  smalltree->Branch("EleCand_py",EleCand_py,"EleCand_py[nEleCand]/D");
  smalltree->Branch("EleCand_pz",EleCand_pz,"EleCand_pz[nEleCand]/D");
  smalltree->Branch("EleCand_e",EleCand_e,"EleCand_e[nEleCand]/D");
  smalltree->Branch("nJetCand",&nJetCand,"nJetCand/I");
  smalltree->Branch("JetCand_px",JetCand_px,"JetCand_px[nJetCand]/D");
  smalltree->Branch("JetCand_py",JetCand_py,"JetCand_py[nJetCand]/D");
  smalltree->Branch("JetCand_pz",JetCand_pz,"JetCand_pz[nJetCand]/D");
  smalltree->Branch("JetCand_e",JetCand_e,"JetCand_e[nJetCand]/D");
  smalltree->Branch("nMuonCand",&nMuonCand,"nMuonCand/I");
  smalltree->Branch("MuonCand_px",MuonCand_px,"MuonCand_px[nMuonCand]/D");
  smalltree->Branch("MuonCand_py",MuonCand_py,"MuonCand_py[nMuonCand]/D");
  smalltree->Branch("MuonCand_pz",MuonCand_pz,"MuonCand_pz[nMuonCand]/D");
  smalltree->Branch("MuonCand_p",MuonCand_p,"MuonCand_p[nMuonCand]/D");
}


ExampleAnalysis::~ExampleAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ExampleAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // this analyzer produces a small root file with basic candidates and some MC information
   // some additional print statements
   nEvt++;
   if((nEvt%10==0 && nEvt<=100)||(nEvt%100==0 && nEvt>100))
     std::cout << "reading event " << nEvt << std::endl;
   

   // step 1: fill some basic MC information into the root tree
   Handle<HepMCProduct> evt;
   iEvent.getByLabel("VtxSmeared", evt);
   
   nMCPar=0;
   HepMC::GenEvent * myGenEvent = new  HepMC::GenEvent(*(evt->GetEvent()));
   for ( HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin();
	 p != myGenEvent->particles_end() && nMCPar<MCPARMAX ; ++p ) {
     if ( abs((*p)->pdg_id()) !=0 && abs((*p)->pdg_id())<30 && nMCPar<MCPARMAX){
       MCPar_px[nMCPar]=(*p)->momentum().x();
       MCPar_py[nMCPar]=(*p)->momentum().y();
       MCPar_pz[nMCPar]=(*p)->momentum().z();
       MCPar_e[nMCPar]=(*p)->momentum().mag();//(*p)->energy() does not work!!;
       MCPar_pdgid[nMCPar]=(*p)->pdg_id();
       nMCPar++;
     }
   }
   // step 2: get some reconstructed objects... and fill them into the root tree too!
   //2a: electrons
   Handle<reco::ElectronCollection> electrons;
   iEvent.getByLabel("pixelMatchElectrons",electrons);
   reco::ElectronCollection::const_iterator ele;
   nEleCand=0;
   for( ele = electrons->begin(); ele != electrons->end() && nEleCand<ELEMAX; ++ ele ) {
     EleCand_e[nEleCand]=ele->energy();
     EleCand_px[nEleCand]=ele->px();
     EleCand_py[nEleCand]=ele->py();
     EleCand_pz[nEleCand]=ele->pz();
     nEleCand++;
   }
   //2b: muons
   
   Handle<reco::TrackCollection> muons;
   iEvent.getByLabel("globalMuons",muons);
   reco::TrackCollection::const_iterator muon;
   nMuonCand=0;
   for( muon = muons->begin(); muon != muons->end() && nMuonCand<MUONMAX; ++ muon ) {
     MuonCand_p[nMuonCand]=muon->p();
     MuonCand_px[nMuonCand]=muon->px();
     MuonCand_py[nMuonCand]=muon->py();
     MuonCand_pz[nMuonCand]=muon->pz();
     MuonCand_charge[nMuonCand]=muon->charge();
     nMuonCand++;
   }
   //2c: jets
   Handle<reco::CaloJetCollection> jets;
   iEvent.getByLabel("iterativeCone5CaloJets",jets);
   reco::CaloJetCollection::const_iterator jet;
   nJetCand=0;
   for( jet = jets->begin(); jet != jets->end() && nJetCand<JETMAX; ++ jet ) {
     JetCand_e[nJetCand]=jet->energy();
     JetCand_px[nJetCand]=jet->px();
     JetCand_py[nJetCand]=jet->py();
     JetCand_pz[nJetCand]=jet->pz();
     nJetCand++;
   }

   
   // calculate something and fill a histogram:
   for(int imu=0; imu<nMuonCand; imu++){
     for(int jmu=imu+1; jmu<nMuonCand; jmu++){
       if(MuonCand_charge[imu]*MuonCand_charge[jmu]<0){//opposite charge muons
	 double mass = pow(MuonCand_p[imu]+MuonCand_p[jmu],2);
	 mass-=pow(MuonCand_px[imu]+MuonCand_px[jmu],2);
	 mass-=pow(MuonCand_py[imu]+MuonCand_py[jmu],2);
	 mass-=pow(MuonCand_pz[imu]+MuonCand_pz[jmu],2);
	 zmassmuons->Fill(sqrt(mass));
       }
     }
   }
   // and last action in the event: save the event in the small tree:
   smalltree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
ExampleAnalysis::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExampleAnalysis::endJob() {
  std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "analyzed " << nEvt << " events " << std::endl;
  std::cout << "writing information into file: " << thefile->GetName() << std::endl;
  thefile->Write();
  thefile->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExampleAnalysis)
