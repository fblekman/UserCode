// -*- C++ -*-
//
// Package:    CombinedNtupler
// Class:      CombinedNtupler
// 
/**\class CombinedNtupler CombinedNtupler.cc EWK/CombinedNtupler/src/CombinedNtupler.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Aleko Khukhunaishvili,6 R-029,+41227678914,
//         Created:  Thu Apr 15 12:00:05 CEST 2010
// $Id: CombinedNtupler.cc,v 1.1 2010/06/03 16:56:21 fblekman Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "RecoJets/JetAlgorithms/interface/JetIDHelper.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/PatCandidates/interface/MHT.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
//#include "FWCore/Framework/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"


#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TMath.h>

using namespace std;
using namespace edm;

struct GenInfo {
  int pid;
  float pthat;
  float alphaQCD;
  float alphaQED;
  float scalePDF;
  int id1;
  int id2;
  float x1;
  float x2;
  float xPDF1;
  float xPDF2;
};

struct RecoElectron{
  int charge;
  int  isEB;
  int  isEE;
  float p;
  float pt;
  float eta;
  float phi;
  float dr04TkSumPt;
  float dr04EcalRecHitSumEt;
  float dr04HcalTowerSumEt;
  float hadronicOverEm;
  float scSigmaEtaEta;
  float deltaEtaSuperClusterTrackAtVtx;
  float deltaPhiSuperClusterTrackAtVtx;
  float eSuperClusterOverP;
};

struct MET{
  float et;
  float phi;
  float sig;
  float sumEt;
  float emf;
  float prob;
};

struct MHT{
  int njets;
  float et;
  float phi;
  float ht;
  float sig;
  float rawmet;
  float rawmetsig;
  float jecmet;
  float jecmetsig;
};

struct Jet{
  float e;
  float et;
  float eta;
  float phi;
  float dpt;
  float dphi;
  float emf;
  float l2l3;
  float fHPD;
  int n90hits;
};

const int nElectrons=2;
const int nJets=2;

//
// class declaration
//

class CombinedNtupler : public edm::EDFilter {
public:
  explicit CombinedNtupler(const edm::ParameterSet&);
  ~CombinedNtupler();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------

  bool	    isMC_;
  string	    OutputFileName_;

  edm::InputTag   TriggerResultsTag_;
  string	    TriggerPath_;
  TriggerNames    triggerNames_;
    

  edm::InputTag   genmetTag_;
  edm::InputTag   calometTag_;
  edm::InputTag   pfmetTag_;
  edm::InputTag   tcmetTag_;
  edm::InputTag   calometoldTag_;
  edm::InputTag   pfmetoldTag_;

  edm::InputTag   t1met05Tag_;
  edm::InputTag   t1met10Tag_;
  edm::InputTag   t1met20Tag_;
  edm::InputTag   t2met05Tag_;
  edm::InputTag   t2met10Tag_;
  edm::InputTag   t2met20Tag_;


  edm::InputTag   mht05Tag_;
  edm::InputTag   mht10Tag_;
  edm::InputTag   mht20Tag_;
  edm::InputTag   mht05emfTag_;
  edm::InputTag   mht10emfTag_;
  edm::InputTag   mht20emfTag_;
  edm::InputTag   l2l3mht05Tag_;
  edm::InputTag   l2l3mht10Tag_;
  edm::InputTag   l2l3mht20Tag_;
  edm::InputTag   l2l3mht05newTag_;
  edm::InputTag   l2l3mht10newTag_;
  edm::InputTag   l2l3mht20newTag_;
  edm::InputTag   l2l3mht05emfTag_;
  edm::InputTag   l2l3mht10emfTag_;
  edm::InputTag   l2l3mht20emfTag_;
  edm::InputTag   jptmht05emfTag_;
  edm::InputTag   jptmht10emfTag_;
  edm::InputTag   jptmht20emfTag_;

  edm::InputTag   electronTag_;
  edm::InputTag     patElectronTag_;

  edm::InputTag   jetTag_;
  edm::InputTag patJetTag_;
  edm::InputTag   jetIDTag_;
  string	    jetCorrectorTag_;

  edm::InputTag triggerTag_;


  edm::InputTag   genjetTag_;



  //    TFile		*OutFile__file;
  TTree		*results_tree;	 
  TH1D *eventcount;// keeps track of processed events...

  //variables to save in the ntuple
  // trigger info, is converted into human-friendly format when written to tree...
  std::vector<std::string> triggernames_;/// gets read from event configuration. Uses whatever is in triggerTag_ (src for triggers)
  std::vector<std::string> hltnamesworker_; // worker vector.
  std::map<std::string,int> triggervals_; // values of triggers, one-to-on mapping with triggernames_.
  // other:
  GenInfo geninfo;
  int run, event, lumi;
  bool hlt;
  MET genmet, calomet, pfmet, tcmet;
  MET calometold, pfmetold;
  MET t1met05, t1met10, t1met20;
  MET t2met05, t2met10, t2met20;
  MHT mht05, mht10, mht20,
    mht05emf, mht10emf, mht20emf,
    l2l3mht05, l2l3mht10, l2l3mht20, 
    l2l3mht05new, l2l3mht10new, l2l3mht20new, 
    l2l3mht05emf, l2l3mht10emf, l2l3mht20emf,
    jptmht05emf, jptmht10emf, jptmht20emf;
  RecoElectron e[2];
  Jet calojet[2];
  Jet patjet[2];
  Jet genjet[2];

  bool filter_;

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
CombinedNtupler::CombinedNtupler(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  filter_		= iConfig.getUntrackedParameter<bool>("filter");
  isMC_		= iConfig.getParameter<bool>("isMC");

  //trigger
  TriggerResultsTag_	= iConfig.getUntrackedParameter<edm::InputTag>("TriggerResultsTag");
  TriggerPath_	= iConfig.getUntrackedParameter<string>("TriggerPath");
  triggernames_ = iConfig.getUntrackedParameter<std::vector<std::string> >("TriggerNamesForNtuple");
  for(size_t ii=0; ii<triggernames_.size(); ii++){
    triggervals_[triggernames_[ii]]=-5;
  }
  //met's
  genmetTag_		= iConfig.getParameter<edm::InputTag>("genmetTag");
  calometTag_		= iConfig.getParameter<edm::InputTag>("calometTag");
  pfmetTag_		= iConfig.getParameter<edm::InputTag>("pfmetTag");
  tcmetTag_		= iConfig.getParameter<edm::InputTag>("tcmetTag");
  calometoldTag_		= iConfig.getUntrackedParameter<edm::InputTag>("calometoldTag");
  pfmetoldTag_		= iConfig.getUntrackedParameter<edm::InputTag>("pfmetoldTag");

  t1met05Tag_		= iConfig.getUntrackedParameter<edm::InputTag>("t1met05Tag");
  t1met10Tag_		= iConfig.getUntrackedParameter<edm::InputTag>("t1met10Tag");
  t1met20Tag_		= iConfig.getUntrackedParameter<edm::InputTag>("t1met20Tag");
  t2met05Tag_		= iConfig.getUntrackedParameter<edm::InputTag>("t2met05Tag");
  t2met10Tag_		= iConfig.getUntrackedParameter<edm::InputTag>("t2met10Tag");
  t2met20Tag_		= iConfig.getUntrackedParameter<edm::InputTag>("t2met20Tag");


  //mht's
  mht05Tag_		= iConfig.getUntrackedParameter<edm::InputTag>("mht05Tag");
  mht10Tag_		= iConfig.getUntrackedParameter<edm::InputTag>("mht10Tag");
  mht20Tag_		= iConfig.getUntrackedParameter<edm::InputTag>("mht20Tag");
  mht05emfTag_	= iConfig.getUntrackedParameter<edm::InputTag>("mht05emfTag");
  mht10emfTag_	= iConfig.getUntrackedParameter<edm::InputTag>("mht10emfTag");
  mht20emfTag_	= iConfig.getUntrackedParameter<edm::InputTag>("mht20emfTag");
  l2l3mht05Tag_	= iConfig.getUntrackedParameter<edm::InputTag>("l2l3mht05Tag");
  l2l3mht10Tag_	= iConfig.getUntrackedParameter<edm::InputTag>("l2l3mht10Tag");
  l2l3mht20Tag_	= iConfig.getUntrackedParameter<edm::InputTag>("l2l3mht20Tag");
  l2l3mht05newTag_	= iConfig.getUntrackedParameter<edm::InputTag>("l2l3mht05newTag");
  l2l3mht10newTag_	= iConfig.getUntrackedParameter<edm::InputTag>("l2l3mht10newTag");
  l2l3mht20newTag_	= iConfig.getUntrackedParameter<edm::InputTag>("l2l3mht20newTag");
  l2l3mht05emfTag_	= iConfig.getUntrackedParameter<edm::InputTag>("l2l3mht05emfTag");
  l2l3mht10emfTag_	= iConfig.getUntrackedParameter<edm::InputTag>("l2l3mht10emfTag");
  l2l3mht20emfTag_	= iConfig.getUntrackedParameter<edm::InputTag>("l2l3mht20emfTag");
  jptmht05emfTag_	= iConfig.getUntrackedParameter<edm::InputTag>("jptmht05emfTag");
  jptmht10emfTag_	= iConfig.getUntrackedParameter<edm::InputTag>("jptmht10emfTag");
  jptmht20emfTag_	= iConfig.getUntrackedParameter<edm::InputTag>("jptmht20emfTag");

  //electrons
  electronTag_	= iConfig.getUntrackedParameter<edm::InputTag>("electronTag");
  patElectronTag_     = iConfig.getParameter<edm::InputTag>("patElectronTag");
  //calojets
  jetTag_		= iConfig.getParameter<edm::InputTag>("jetTag");
  patJetTag_          = iConfig.getParameter<edm::InputTag>("patJetTag");
  jetIDTag_		= iConfig.getParameter<edm::InputTag>("jetIDTag");
  //  jetCorrectorTag_	= iConfig.getUntrackedParameter<string>("jetCorrectorTag");

  genjetTag_		= iConfig.getParameter<edm::InputTag>("genjetTag");
}

CombinedNtupler::~CombinedNtupler(){;}


// ------------ method called on each new Event  ------------
bool
CombinedNtupler::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  eventcount->Fill(0);// this histogram counts the events.

  run = iEvent.id().run();
  event = iEvent.id().event();
  lumi = iEvent.id().luminosityBlock();

  //some initialization to identify if there was such object or not.
  for(int i=0; i<nElectrons; i++) e[i].charge=0;
  for(int i=0; i<nJets; i++) {calojet[i].e = -1.; genjet[i].e = -1.;patjet[i].e=-1;}

  if(isMC_){
    Handle<GenEventInfoProduct> gi;
    iEvent.getByLabel("generator", gi);

    geninfo.pid = (int)gi->signalProcessID();
    geninfo.pthat = (float)gi->qScale();
    geninfo.alphaQCD = (float)gi->alphaQCD();
    geninfo.alphaQED = (float)gi->alphaQED();
    if(gi->hasPDF()){
      geninfo.scalePDF = gi->pdf()->scalePDF;
      geninfo.id1 = gi->pdf()->id.first;
      geninfo.id2 = gi->pdf()->id.second;
      geninfo.x1 = gi->pdf()->x.first;
      geninfo.x2 = gi->pdf()->x.second;
      geninfo.xPDF1 = gi->pdf()->xPDF.first;
      geninfo.xPDF2 = gi->pdf()->xPDF.second;
    }
  }

    
  //
  
  edm::Handle<TriggerResults> hltResults;
  iEvent.getByLabel(TriggerResultsTag_, hltResults);
  if(hltResults.product()->wasrun()){
    edm::TriggerNames triggerNames;
    triggerNames.init( *hltResults );

    hltnamesworker_=triggerNames.triggerNames();

    for(int itrig=0; itrig< (int)hltnamesworker_.size();++itrig){ 
      const bool accept(hltResults->accept(itrig));

      int trigbit=0;
      if(accept)
	trigbit=1;
      if(!hltResults->wasrun(itrig))
	trigbit=-1;
      if(hltResults->error(itrig))
	trigbit=-2;

      for(size_t ii=0; ii<triggernames_.size(); ii++){
	if(hltnamesworker_[itrig]==triggernames_[ii])
	  triggervals_[triggernames_[ii]]=trigbit;
	
      }
    }      
  }

    
  edm::Handle<edm::View<reco::GsfElectron> > electronHandle;
  iEvent.getByLabel(electronTag_, electronHandle);
  edm::View<reco::GsfElectron> electrons = *electronHandle;
  edm::View<reco::GsfElectron>::const_iterator ele;

  //get the two highets pt electrons
  //those are already sorted. Well, does not hurt...

  if(filter_ && electrons.size()==0) return false;

  unsigned int idx[2];
  idx[0]=999999; idx[1]=999999;
  float pmax1=0, pmax2=0;
  for (unsigned int i=0; i<electrons.size(); ++i){
    float pt=electrons[i].pt();
    if(pt>pmax1 && pt>pmax2){
      idx[0]=i;
      pmax1= pt;
    }
    else if (pt<pmax1 && pt>pmax2){
      idx[1] = i;
      pmax2 = pt;
    }
  }

  for(int i=0; i<2; i++){
    if (idx[i] == 999999) continue;
    reco::GsfElectron electron = electrons[idx[i]];
    e[i].charge = electron.charge();
    e[i].isEB = (int)electron.isEB();
    e[i].isEE = (int)electron.isEE();
    e[i].p = electron.p();
    e[i].pt = electron.pt();
    e[i].eta = electron.eta(); 
    e[i].phi = electron.phi(); 
    //e[i].dr04TkSumPt = electron.dr04TkSumPt(); 
    e[i].dr04TkSumPt = electron.dr03TkSumPt();  //!!!!!!!!!!!!!!!!!!!!!! change this later !!
    e[i].dr04EcalRecHitSumEt = electron.dr03EcalRecHitSumEt();
    e[i].dr04HcalTowerSumEt = electron.dr03HcalTowerSumEt();
    e[i].hadronicOverEm = electron.hadronicOverEm(); 
    e[i].scSigmaEtaEta = electron.scSigmaEtaEta();
    e[i].deltaEtaSuperClusterTrackAtVtx = electron.deltaEtaSuperClusterTrackAtVtx(); 
    e[i].deltaPhiSuperClusterTrackAtVtx = electron.deltaPhiSuperClusterTrackAtVtx(); 
    e[i].eSuperClusterOverP = electron.eSuperClusterOverP();

  }

  //genmet
  if(isMC_){
    edm::Handle<edm::View<reco::GenMET> > genmetHandle;
    iEvent.getByLabel(genmetTag_,genmetHandle);
    edm::View<reco::GenMET> genmets = *genmetHandle;
    reco::GenMET gencalomet = genmets[0];
    genmet.et	    = gencalomet.pt();
    genmet.phi	    = gencalomet.phi();
    genmet.sig	    = -1.;
    genmet.sumEt	    = gencalomet.sumEt();
    genmet.emf	    = -1.;
  }

  //calomet
  edm::Handle<edm::View<reco::CaloMET> > calometHandle;
  iEvent.getByLabel(calometTag_,calometHandle);
  edm::View<reco::CaloMET> calomets = *calometHandle;
  reco::CaloMET calometiter = calomets[0];
  calomet.et	    = calometiter.pt();
  calomet.phi	    = calometiter.phi();
  calomet.sig	    = calometiter.significance();
  calomet.sumEt   = calometiter.sumEt();
  calomet.emf     = calometiter.emEtFraction();
  calomet.prob = TMath::Prob(calomet.sig,2);
  //pfmet
  edm::Handle<edm::View<reco::MET> > pfmetHandle;
  iEvent.getByLabel(pfmetTag_,pfmetHandle);
  edm::View<reco::MET> pfmets = *pfmetHandle;
  reco::MET pfmetiter = pfmets[0];
  pfmet.et	    = pfmetiter.pt();
  pfmet.phi	    = pfmetiter.phi();
  pfmet.sig	    = pfmetiter.significance();
  pfmet.sumEt     = pfmetiter.sumEt();
  pfmet.emf       = -1.;
  pfmet.prob = TMath::Prob(pfmet.sig,2);

  //tcmet
  edm::Handle<edm::View<reco::MET> > tcmetHandle;
  iEvent.getByLabel(tcmetTag_,tcmetHandle);
  edm::View<reco::MET> tcmets = *tcmetHandle;
  reco::MET tcmetiter = tcmets[0];
  tcmet.et	    = tcmetiter.pt();
  tcmet.phi	    = tcmetiter.phi();
  tcmet.sig	    = tcmetiter.significance();
  tcmet.sumEt     = tcmetiter.sumEt();
  tcmet.emf       = -1.;
  tcmet.prob = -1;

  //calomet old without thresholds
  edm::Handle<edm::View<reco::CaloMET> > calometoldHandle;
  iEvent.getByLabel(calometoldTag_,calometoldHandle);
  edm::View<reco::CaloMET> calometolds = *calometoldHandle;
  reco::CaloMET calometolditer = calometolds[0];
  calometold.et	    = calometolditer.pt();
  calometold.phi	    = calometolditer.phi();
  calometold.sig	    = calometolditer.significance();
  calometold.sumEt   = calometolditer.sumEt();
  calometold.emf     = calometolditer.emEtFraction();
  calometold.prob = TMath::Prob(calometold.sig,2);
  
  //pfmet wthouth thresholds
  edm::Handle<edm::View<reco::MET> > pfmetoldHandle;
  iEvent.getByLabel(pfmetoldTag_,pfmetoldHandle);
  edm::View<reco::MET> pfmetolds = *pfmetoldHandle;
  reco::MET pfmetolditer = pfmetolds[0];
  pfmetold.et	       = pfmetolditer.pt();
  pfmetold.phi       = pfmetolditer.phi();
  pfmetold.sig       = pfmetolditer.significance();
  pfmetold.sumEt     = pfmetolditer.sumEt();
  pfmetold.emf       = -1.;
  pfmetold.prob = TMath::Prob(pfmetold.sig,2);

  //t1 mets
 //  edm::Handle<edm::View<reco::CaloMET> > t1met05Handle;
//   iEvent.getByLabel(t1met05Tag_,t1met05Handle);
//   edm::View<reco::CaloMET> t1met05s = *t1met05Handle;
//   reco::CaloMET t1met05iter = t1met05s[0];
//   t1met05.et	    = t1met05iter.pt();
//   t1met05.phi	    = t1met05iter.phi();
//   t1met05.sig	    = t1met05iter.significance();
//   t1met05.sumEt   = t1met05iter.sumEt();
//   t1met05.emf     = t1met05iter.emEtFraction();

//   edm::Handle<edm::View<reco::CaloMET> > t1met10Handle;
//   iEvent.getByLabel(t1met10Tag_,t1met10Handle);
//   edm::View<reco::CaloMET> t1met10s = *t1met10Handle;
//   reco::CaloMET t1met10iter = t1met10s[0];
//   t1met10.et	    = t1met10iter.pt();
//   t1met10.phi	    = t1met10iter.phi();
//   t1met10.sig	    = t1met10iter.significance();
//   t1met10.sumEt   = t1met10iter.sumEt();
//   t1met10.emf     = t1met10iter.emEtFraction();

//   edm::Handle<edm::View<reco::CaloMET> > t1met20Handle;
//   iEvent.getByLabel(t1met20Tag_,t1met20Handle);
//   edm::View<reco::CaloMET> t1met20s = *t1met20Handle;
//   reco::CaloMET t1met20iter = t1met20s[0];
//   t1met20.et	    = t1met20iter.pt();
//   t1met20.phi	    = t1met20iter.phi();
//   t1met20.sig	    = t1met20iter.significance();
//   t1met20.sumEt   = t1met20iter.sumEt();
//   t1met20.emf     = t1met20iter.emEtFraction();

//   //t2 mets
//   edm::Handle<edm::View<reco::CaloMET> > t2met05Handle;
//   iEvent.getByLabel(t2met05Tag_,t2met05Handle);
//   edm::View<reco::CaloMET> t2met05s = *t2met05Handle;
//   reco::CaloMET t2met05iter = t2met05s[0];
//   t2met05.et	    = t2met05iter.pt();
//   t2met05.phi	    = t2met05iter.phi();
//   t2met05.sig	    = t2met05iter.significance();
//   t2met05.sumEt   = t2met05iter.sumEt();
//   t2met05.emf     = t2met05iter.emEtFraction();

//   edm::Handle<edm::View<reco::CaloMET> > t2met10Handle;
//   iEvent.getByLabel(t2met10Tag_,t2met10Handle);
//   edm::View<reco::CaloMET> t2met10s = *t2met10Handle;
//   reco::CaloMET t2met10iter = t2met10s[0];
//   t2met10.et	    = t2met10iter.pt();
//   t2met10.phi	    = t2met10iter.phi();
//   t2met10.sig	    = t2met10iter.significance();
//   t2met10.sumEt   = t2met10iter.sumEt();
//   t2met10.emf     = t2met10iter.emEtFraction();

//   edm::Handle<edm::View<reco::CaloMET> > t2met20Handle;
//   iEvent.getByLabel(t2met20Tag_,t2met20Handle);
//   edm::View<reco::CaloMET> t2met20s = *t2met20Handle;
//   reco::CaloMET t2met20iter = t2met20s[0];
//   t2met20.et	    = t2met20iter.pt();
//   t2met20.phi	    = t2met20iter.phi();
//   t2met20.sig	    = t2met20iter.significance();
//   t2met20.sumEt   = t2met20iter.sumEt();
//   t2met20.emf     = t2met20iter.emEtFraction();


  //Calo MHT's

  //Uncorrected MHT's This is only for jb sig use calotowers for jets with emf<0.9 
  //    edm::Handle<edm::View<pat::MHT> > mht05Handle;
  //     iEvent.getByLabel(mht05Tag_,mht05Handle);
  //     edm::View<pat::MHT> mht05s = *mht05Handle;
  //     mht05.et = -1.;
  //     pat::MHT mht05iter	    = mht05s[0];
  //     mht05.njets		    = mht05iter.getNumberOfJets();
  //     mht05.et		    = mht05iter.pt();
  //     mht05.phi		    = mht05iter.phi();
  //     mht05.ht		    = mht05iter.ht();
  //     mht05.sig		    = mht05iter.significance();
  //     mht05.rawmet	    = mht05iter.getUncorMET();
  //     mht05.rawmetsig	    = mht05iter.getUncorMETsignificance();
  //     mht05.jecmet	    = mht05iter.getMET();
  //     mht05.jecmetsig	    = mht05iter.getMETsignificance();

  //     edm::Handle<edm::View<pat::MHT> > mht10Handle;
  //     iEvent.getByLabel(mht10Tag_,mht10Handle);
  //     edm::View<pat::MHT> mht10s = *mht10Handle;
  //     mht10.et = -1.;
  //     pat::MHT mht10iter  = mht10s[0];
  //     mht10.njets	    = mht10iter.getNumberOfJets();
  //     mht10.et	    = mht10iter.pt();
  //     mht10.phi	    = mht10iter.phi();
  //     mht10.ht	    = mht10iter.ht();
  //     mht10.sig	    = mht10iter.significance();
  //     mht10.rawmet	    = mht10iter.getUncorMET();
  //     mht10.rawmetsig	    = mht10iter.getUncorMETsignificance();
  //     mht10.jecmet	    = mht10iter.getMET();
  //     mht10.jecmetsig	    = mht10iter.getMETsignificance();


  //     edm::Handle<edm::View<pat::MHT> > mht20Handle;
  //     iEvent.getByLabel(mht20Tag_,mht20Handle);
  //     edm::View<pat::MHT> mht20s = *mht20Handle;
  //     mht20.et = -1.;
  //     pat::MHT mht20iter  = mht20s[0];
  //     mht20.njets	    = mht20iter.getNumberOfJets();
  //     mht20.et	    = mht20iter.pt();
  //     mht20.phi	    = mht20iter.phi();
  //     mht20.ht	    = mht20iter.ht();
  //     mht20.sig	    = mht20iter.significance();
  //     mht20.rawmet	    = mht20iter.getUncorMET();
  //     mht20.rawmetsig	    = mht20iter.getUncorMETsignificance();
  //     mht20.jecmet	    = mht20iter.getMET();
  //     mht20.jecmetsig	    = mht20iter.getMETsignificance();

  //Uncorrected MHT's This is only for jb sig. use jet resolutions for all jets
  //     edm::Handle<edm::View<pat::MHT> > mht05emfHandle;
  //     iEvent.getByLabel(mht05emfTag_,mht05emfHandle);
  //     edm::View<pat::MHT> mht05emfs = *mht05emfHandle;
  //     mht05emf.et = -1.;
  //     pat::MHT mht05emfiter  = mht05emfs[0];
  //     mht05emf.njets	    = mht05emfiter.getNumberOfJets();
  //     mht05emf.et	    = mht05emfiter.pt();
  //     mht05emf.ht	    = mht05emfiter.ht();
  //     mht05emf.phi	    = mht05emfiter.phi();
  //     mht05emf.sig	    = mht05emfiter.significance();
  //     mht05emf.rawmet	    = mht05emfiter.getUncorMET();
  //     mht05emf.rawmetsig	    = mht05emfiter.getUncorMETsignificance();
  //     mht05emf.jecmet	    = mht05emfiter.getMET();
  //     mht05emf.jecmetsig	    = mht05emfiter.getMETsignificance();

  //    edm::Handle<edm::View<pat::MHT> > mht10emfHandle;
  //     iEvent.getByLabel(mht10emfTag_,mht10emfHandle);
  //     edm::View<pat::MHT> mht10emfs = *mht10emfHandle;
  //     mht10emf.et = -1.;
  //     pat::MHT mht10emfiter  = mht10emfs[0];
  //     mht10emf.njets	    = mht10emfiter.getNumberOfJets();
  //     mht10emf.et	    = mht10emfiter.pt();
  //     mht10emf.ht	    = mht10emfiter.ht();
  //     mht10emf.phi	    = mht10emfiter.phi();
  //     mht10emf.sig	    = mht10emfiter.significance();
  //     mht10emf.rawmet	    = mht10emfiter.getUncorMET();
  //     mht10emf.rawmetsig	    = mht10emfiter.getUncorMETsignificance();
  //     mht10emf.jecmet	    = mht10emfiter.getMET();
  //     mht10emf.jecmetsig	    = mht10emfiter.getMETsignificance();


  //     edm::Handle<edm::View<pat::MHT> > mht20emfHandle;
  //     iEvent.getByLabel(mht20emfTag_,mht20emfHandle);
  //     edm::View<pat::MHT> mht20emfs = *mht20emfHandle;
  //     mht20emf.et = -1.;
  //     pat::MHT mht20emfiter  = mht20emfs[0];
  //     mht20emf.njets	    = mht20emfiter.getNumberOfJets();
  //     mht20emf.et	    = mht20emfiter.pt();
  //     mht20emf.ht	    = mht20emfiter.ht();
  //     mht20emf.phi	    = mht20emfiter.phi();
  //     mht20emf.sig	    = mht20emfiter.significance();
  //     mht20emf.rawmet	    = mht20emfiter.getUncorMET();
  //     mht20emf.rawmetsig	    = mht20emfiter.getUncorMETsignificance();
  //     mht20emf.jecmet	    = mht20emfiter.getMET();
  //     mht20emf.jecmetsig	    = mht20emfiter.getMETsignificance();



  //     //L2L3 05 MHT's This is only for T1 MET ::: Look into config file for more acurate comments
  //     edm::Handle<edm::View<pat::MHT> > l2l3mht05Handle;
  //     iEvent.getByLabel(l2l3mht05Tag_,l2l3mht05Handle);
  //     edm::View<pat::MHT> l2l3mht05s = *l2l3mht05Handle;
  //     l2l3mht05.et = -1.;
  //     pat::MHT l2l3mht05iter  = l2l3mht05s[0];
  //     l2l3mht05.njets	    = l2l3mht05iter.getNumberOfJets();
  //     l2l3mht05.et	    = l2l3mht05iter.pt();
  //     l2l3mht05.ht	    = l2l3mht05iter.ht();
  //     l2l3mht05.phi	    = l2l3mht05iter.phi();
  //     l2l3mht05.sig	    = l2l3mht05iter.significance();
  //     l2l3mht05.rawmet	    = l2l3mht05iter.getUncorMET();
  //     l2l3mht05.rawmetsig	    = l2l3mht05iter.getUncorMETsignificance();
  //     l2l3mht05.jecmet	    = l2l3mht05iter.getMET();
  //     l2l3mht05.jecmetsig	    = l2l3mht05iter.getMETsignificance();

  //L2L3 10 MHT's
//   edm::Handle<edm::View<pat::MHT> > l2l3mht10Handle;
//   iEvent.getByLabel(l2l3mht10Tag_,l2l3mht10Handle);
//   edm::View<pat::MHT> l2l3mht10s = *l2l3mht10Handle;
//   l2l3mht10.et = -1.;
//   pat::MHT l2l3mht10iter  = l2l3mht10s[0];
//   l2l3mht10.njets	    = l2l3mht10iter.getNumberOfJets();
//   l2l3mht10.et	    = l2l3mht10iter.pt();
//   l2l3mht10.ht	    = l2l3mht10iter.ht();
//   l2l3mht10.phi	    = l2l3mht10iter.phi();
//   l2l3mht10.sig	    = l2l3mht10iter.significance();
//   l2l3mht10.rawmet	    = l2l3mht10iter.getUncorMET();
//   l2l3mht10.rawmetsig	    = l2l3mht10iter.getUncorMETsignificance();
//   l2l3mht10.jecmet	    = l2l3mht10iter.getMET();
//   l2l3mht10.jecmetsig	    = l2l3mht10iter.getMETsignificance();

//   //L2L3 20 MHT's
//   edm::Handle<edm::View<pat::MHT> > l2l3mht20Handle;
//   iEvent.getByLabel(l2l3mht20Tag_,l2l3mht20Handle);
//   edm::View<pat::MHT> l2l3mht20s = *l2l3mht20Handle;
//   l2l3mht20.et = -1.;
//   pat::MHT l2l3mht20iter  = l2l3mht20s[0];
//   l2l3mht20.njets	    = l2l3mht20iter.getNumberOfJets();
//   l2l3mht20.et	    = l2l3mht20iter.pt();
//   l2l3mht20.ht	    = l2l3mht20iter.ht();
//   l2l3mht20.phi	    = l2l3mht20iter.phi();
//   l2l3mht20.sig	    = l2l3mht20iter.significance();
//   l2l3mht20.rawmet	    = l2l3mht20iter.getUncorMET();
//   l2l3mht20.rawmetsig	    = l2l3mht20iter.getUncorMETsignificance();
//   l2l3mht20.jecmet	    = l2l3mht20iter.getMET();
//   l2l3mht20.jecmetsig	    = l2l3mht20iter.getMETsignificance();


  //New resolutions
//   edm::Handle<edm::View<pat::MHT> > l2l3mht05newHandle;
//   iEvent.getByLabel(l2l3mht05newTag_,l2l3mht05newHandle);
//   edm::View<pat::MHT> l2l3mht05news = *l2l3mht05newHandle;
//   l2l3mht05new.et = -1.;
//   pat::MHT l2l3mht05newiter  = l2l3mht05news[0];
//   l2l3mht05new.njets	    = l2l3mht05newiter.getNumberOfJets();
//   l2l3mht05new.et	    = l2l3mht05newiter.pt();
//   l2l3mht05new.ht	    = l2l3mht05newiter.ht();
//   l2l3mht05new.phi	    = l2l3mht05newiter.phi();
//   l2l3mht05new.sig	    = l2l3mht05newiter.significance();
//   l2l3mht05new.rawmet	    = l2l3mht05newiter.getUncorMET();
//   l2l3mht05new.rawmetsig  = l2l3mht05newiter.getUncorMETsignificance();
//   l2l3mht05new.jecmet	    = l2l3mht05newiter.getMET();
//   l2l3mht05new.jecmetsig  = l2l3mht05newiter.getMETsignificance();

  //L2L3 10 MHT's
//   edm::Handle<edm::View<pat::MHT> > l2l3mht10newHandle;
//   iEvent.getByLabel(l2l3mht10newTag_,l2l3mht10newHandle);
//   edm::View<pat::MHT> l2l3mht10news = *l2l3mht10newHandle;
//   l2l3mht10new.et = -1.;
//   pat::MHT l2l3mht10newiter  = l2l3mht10news[0];
//   l2l3mht10new.njets	    = l2l3mht10newiter.getNumberOfJets();
//   l2l3mht10new.et	    = l2l3mht10newiter.pt();
//   l2l3mht10new.ht	    = l2l3mht10newiter.ht();
//   l2l3mht10new.phi	    = l2l3mht10newiter.phi();
//   l2l3mht10new.sig	    = l2l3mht10newiter.significance();
//   l2l3mht10new.rawmet	    = l2l3mht10newiter.getUncorMET();
//   l2l3mht10new.rawmetsig	    = l2l3mht10newiter.getUncorMETsignificance();
//   l2l3mht10new.jecmet	    = l2l3mht10newiter.getMET();
//   l2l3mht10new.jecmetsig	    = l2l3mht10newiter.getMETsignificance();

  //L2L3 20 MHT's
//   edm::Handle<edm::View<pat::MHT> > l2l3mht20newHandle;
//   iEvent.getByLabel(l2l3mht20newTag_,l2l3mht20newHandle);
//   edm::View<pat::MHT> l2l3mht20news = *l2l3mht20newHandle;
//   l2l3mht20new.et = -1.;
//   pat::MHT l2l3mht20newiter  = l2l3mht20news[0];
//   l2l3mht20new.njets	    = l2l3mht20newiter.getNumberOfJets();
//   l2l3mht20new.et	    = l2l3mht20newiter.pt();
//   l2l3mht20new.ht	    = l2l3mht20newiter.ht();
//   l2l3mht20new.phi	    = l2l3mht20newiter.phi();
//   l2l3mht20new.sig	    = l2l3mht20newiter.significance();
//   l2l3mht20new.rawmet	    = l2l3mht20newiter.getUncorMET();
//   l2l3mht20new.rawmetsig	    = l2l3mht20newiter.getUncorMETsignificance();
//   l2l3mht20new.jecmet	    = l2l3mht20newiter.getMET();
//   l2l3mht20new.jecmetsig	    = l2l3mht20newiter.getMETsignificance();


  //
  //L2L3 05 MHT's, JB Signficance and MHT
//   edm::Handle<edm::View<pat::MHT> > l2l3mht05emfHandle;
//   iEvent.getByLabel(l2l3mht05emfTag_,l2l3mht05emfHandle);
//   edm::View<pat::MHT> l2l3mht05emfs = *l2l3mht05emfHandle;
//   l2l3mht05emf.et = -1.;
//   pat::MHT l2l3mht05emfiter  = l2l3mht05emfs[0];
//   l2l3mht05emf.njets	    = l2l3mht05emfiter.getNumberOfJets();
//   l2l3mht05emf.et	    = l2l3mht05emfiter.pt();
//   l2l3mht05emf.ht	    = l2l3mht05emfiter.ht();
//   l2l3mht05emf.phi	    = l2l3mht05emfiter.phi();
//   l2l3mht05emf.sig	    = l2l3mht05emfiter.significance();
//   l2l3mht05emf.rawmet	    = l2l3mht05emfiter.getUncorMET();
//   l2l3mht05emf.rawmetsig  = l2l3mht05emfiter.getUncorMETsignificance();
//   l2l3mht05emf.jecmet	    = l2l3mht05emfiter.getMET();
//   l2l3mht05emf.jecmetsig  = l2l3mht05emfiter.getMETsignificance();

  //L2L3 10 MHT's
//   edm::Handle<edm::View<pat::MHT> > l2l3mht10emfHandle;
//   iEvent.getByLabel(l2l3mht10emfTag_,l2l3mht10emfHandle);
//   edm::View<pat::MHT> l2l3mht10emfs = *l2l3mht10emfHandle;
//   l2l3mht10emf.et = -1.;
//   pat::MHT l2l3mht10emfiter  = l2l3mht10emfs[0];
//   l2l3mht10emf.njets	    = l2l3mht10emfiter.getNumberOfJets();
//   l2l3mht10emf.et	    = l2l3mht10emfiter.pt();
//   l2l3mht10emf.ht	    = l2l3mht10emfiter.ht();
//   l2l3mht10emf.phi	    = l2l3mht10emfiter.phi();
//   l2l3mht10emf.sig	    = l2l3mht10emfiter.significance();
//   l2l3mht10emf.rawmet	    = l2l3mht10emfiter.getUncorMET();
//   l2l3mht10emf.rawmetsig	    = l2l3mht10emfiter.getUncorMETsignificance();
//   l2l3mht10emf.jecmet	    = l2l3mht10emfiter.getMET();
//   l2l3mht10emf.jecmetsig	    = l2l3mht10emfiter.getMETsignificance();

  //L2L3 20 MHT's
//   edm::Handle<edm::View<pat::MHT> > l2l3mht20emfHandle;
//   iEvent.getByLabel(l2l3mht20emfTag_,l2l3mht20emfHandle);
//   edm::View<pat::MHT> l2l3mht20emfs = *l2l3mht20emfHandle;
//   l2l3mht20emf.et = -1.;
//   pat::MHT l2l3mht20emfiter  = l2l3mht20emfs[0];
//   l2l3mht20emf.njets	    = l2l3mht20emfiter.getNumberOfJets();
//   l2l3mht20emf.et	    = l2l3mht20emfiter.pt();
//   l2l3mht20emf.ht	    = l2l3mht20emfiter.ht();
//   l2l3mht20emf.phi	    = l2l3mht20emfiter.phi();
//   l2l3mht20emf.sig	    = l2l3mht20emfiter.significance();
//   l2l3mht20emf.rawmet	    = l2l3mht20emfiter.getUncorMET();
//   l2l3mht20emf.rawmetsig	    = l2l3mht20emfiter.getUncorMETsignificance();
//   l2l3mht20emf.jecmet	    = l2l3mht20emfiter.getMET();
//   l2l3mht20emf.jecmetsig	    = l2l3mht20emfiter.getMETsignificance();



  //JPT MHT's

  //    edm::Handle<edm::View<pat::MHT> > jptmht05emfHandle;
  //     iEvent.getByLabel(jptmht05emfTag_,jptmht05emfHandle);
  //     edm::View<pat::MHT> jptmht05emfs = *jptmht05emfHandle;
  //     jptmht05emf.et = -1.;
  //     pat::MHT jptmht05emfiter  = jptmht05emfs[0];
  //     jptmht05emf.njets	    = jptmht05emfiter.getNumberOfJets();
  //     jptmht05emf.et	    = jptmht05emfiter.pt();
  //     jptmht05emf.ht	    = jptmht05emfiter.ht();
  //     jptmht05emf.phi	    = jptmht05emfiter.phi();
  //     jptmht05emf.sig	    = jptmht05emfiter.significance();
  //     jptmht05emf.rawmet	    = jptmht05emfiter.getUncorMET();
  //     jptmht05emf.rawmetsig  = jptmht05emfiter.getUncorMETsignificance();
  //     jptmht05emf.jecmet	    = jptmht05emfiter.getMET();
  //     jptmht05emf.jecmetsig  = jptmht05emfiter.getMETsignificance();

  //     //L2L3 10 MHT's
  //     edm::Handle<edm::View<pat::MHT> > jptmht10emfHandle;
  //     iEvent.getByLabel(jptmht10emfTag_,jptmht10emfHandle);
  //     edm::View<pat::MHT> jptmht10emfs = *jptmht10emfHandle;
  //     jptmht10emf.et = -1.;
  //     pat::MHT jptmht10emfiter  = jptmht10emfs[0];
  //     jptmht10emf.njets	    = jptmht10emfiter.getNumberOfJets();
  //     jptmht10emf.et	    = jptmht10emfiter.pt();
  //     jptmht10emf.ht	    = jptmht10emfiter.ht();
  //     jptmht10emf.phi	    = jptmht10emfiter.phi();
  //     jptmht10emf.sig	    = jptmht10emfiter.significance();
  //     jptmht10emf.rawmet	    = jptmht10emfiter.getUncorMET();
  //     jptmht10emf.rawmetsig	    = jptmht10emfiter.getUncorMETsignificance();
  //     jptmht10emf.jecmet	    = jptmht10emfiter.getMET();
  //     jptmht10emf.jecmetsig	    = jptmht10emfiter.getMETsignificance();

  //     //L2L3 20 MHT's
  //     edm::Handle<edm::View<pat::MHT> > jptmht20emfHandle;
  //     iEvent.getByLabel(jptmht20emfTag_,jptmht20emfHandle);
  //     edm::View<pat::MHT> jptmht20emfs = *jptmht20emfHandle;
  //     jptmht20emf.et = -1.;
  //     pat::MHT jptmht20emfiter  = jptmht20emfs[0];
  //     jptmht20emf.njets	    = jptmht20emfiter.getNumberOfJets();
  //     jptmht20emf.et	    = jptmht20emfiter.pt();
  //     jptmht20emf.ht	    = jptmht20emfiter.ht();
  //     jptmht20emf.phi	    = jptmht20emfiter.phi();
  //     jptmht20emf.sig	    = jptmht20emfiter.significance();
  //     jptmht20emf.rawmet	    = jptmht20emfiter.getUncorMET();
  //     jptmht20emf.rawmetsig   = jptmht20emfiter.getUncorMETsignificance();
  //     jptmht20emf.jecmet	    = jptmht20emfiter.getMET();
  //     jptmht20emf.jecmetsig   = jptmht20emfiter.getMETsignificance();

  edm::Handle<std::vector<pat::Jet> > patJetHandle;
  iEvent.getByLabel(patJetTag_,patJetHandle);
    
  size_t patJetCounter=0;
  for(std::vector<pat::Jet>::const_iterator jet_iter = patJetHandle->begin(); jet_iter!=patJetHandle->end() && patJetCounter<2; ++jet_iter){
    patjet[patJetCounter].e = jet_iter->energy();
    patjet[patJetCounter].et = jet_iter->pt();
    patjet[patJetCounter].eta = jet_iter->eta();
    patjet[patJetCounter].phi = jet_iter->phi();
    patjet[patJetCounter].emf=jet_iter->emEnergyFraction();
    patjet[patJetCounter].fHPD = jet_iter->jetID().fHPD;
    patjet[patJetCounter].n90hits = jet_iter->jetID().n90Hits;

    patJetCounter++;
  }
  //calo jets
  edm::Handle<edm::View<reco::CaloJet > > jetHandle;
  iEvent.getByLabel(jetTag_, jetHandle);
  edm::View<reco::CaloJet> jets = *jetHandle;
  idx[0]=999999; idx[1]=999999;
  pmax1=0; pmax2=0;
  for (unsigned int i=0; i<jets.size(); ++i){
    float pt=jets[i].pt();
    if(pt>pmax1 && pt>pmax2){
      idx[0]=i;
      pmax1= pt;
    }
    else if (pt<pmax1 && pt>pmax2){
      idx[1] = i;
      pmax2 = pt;
    }
  }
//   const JetCorrector* l2l3 = JetCorrector::getJetCorrector(jetCorrectorTag_,iSetup); 
//   edm::Handle<reco::JetIDValueMap> hJetIDMap;
//   iEvent.getByLabel(jetIDTag_, hJetIDMap);
//   for (int i=0; i<2; i++){
//     if(idx[i]==999999) continue;
//     reco::CaloJet jet = jets[idx[i]];
//     edm::RefToBase<reco::CaloJet> jetRef = jetHandle->refAt(idx[i]);
//     reco::JetID jetId = (*hJetIDMap)[jetRef];
//     calojet[i].e = jet.energy();
//     calojet[i].et = jet.pt();
//     calojet[i].eta = jet.eta();
//     calojet[i].phi = jet.phi();
//     calojet[i].emf = jet.emEnergyFraction();
//     calojet[i].l2l3 = l2l3->correction(jet.p4());
//     calojet[i].fHPD = jetId.fHPD;
//     calojet[i].n90hits = jetId.n90Hits;
//     calojet[i].dpt = 0.;
//     calojet[i].dphi = 0.;
//   } 

  //gen jets
  if(isMC_){
    edm::Handle<edm::View<reco::GenJet > > genjetHandle;
    iEvent.getByLabel(genjetTag_, genjetHandle);
    edm::View<reco::GenJet> genjets = *genjetHandle;
    idx[0]=999999; idx[1]=999999;
    pmax1=0; pmax2=0;
    for (unsigned int i=0; i<genjets.size(); ++i){
      float pt=genjets[i].pt();
      if(pt>pmax1 && pt>pmax2){
	idx[0]=i;
	pmax1= pt;
      }
      else if (pt<pmax1 && pt>pmax2){
	idx[1] = i;
	pmax2 = pt;
      }
    }
    for (int i=0; i<2; i++){
      if(idx[i]==999999) continue;
      reco::GenJet jet = genjets[idx[i]];
      genjet[i].e = jet.energy();
      genjet[i].et = jet.pt();
      genjet[i].eta = jet.eta();
      genjet[i].phi = jet.phi();
      genjet[i].emf = -1;
      genjet[i].l2l3 = -1;
      genjet[i].fHPD = -1;
      genjet[i].n90hits = -1;
      genjet[i].dpt = -1;
      genjet[i].dphi = -1;
    } 
  }

  results_tree -> Fill();
  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
CombinedNtupler::beginJob()
{

  edm::Service<TFileService> fs;
  eventcount = fs->make<TH1D>("eventcount","events processed",1,-0.5,0.5);
  results_tree = new TTree("events", "events");
  results_tree -> Branch("run", &run, "run/I");
  results_tree -> Branch("lumi", &lumi, "lumi/I");
  results_tree -> Branch("event", &event, "event/I");
  results_tree -> Branch("generator", &geninfo, "pid/I:pthat/F:alphaQCD/F:alphaQED/F:scalePDF/F:id1/I:id2/I:x1/F:x2/F:xPDF1/F:xPDF2/F");
  

  results_tree -> Branch("hlt", &hlt, "hlt/O");
  for(std::map<std::string,int>::iterator iter=triggervals_.begin(); iter!=triggervals_.end();++iter){
    TString name2="trigger_";
    name2+=iter->first;
    name2.ReplaceAll(" ","");
    TString secondname2=name2;
    secondname2+="/I";
    std::cout << "booking branch: " << secondname2 << std::endl;
    results_tree->Branch(name2.Data(),&(iter->second),secondname2.Data());
  }

  //electrons
  for(int i=0; i<nElectrons; i++){
    stringstream ss;
    ss<<i;
    TString name=("e" + ss.str()).c_str(); 
    results_tree -> Branch(name, &e[i], "charge/I:isEB/I:isEE/I:p/F:pt/F:eta/F:phi/F:dr04TkSumPt/F:dr04EcalRecHitSumEt/F:dr04HcalTowerSumEt/F:hadronicOverEm/F:scSigmaEtaEta/F:deltaEtaSuperClusterTrackAtVtx/F:deltaPhiSuperClusterTrackAtVtx/F:eSuperClusterOverP/F");
  }

  //calojets
  for(int i=0; i<nJets; i++){
    stringstream ss;
    ss<<i;
    TString name=("calojet" + ss.str()).c_str(); 
    results_tree -> Branch(name, &calojet[i], "e/F:et:eta:phi:dpt:dphi:emf:l2l3:fHPD:n90hits/I");
  }
  
  //patjets
  for(int i=0; i<nJets; i++){
    stringstream ss;
    ss<<i;
    TString name=("patjet" + ss.str()).c_str(); 
    results_tree -> Branch(name, &patjet[i], "e/F:et:eta:phi:dpt:dphi:emf:l2l3:fHPD:n90hits/I");
  }

  //genjets
  for(int i=0; i<nJets; i++){
    stringstream ss;
    ss<<i;
    TString name=("genjet" + ss.str()).c_str(); 
    results_tree -> Branch(name, &genjet[i], "e/F:et:eta:phi:dpt:dphi:emf:l2l3:fHPD:n90hits/I");
  }

  //met's, mht's
  results_tree -> Branch("genmet", &genmet, "et/F:phi:sig:sumEt:emf");
  results_tree -> Branch("calomet", &calomet, "et/F:phi:sig:sumEt:emf");
  results_tree -> Branch("pfmet", &pfmet, "et/F:phi:sig:sumEt:emf");
  results_tree -> Branch("tcmet", &tcmet, "et/F:phi:sig:sumEt:emf");
  results_tree -> Branch("calometold", &calometold, "et/F:phi:sig:sumEt:emf");
  results_tree -> Branch("pfmetold", &pfmetold, "et/F:phi:sig:sumEt:emf");

  results_tree -> Branch("t1met05", &t1met05, "et/F:phi:sig:sumEt:emf");
  results_tree -> Branch("t1met10", &t1met10, "et/F:phi:sig:sumEt:emf");
  results_tree -> Branch("t1met20", &t1met20, "et/F:phi:sig:sumEt:emf");
  results_tree -> Branch("t2met05", &t2met05, "et/F:phi:sig:sumEt:emf");
  results_tree -> Branch("t2met10", &t2met10, "et/F:phi:sig:sumEt:emf");
  results_tree -> Branch("t2met20", &t2met20, "et/F:phi:sig:sumEt:emf");

  results_tree -> Branch("mht05", &mht05, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("mht10", &mht10, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("mht20", &mht20, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");

  results_tree -> Branch("mht05emf", &mht05emf, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("mht10emf", &mht10emf, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("mht20emf", &mht20emf, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");

  results_tree -> Branch("l2l3mht05", &l2l3mht05, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("l2l3mht10", &l2l3mht10, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("l2l3mht20", &l2l3mht20, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");

  results_tree -> Branch("l2l3mht05new", &l2l3mht05new, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("l2l3mht10new", &l2l3mht10new, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("l2l3mht20new", &l2l3mht20new, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");

  results_tree -> Branch("l2l3mht05emf", &l2l3mht05emf, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("l2l3mht10emf", &l2l3mht10emf, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("l2l3mht20emf", &l2l3mht20emf, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");

  results_tree -> Branch("jptmht05emf", &jptmht05emf, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("jptmht10emf", &jptmht10emf, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");
  results_tree -> Branch("jptmht20emf", &jptmht20emf, "njets/I:et/F:phi:ht:sig:rawmet:rawmetsig:jecmet:jecmetsig");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
CombinedNtupler::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CombinedNtupler);
