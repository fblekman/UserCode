// -*- C++ -*-
//
// Package:    ZmumuTauRecoil
// Class:      ZmumuTauRecoil
// 
/**\class ZmumuTauRecoil ZmumuTauRecoil.cc myana/ZmumuTauRecoil/src/ZmumuTauRecoil.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Freya Blekman
//         Created:  Fri Apr 18 23:57:53 CEST 2008
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
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/TauReco/interface/BaseTau.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/IsolatedTauTagInfo.h"

#include "TH1D.h"
#include "TLorentzVector.h"
#include <map>
#include "TTree.h"
#include <iostream>
#include <vector>
//
// class decleration
//

class ZmumuTauRecoil : public edm::EDAnalyzer {
   public:
      explicit ZmumuTauRecoil(const edm::ParameterSet&);
      ~ZmumuTauRecoil();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  std::map<std::string,TH1D*> histos_;
  TTree *atree_;
  TLorentzVector *lv_z;
  TLorentzVector *lv_j;
  TLorentzVector *lv_gj;
  TLorentzVector *lv_tau;
  float passedtautags[5];
  
  edm::InputTag tauSrc_;
  float muPtcut;
  float jetPtcut;
  float tauPtcut;
  float tauDeltaRcut;
  float MZ;
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
ZmumuTauRecoil::ZmumuTauRecoil(const edm::ParameterSet& iConfig):muPtcut(20),jetPtcut(30),tauPtcut(30),tauDeltaRcut(0.4),MZ(91.1876),
								 tauSrc_(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc"))

{
   //now do what ever initialization is needed

}


ZmumuTauRecoil::~ZmumuTauRecoil()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ZmumuTauRecoil::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   // first: get muons. Check that there's two and reconstruct the z-mass:

   edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByLabel("selectedLayer1Muons",muonHandle);
   std::vector<pat::Muon> muons = *muonHandle;
   
   edm::Handle<std::vector<pat::Jet> > jetHandle;
   iEvent.getByLabel("selectedLayer1Jets",jetHandle);
   std::vector<pat::Jet> jets = *jetHandle;

   edm::Handle<std::vector<pat::Tau> > tauHandle;
   iEvent.getByLabel("selectedLayer1Taus",tauHandle);
   std::vector<pat::Tau> taus = *tauHandle;

   TLorentzVector zcand, bestzcand,mu1, mu2, jet, bestjet, mcjet, mcpar, lvworker,lvworker2;
   size_t nzcand=0;
   size_t ngoodmus=0;
   size_t cutstep=0;
   histos_["hist_totalevents"]->Fill(cutstep); cutstep++;
   for(size_t i=0; i<muons.size(); ++i){
     
     histos_["hist_muPt_precut"]->Fill(muons[i].pt());
     histos_["hist_muEta_precut"]->Fill(muons[i].eta());
     if(muons[i].pt()<muPtcut)
       continue;
     ngoodmus++;
     histos_["hist_muPt_postcut"]->Fill(muons[i].pt());
     histos_["hist_muEta_postcut"]->Fill(muons[i].eta());
     for(size_t j=i+1; j<muons.size(); ++j){
       if(muons[j].pt()<muPtcut)
	 continue;
       //       std::cout << muons[i].charge() << " " << muons[j].charge() << std::endl;
       if(muons[i].charge()*muons[j].charge()<0){
	 mu1.SetPtEtaPhiE(muons[i].pt(),muons[i].eta(),muons[i].phi(),muons[i].energy());
	 mu2.SetPtEtaPhiE(muons[j].pt(),muons[j].eta(),muons[j].phi(),muons[j].energy());
	 zcand=mu1+mu2;
	
	 if(nzcand=0){
	   bestzcand=zcand;
	   
	 }
	 else
	   if(fabs(zcand.M()-MZ)<fabs(bestzcand.M()-MZ)){
	     bestzcand=zcand;
	   }
	 nzcand++;
       }
     }
   }
   
   histos_["hist_Nmu_precut"]->Fill(muons.size());
   histos_["hist_Nmu_postcut"]->Fill(ngoodmus);
   if(nzcand==0)
     return; 
   histos_["hist_totalevents"]->Fill(cutstep); cutstep++;

   histos_["hist_muZPt"]->Fill(bestzcand.Pt());
   histos_["hist_muZM"]->Fill(bestzcand.M());
   histos_["hist_muZPhi"]->Fill(bestzcand.Phi());
   std::cout << "number of muons = " << muons.size() << std::endl;
   // now loop over the jets:
   size_t ngoodjet=0;
   size_t goodjetindex=0;
   size_t leadingjetindex=0;
   for(size_t i=0;i< jets.size(); i++){
     
     histos_["hist_jetPt_precut"]->Fill(jets[i].pt());
     histos_["hist_jetEta_precut"]->Fill(jets[i].eta());
     if(jets[i].pt()<jetPtcut)
       continue;
     histos_["hist_jetPt_postcut"]->Fill(jets[i].pt());
     histos_["hist_jetEta_postcut"]->Fill(jets[i].eta());
     jet.SetPtEtaPhiE(jets[i].pt(),jets[i].eta(),jets[i].phi(),jets[i].energy());
     if(jets[i].pt()>jets[leadingjetindex].pt()){
       leadingjetindex=i;
     }
     if(ngoodjet=0){
       bestjet=jet;
       goodjetindex=i;
     }
     else{
       if(fabs(jet.DeltaR(bestzcand)-TMath::Pi())<fabs(bestjet.DeltaR(bestzcand)-TMath::Pi())){
	 bestjet=jet;
	 goodjetindex=i;
       }
     }
     
     ngoodjet++;
   }
   
   if(ngoodjet==0)
     return;
  
   histos_["hist_totalevents"]->Fill(cutstep); cutstep++;
   
   bestjet.SetPtEtaPhiE(jets[goodjetindex].pt(),jets[goodjetindex].eta(),jets[goodjetindex].phi(),jets[goodjetindex].energy());
   const reco::GenJet *genjet = jets[goodjetindex].genJet();
   mcjet.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(),genjet->energy());
   mcpar.SetPtEtaPhiE(jets[goodjetindex].genParton()->pt(),jets[goodjetindex].genParton()->eta(),jets[goodjetindex].genParton()->phi(),jets[goodjetindex].genParton()->energy());
   histos_["hist_jetflavor"]->Fill(jets[goodjetindex].partonFlavour());
   histos_["hist_genjetEta"]->Fill(mcjet.Eta());
   histos_["hist_genjetPt"]->Fill(mcjet.Pt());
   histos_["hist_jetpartonEta"]->Fill(mcpar.Eta());
   histos_["hist_jetpartonPt"]->Fill(mcpar.Pt());
   histos_["hist_dPhi_jetZ"]->Fill(bestjet.DeltaR(bestzcand)); 

  
   // and now do tau tagging:
   TLorentzVector tlv(0,0,0,0);
   size_t ngoodtaus=0;
   size_t ntaus=0;
   size_t nmatchedtaus=0;
   size_t ntaggedtaus=0;
  
   if(1){
     Handle<reco::IsolatedTauTagInfoCollection> tauTagInfoHandle;
     iEvent.getByLabel(tauSrc_, tauTagInfoHandle);
     
      
     float taucuts1[6]={0.1,0.07,0.45,6.,1.,2.};//Rm, Rs, Ri, pTltrCut, pTiCut,dZtr
     float taucuts2[6]={0.1,0.14,0.45,6.,1.,2.};//Rm, Rs, Ri, pTltrCut, pTiCut,dZtr
     float taucuts3[6]={0.12,0.14,0.45,6.,1.,2.};//Rm, Rs, Ri, pTltrCut, pTiCut,dZtr
     for(size_t ii=0; ii<5; ii++)
       passedtautags[ii]=0;
     for (reco::IsolatedTauTagInfoCollection::const_iterator itau = tauTagInfoHandle->begin(); itau != tauTagInfoHandle->end(); ++itau) {
       
       ntaus++;
       const reco::CaloJet* cj = dynamic_cast<reco::CaloJet*>( const_cast<reco::Jet*>( itau->jet().get() ) );
       tlv.SetXYZT(0,0,0,0);
       TLorentzVector recoTauJet(itau->jet()->px(), itau->jet()->py(),itau->jet()->pz(),itau->jet()->energy());
       histos_["hist_tauPt_precut"]->Fill(recoTauJet.Pt());
       histos_["hist_tauEta_precut"]->Fill(recoTauJet.Eta());
       if(recoTauJet.Pt()<tauPtcut)
	 continue;     
       histos_["hist_tauPt_postcut"]->Fill(recoTauJet.Pt());
       histos_["hist_tauEta_postcut"]->Fill(recoTauJet.Eta());
       ngoodtaus++;
       if(recoTauJet.DeltaR(bestjet)>tauDeltaRcut)
	 continue;
       nmatchedtaus++;
       histos_["hist_tauPt_postmatch"]->Fill(recoTauJet.Pt());
       histos_["hist_tauEta_postmatch"]->Fill(recoTauJet.Eta());
       tlv=recoTauJet;
       passedtautags[0]=1;
       if(itau->discriminator())
	 passedtautags[1]=1;
       if(itau->discriminator(taucuts1[0],taucuts1[1],taucuts1[2],taucuts1[3],taucuts1[4],taucuts1[5]))
	 passedtautags[2]=1;
       if(itau->discriminator(taucuts2[0],taucuts2[1],taucuts2[2],taucuts2[3],taucuts2[4],taucuts2[5]))
	 passedtautags[3]=1;
       if(itau->discriminator(taucuts3[0],taucuts3[1],taucuts3[2],taucuts3[3],taucuts3[4],taucuts3[5]))
	 passedtautags[4]=1;
       
     }
   }
   histos_["hist_Ntau_precut"]->Fill(ntaus);
   histos_["hist_Ntau_postcut"]->Fill(ngoodtaus);
   histos_["hist_Ntau_postmatch"]->Fill(nmatchedtaus);
   
   if(nmatchedtaus>0)
     histos_["hist_totalevents"]->Fill(cutstep); cutstep++;
   // and fill ntuple
   lv_z->SetVectM(bestzcand.Vect(),bestzcand.M());
  
   lv_j->SetVectM(bestjet.Vect(),bestjet.M());
   lv_gj->SetVectM(mcjet.Vect(),mcjet.M());
   lv_tau->SetVectM(tlv.Vect(),tlv.M());
   atree_->Fill();

   
//    for(size_t i=0; i<taus.size(); i++){
//      histos_["hist_tauPt_precut"]->Fill(taus[i].pt());
//      histos_["hist_tauEta_precut"]->Fill(taus[i].eta());
//      if(taus[i].pt()<tauPtcut)
//        continue;
//      histos_["hist_tauPt_postcut"]->Fill(taus[i].pt());
//      histos_["hist_tauEta_postcut"]->Fill(taus[i].eta());
//      lvworker.SetPtEtaPhiE(taus[i].pt(),taus[i].eta(),taus[i].phi(),taus[i].energy());
//      ngoodtaus++;
//      if(lvworker.DeltaR(bestjet)<tauDeltaRcut){
//        tlv=lvworker;
//        goodtauindex=i;

//               1.  "standard PTDR": Rm=0.1, Rs=0.07, Ri=0.45, pTi=1 GeV,
//                                    Ptltr=6 GeV, DZtr=2mm
//               2.   relaxed: Rs=0.14, the rest is the same
//               3.   more relaxed: Rs= 0.14, Rm=0.12, the rest is the same

//        std::cout << "tau info: " <<std::endl;
//        std::cout << "emEnergyFraction" << taus[i].emEnergyFraction() << std::endl;
//        std::cout << "eOverP" << taus[i].eOverP() << std::endl;
//        std::cout << "leadEoverP" << taus[i].leadEoverP() << std::endl;
//        std::cout << "hHotOverP" << taus[i].hHotOverP() << std::endl;
//        std::cout << "hTotOverP" << taus[i].hTotOverP() << std::endl;
//      }
//    }
}


// ------------ method called once each job just before starting event loop  ------------
void 
ZmumuTauRecoil::beginJob(const edm::EventSetup&)
{  
  edm::Service<TFileService> fs;
  TFileDirectory muDir = fs->mkdir("muonInfo");
  histos_["hist_Nmu_precut"]=muDir.make<TH1D>("hist_Nmu_precut","PAT muons",10,0,10);
  histos_["hist_Nmu_postcut"]=muDir.make<TH1D>("hist_Nmu_postcut","PAT muons",10,0,10);
  histos_["hist_muPt_precut"]=muDir.make<TH1D>("hist_muPt_precut","PAT muons",40,0,400);
  histos_["hist_muPt_postcut"]=muDir.make<TH1D>("hist_muPt_postcut","PAT muons",40,0,400);
  histos_["hist_muEta_precut"]=muDir.make<TH1D>("hist_muEta_precut","PAT muons",20,-2.5,2.5);
  histos_["hist_muEta_postcut"]=muDir.make<TH1D>("hist_muEta_postcut","PAT muons",20,-2.5,2.5);

  TFileDirectory tauDir = fs->mkdir("tauonInfo");
  histos_["hist_Ntau_precut"]=tauDir.make<TH1D>("hist_Ntau_precut","PAT taus",10,0,10);
  histos_["hist_Ntau_postcut"]=tauDir.make<TH1D>("hist_Ntau_postcut","PAT taus",10,0,10);
  histos_["hist_Ntau_postmatch"]=tauDir.make<TH1D>("hist_Ntau_postmatch","PAT taus",10,0,10);
  histos_["hist_tauPt_precut"]=tauDir.make<TH1D>("hist_tauPt_precut","PAT taus",40,0,400);
  histos_["hist_tauPt_postcut"]=tauDir.make<TH1D>("hist_tauPt_postcut","PAT taus",40,0,400);
  histos_["hist_tauPt_postmatch"]=tauDir.make<TH1D>("hist_tauPt_postmatch","PAT taus",40,0,400);
  histos_["hist_tauEta_precut"]=tauDir.make<TH1D>("hist_tauEta_precut","PAT taus",20,-2.5,2.5);
  histos_["hist_tauEta_postcut"]=tauDir.make<TH1D>("hist_tauEta_postcut","PAT taus",20,-2.5,2.5);
  histos_["hist_tauEta_postmatch"]=tauDir.make<TH1D>("hist_tauEta_postmatch","PAT taus",20,-2.5,2.5);

  TFileDirectory jetDir =  fs->mkdir("jetInfo");
  
  histos_["hist_jetflavor"]=jetDir.make<TH1D>("hist_jetflavor","PAT jets",40,-20,20);
  histos_["hist_jetPt_precut"]=jetDir.make<TH1D>("hist_jetPt_precut","PAT jets",40,0,400);
  histos_["hist_jetPt_postcut"]=jetDir.make<TH1D>("hist_jetPt_postcut","PAT jets",40,0,400);
  histos_["hist_jetEta_precut"]=jetDir.make<TH1D>("hist_jetEta_precut","PAT jets",20,-2.5,2.5);
  histos_["hist_jetEta_postcut"]=jetDir.make<TH1D>("hist_jetEta_postcut","PAT jets",20,-2.5,2.5);
  histos_["hist_genjetEta"]=jetDir.make<TH1D>("hist_genjetEta","PAT jets",20,-2.5,2.5);
  histos_["hist_genjetPt"]=jetDir.make<TH1D>("hist_genjetPt","PAT jets",40,0,400);
  histos_["hist_jetpartonEta"]=jetDir.make<TH1D>("hist_jetpartonEta","PAT jets",20,-2.5,2.5);
  histos_["hist_jetpartonPt"]=jetDir.make<TH1D>("hist_jetpartonPt","PAT jets",40,0,400);
  
  histos_["hist_totalevents"]=fs->make<TH1D>("hist_totalevents","numbers of events passing each cut",10,0,10);
  histos_["hist_dPhi_jetZ"]=fs->make<TH1D>("hist_dPhi_jetZ","Z+jet #rightarrow #mu #mu jet",20, 0, TMath::Pi());
  histos_["hist_muZPt"]=fs->make<TH1D>("hist_muZPt", "Z #rightarrow #mu #mu candidates", 20, 0, 200 );
  histos_["hist_muZPhi"]=fs->make<TH1D>("hist_muZPhi", "Z #rightarrow #mu #mu candidates", 20, -TMath::Pi(), TMath::Pi());
  histos_["hist_muZM"]=fs->make<TH1D>("hist_muZM","Z #rightarrow #mu #mu candidates", 20, 0, 200 );
  //  histos_[]=fs->make<TH1D>("","",,,);
   atree_= new TTree("t","t");
   lv_z = new TLorentzVector(0,0,0,0);
   atree_->Branch("z",&lv_z);
   lv_j = new TLorentzVector(0,0,0,0);
   atree_->Branch("jet",&lv_j);
   lv_gj= new TLorentzVector(0,0,0,0);
   atree_->Branch("genjet",&lv_gj);
   lv_tau = new TLorentzVector(0,0,0,0);
   atree_->Branch("tau",&lv_tau);
   atree_->Branch("tautags",passedtautags,"tautags[5]/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZmumuTauRecoil::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZmumuTauRecoil);
