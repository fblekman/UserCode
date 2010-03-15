#include <iostream>
#include <fstream>
#include <TMatrix.h>
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
#include <TMatrix.h>
#include <string>
#include "SigInputObj.h"
#include "significanceAlgo.h"
#include "DataFormats/Math/interface/deltaR.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
// INCLUDE ALL CMSSW FORMATS HERE!
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/MET.h" 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"


#endif // end of INCLUDE ALL CMSSW FORMATS HERE!

void makeplot_fwlite(TString fileliststring, TString outfilestring){
  std::vector<std::pair<uint64_t,uint64_t> > badevents;
  ifstream badeventlist("runlist_ecalselection.txt");
  if(badeventlist.is_open()){
    std::cout << "adding bad event list" << std::endl;
    while(! badeventlist.eof() ){
      uint64_t run,evt;
      badeventlist >> run >> evt ;
      //      std::cout << run << " " << evt << std::endl;
      std::pair<uint64_t,uint64_t> badone(run,evt);
      badevents.push_back(badone);
    }
  }
  badeventlist.close();
 
  std::string emptystring="";
  std::vector<metsig::SigInputObj> EventVec1;
  std::vector<metsig::SigInputObj> EventVec2;
  std::vector<metsig::SigInputObj> EventVec3;

  std::vector<uint64_t> badruns;
 // good runs:
  std::vector<uint64_t> goodruns;
  ifstream goodrunlist("goodrunlist.txt");
  if(goodrunlist.is_open()){
    std::cout << "adding good run list" << std::endl;
    while(!goodrunlist.eof()){
      uint64_t run;
      goodrunlist >> run;
      goodruns.push_back(run);
    }
  }
  goodrunlist.close();
 
  std::vector<std::string> fileNames;
  ifstream filelist(fileliststring);
  if(filelist.is_open()){
    std::cout << "adding files..." << std::endl;
    while(!filelist.eof()){
      std::string tempname;
      filelist >> tempname;
      TString command="stager_get -M ";
      command+=tempname;
      command+=" &";
      command.ReplaceAll("rfio:","");
      //      std::cout << command << std::endl;
      if(tempname.size()>0 && command.Contains(".root")){
	//	gSystem->Exec(command.Data());     
	fileNames.push_back(tempname);
      }
    }
  }


  std::cout << "reading " << fileNames.size() << " files." << std::endl;
  filelist.close();
 
  fwlite::ChainEvent ev(fileNames);  

  fwlite::Handle<L1GlobalTriggerReadoutRecord>  handle_techbits;
  fwlite::Handle<std::vector<reco::CaloMET> > handle_mets;
  fwlite::Handle<std::vector<reco::PFMET> > handle_pfmets;
  fwlite::Handle<std::vector<reco::Track> > handle_tracks;
  fwlite::Handle<std::vector<reco::Vertex> > handle_vtxs;
  fwlite::Handle<reco::PFCandidateCollection> handle_pfcands;
  fwlite::Handle<CaloTowerCollection> handle_ct;

  TFile *outputfile = new TFile(outfilestring,"recreate");
  outputfile->cd();
//   TH1D *h_vtx_z = new TH1D("h_vtx_z","",100,-20,20);
//   h_vtx_z->SetXTitle("PV z (cm)");
//   h_vtx_z->SetYTitle("events");
 
//   TH1D *h_vtx_ndof = new TH1D("h_vtx_ndof","",100,0,100);
//   h_vtx_ndof->SetXTitle("PV ndof");
//   h_vtx_ndof->SetYTitle("events");

//   TH1D *h_fractrackgood = new TH1D("h_fractrackgood","",100,0,1);
//   h_fractrackgood->SetXTitle("#frac{good tracks}{all tracks}");
//   h_fractrackgood->SetYTitle("events/0.01");

//   TH1D *h_met = new TH1D("h_met","",100,0,10);
//   h_met->SetXTitle("missing E_{T} (GeV)");
//   h_met->SetYTitle("events");
//   //  outputfile->Add(h_met);

//   TH1D *h_metsignif = new TH1D("h_metsignif","",104,-2,50);
//   h_metsignif->SetXTitle("missing E_{T} significance");
//   h_metsignif->SetYTitle("events");

//   TH1D *h_metsum = new TH1D("h_metsum","",100,0,100);
//   h_metsum->SetXTitle("#Sigma E_{T}");
//   h_metsum->SetYTitle("events");

//   TH2D *h_metsignfsumEt = new TH2D("h_metsignfsumEt","",100,0,100,100,0,20);
//   h_metsignfsumEt->SetXTitle("#Sigma E_{T}");
//   h_metsignfsumEt->SetYTitle("missing E_{T} significance");
//   TH2D *h_metsignfsqrtsumEt = new TH2D("h_metsignfsqrtsumEt","",100,0,5,100,0,20);
  
//   h_metsignfsqrtsumEt->SetXTitle("missing E_{T}/sqrt(#Sigma E_{T})");
//   h_metsignfsqrtsumEt->SetYTitle("missing E_{T} significance");

  TTree *littletree = new TTree("t","t");
  struct stuff{
    Float_t met;
    Float_t sig;
    Float_t phi;
    Float_t sumEt;
  };
  stuff container;
  stuff pfcontainer;
  stuff pftrkcontainer;
  stuff pfcalcontainer;
  
  struct eventstuff{
    Int_t runnum;
    Int_t evnum;
  };
  eventstuff evcontainer;

  struct hfstuff{
    Float_t L;
    Float_t S;
    Float_t EM;
    Float_t HAD;
    Float_t alpha;
    Float_t HO;
    Float_t e;
    Float_t et;
    Float_t eta;
    Float_t phi;
    Float_t ieta;
    Float_t iphi;
  };

  hfstuff hfcontainer;
  struct pfstuff{
    Float_t p;
    Float_t pt;
    Float_t e;
    Float_t et;
    Float_t phi;
    Float_t eta;
    Float_t ehad;
    Float_t eem;
    Float_t de;
    Float_t dphi;
    Int_t pid;    
  };
  pfstuff pfccontainer;
  TTree *hftree = new TTree("hft","hft");
  //  hftree->Branch("calotow",&hfcontainer,"L/F:S/F:EM/F:HAD/F:alpha/F:HO/F:e/F:et/F:eta/F:phi/F:ieta/F:iphi/F");
  hftree->Branch("pfc",&pfccontainer,"p/F:pt/F:e/F:et/F:phi/F:eta/F:ehad/F:eem/F:de/F:dphi/F:pid/I");
  littletree->Branch("met",&container,"met/F:sig/F:phi/F:sumEt/F");
  littletree->Branch("ev",&evcontainer,"run/I:ev/I");
  littletree->Branch("pfmet",&pfcontainer,"met/F:sig/F:phi/F:sumEt/F");
  littletree->Branch("pfmettrk",&pftrkcontainer,"met/F:sig/F:phi/F:sumEt/F");
  littletree->Branch("pfmetcal",&pfcalcontainer,"met/F:sig/F:phi/F:sumEt/F");
  
  double trackfrac;
  size_t iev=0;
  size_t nevmax=-1;//100000;
  size_t repeat=1;
  size_t ii,jj,kk;
  std::pair<uint64_t,uint64_t> currentevent;
  bool ismc;
  // various iterators.
  std::vector<reco::Vertex>::const_iterator vtxs;
  std::vector<reco::CaloMET>::const_iterator mets;
  std::vector<reco::PFMET>::const_iterator pfmets;
  std::vector<reco::Track>::const_iterator trk;
  reco::PFCandidateCollection::const_iterator pfc;
  //     #PF:
//     # type 1: charged hadron - essentially tracking resolutions
  //PF_EtResType1 = cms.vdouble(0.05,0,0),
//     PF_PhiResType1 = cms.vdouble(0.002), 
  std::map<int,std::vector<double> > pfvalues;
  std::vector<double> bla(4,0);
  bla[2]=0.007;
  bla[1]=0;
  bla[0]=0;
  bla[3]=0.01;
  pfvalues[1]=bla;
  
//     # type 2: EM with track (electron) - essentially tracking resolution
//     PF_EtResType2 = cms.vdouble(0.05,0,0),
//     PF_PhiResType2 = cms.vdouble(0.002),
  bla[2]=0.007;
  bla[1]=0;
  bla[0]=0;
  bla[3]=0.01;
  pfvalues[2]=bla;
  
//     #type 3: muon
//     PF_EtResType3 = cms.vdouble(0.05,0,0),
//     PF_PhiResType3 = cms.vdouble(0.002),
  bla[2]=0.007;
  bla[1]=0;
  bla[0]=0;
  bla[3]=0.01;
  pfvalues[3]=bla;
  
//     # type 4: EM witout track (photon)
//     PF_EtResType4 = cms.vdouble(0.2,0.03,0.005),
//     PF_PhiResType4 = cms.vdouble(0.005),  
//   bla[0]=0.2;
//   bla[1]=0.03;
//   bla[2]=0.005;
  bla[0]=0.3;
  bla[1]=0;
  bla[2]=0;
  bla[3]=0.005;
  pfvalues[4]=bla;
//     # type 5: hadron without track (all calorimeter)
//     PF_EtResType5 = cms.vdouble(0.,1.22,0.05),
//     PF_PhiResType5 = cms.vdouble(0.02511),
//   bla[0]=0.;
//   bla[1]=1.22;
//   bla[2]=0.05;
  bla[3]=0.02511;
  // following numbers from J. Ballin's thesis
  bla[0]=0.;
  bla[1]=0.93;
  bla[2]=0.15;
  pfvalues[5]=bla;
//     # type 6: hadron without track (Forward HCAL)
//     PF_EtResType6 = cms.vdouble(0.,1.22,0.05),
//     PF_PhiResType6 = cms.vdouble(0.02511),
//   bla[0]=0.;
//   bla[1]=1.22;
//  bla[2]=0.05;
  bla[3]=0.02511;
  // following numbers from J. Ballin's thesis
  bla[0]=0.;
  bla[1]=0.93;
  bla[2]=0.15;
  pfvalues[6]=bla;
//     # type 7: EM without track (Forward HCAL)
//     PF_EtResType7 = cms.vdouble(0.,1.22,0.05),
//     PF_PhiResType7 = cms.vdouble(0.02511)  
  bla[0]=0.;
  bla[1]=1.22;
  bla[2]=0.05;
  bla[3]=0.02511;
  pfvalues[7]=bla;
  
  // some more worker variables:
  int pftype;
  double par[3],et,e,det,phi,dphi;

  CaloTowerCollection::const_iterator ct;
  reco::TrackBase::TrackQuality QUALITY= reco::TrackBase::qualityByName("highPurity");
  for( ev.toBegin();
       ! ev.atEnd() && (iev<nevmax || nevmax<1);
       ++ev) {
    iev++;
    
    if(iev%repeat==0)
      std::cout << iev << " - run " <<ev.id().run()<< ", evt " << ev.id().event() <<  std::endl;
    if(iev/repeat>=10 && repeat < 10000)
      repeat*=10;

    if(ev.id().run()>100){
      ismc=false;
      if( std::find(goodruns.begin(),goodruns.end(),ev.id().run())==goodruns.end() ){
	if(std::find(badruns.begin(),badruns.end(),ev.id().run())==badruns.end() ){
	  badruns.push_back(ev.id().run());
	  //	  std::cout << "run " << ev.id().run() << " is not in good run list" << std::endl;
	}
	continue;
      }
    }
    else
      ismc=true;
    // and check for bad events:
    currentevent.first=ev.id().run();
    currentevent.second=ev.id().event();
    if(std::find(badevents.begin(),badevents.end(),currentevent)!=badevents.end()){ 
      //      std::cout << "run " << ev.id().run() << ", event " << ev.id().event() << " is in bad event list. " << std::endl;
      continue;
    
    }
    //    try {
    //handle_techbits.getByLabel(ev,"gtDigis");
    //    }
    //    catch (cms.Exception &x){
    //      std::cout << std::endl << "Failed with cms::Exception:" << std::endl;
    //      std::cout << x.what() << std::endl;
    //    }
    
    if(iev<10)
      std::cout << "getting met..." << std::endl;
    //    handle_mets.getByLabel(ev,"met","","METSIGNIF");
    //    handle_pfmets.getByLabel(ev,"pfMet");
    //    handle_pfcands.getByLabel(ev,"particleFlow");
    //    handle_mets2.getByLabel(ev,"met","","RERECO");
    if(iev<10)
      std::cout << "getting offlinePrimaryVertices..." << std::endl;
    handle_vtxs.getByLabel(ev,"offlinePrimaryVertices");
    if(iev<10)
      std::cout << "getting tracks... " << std::endl;
    handle_tracks.getByLabel(ev,"generalTracks");
    //    std::cout << "get towers..." << std::endl;
    //    handle_ct.getByLabel(ev,"towerMaker");
    //    std::cout << handle_techbits.ptr()->technicalTriggerWord().at(0) << std::endl;
    if(iev<10)
      std::cout << "checking techinical bits.." << std::endl;
    // make sure bit 0 is on.
    if(!ismc){
      if(handle_techbits.ptr()->technicalTriggerWord().at(0)==0 && !ismc)
	continue;
      // and the beam halo bits are off...
      if(handle_techbits.ptr()->technicalTriggerWord().at(36)==1 && !ismc)
	continue;
      if(handle_techbits.ptr()->technicalTriggerWord().at(37)==1 && !ismc)
	continue;
      if(handle_techbits.ptr()->technicalTriggerWord().at(38)==1 && !ismc)
	continue;
      if(handle_techbits.ptr()->technicalTriggerWord().at(39)==1 && !ismc)
	continue;
      
      // check that bits 40 or 41 are on...
      if(handle_techbits.ptr()->technicalTriggerWord().at(40)==0 && 
	 handle_techbits.ptr()->technicalTriggerWord().at(41)==0)
	continue;
    }
    //    for(int ttr=0; ttr<50; ttr++)
    //      std::cout << handle_techbits.ptr()->technicalTriggerWord().at(ttr) << ",";
    //    std::cout << std::endl;
    if(handle_mets.ptr()->size()<1){
      std::cout << "strange, no met!" << std::endl;
      continue;
    }
    if(handle_pfmets.ptr()!=0){
      if(handle_pfmets.ptr()->size()<1){
	std::cout << "strange, no pfmet!" << std::endl;
	continue;
      }
    }
    else
      continue;
    if(iev<10)
      std::cout << "checking mets and vtxs" << std::endl;
    if(handle_vtxs.ptr()!=0){
      if(handle_vtxs.ptr()->size()<1){
	std::cout << "strange, no pv!" << std::endl;
	continue;
      }
    }
    else
      continue;
    if(iev<10)
      std::cout << "getting vtx, mets" << std::endl;
    vtxs = handle_vtxs.ptr()->begin();
    mets = handle_mets.ptr()->begin();
    pfmets = handle_pfmets.ptr()->begin();

    if(vtxs->isFake())
      continue;
    if(fabs(vtxs->z())>15){
      continue;      
    }
    if(vtxs->ndof()<5)
      continue;
    
    trackfrac=0;
    for( trk=handle_tracks.ptr()->begin(); trk!=handle_tracks.ptr()->end(); ++trk){
      if(trk->quality(QUALITY))
	trackfrac++;
    }

    if(handle_tracks.ptr()->size()>0)
       trackfrac/=(double)handle_tracks.ptr()->size();
    //    h_fractrackgood->Fill(trackfrac);
    if(trackfrac<0.25 && handle_tracks.ptr()->size()>=10)
      continue;

    container.met=mets->et();
    container.sig=mets->significance();
    container.sumEt=mets->sumEt();
    container.phi=mets->phi();
    pfcontainer.met=pfmets->et();
    pfcontainer.sig=pfmets->significance();
    pfcontainer.phi=pfmets->phi();
    pfcontainer.sumEt=pfmets->sumEt();
    //  if(!ismc || pfcontainer.sig<0){
    if(1){// always
      //      std::cout << pfmets->significance() << std::endl;
      for(pfc=handle_pfcands.ptr()->begin(); pfc!=handle_pfcands.ptr()->end(); pfc++){
	//	std::cout << pfc->phi() << " " << pfc->et() << " " << pfc->particleId() << std::endl;
	pftype = pfc->particleId();
	et=pfc->et();
	e=pfc->energy();

// 	if(pftype<=3){
// 	  e=pfc->p();
// 	  et=pfc->pt();
// 	}
	par[0]=pfvalues[pftype][0];
	par[1]=pfvalues[pftype][1];
	par[2]=pfvalues[pftype][2];
	
	det=et*sqrt((par[2]*par[2])+(par[1]*par[1]/e)+(par[0]*par[0]/(e*e)));
	phi = pfc->phi();
	dphi = pfvalues[pftype][3]*et;
	if(pftype<=3){
	  reco::TrackRef trk = pfc->trackRef();
	  if(pfc->trackRef().isAvailable()){
	    det=pfc->trackRef()->ptError();
	    dphi=pfc->trackRef()->phiError();
	    //	    std::cout << "TRK : " << pfc->trackRef()->ptError() << " " << pfc->trackRef()->phiError() << std::endl;
	  }
	  else if(pfc->gsfTrackRef().isAvailable()){
	    std::cout <<"GSF: " <<  pfc->gsfTrackRef()->ptError() << " " << pfc->gsfTrackRef()->phiError() << std::endl;
	    det= pfc->gsfTrackRef()->ptError();
	    dphi=pfc->gsfTrackRef()->phiError();
	  }
	}
	pfccontainer.p=pfc->p();
	pfccontainer.pt=pfc->pt();
	pfccontainer.e=pfc->energy();
	pfccontainer.et=pfc->et();
	pfccontainer.phi=pfc->phi();
	pfccontainer.eta=pfc->eta();
	pfccontainer.ehad=pfc->hcalEnergy();
	pfccontainer.eem=pfc->ecalEnergy();
	pfccontainer.de=det;
	pfccontainer.dphi=dphi;
	pfccontainer.pid=pftype;
	hftree->Fill();
	//	std::cout << et << " "<< det << " " << phi << " " << dphi << std::endl;
	metsig::SigInputObj inobj(emptystring,et,phi,det,dphi);
	EventVec1.push_back(inobj);
	metsig::SigInputObj inobj2(emptystring,et,phi,det,dphi);
	if(pftype<=3)
	  EventVec2.push_back(inobj2);
	else
	  EventVec3.push_back(inobj2);
      }
      metsig::significanceAlgo algo1,algo2,algo3;
      algo1.addObjects(EventVec1);
      algo2.addObjects(EventVec2);
      algo3.addObjects(EventVec3);
      double workmet=0;//pfcontainer.met;
      double workphi=0;//pfcontainer.phi;
      double workset=0;//.qpfcontainer.sumEt;
      //      std::cout << workmet << " " << workphi << " " << workset << " " << pfcontainer.sig << std::endl;
      pfcontainer.sig=algo1.significance(workmet,workphi,workset);
      pfcontainer.met=workmet;
      pfcontainer.phi=workphi;
      pfcontainer.sumEt=workset;
      //      std::cout << workmet << " " << workphi << " " << workset << " " << pfcontainer.sig << std::endl;
      workmet=workphi=workset=0;
      pftrkcontainer.sig=algo2.significance(workmet,workphi,workset);
      pftrkcontainer.met=workmet;
      pftrkcontainer.phi=workphi;
      pftrkcontainer.sumEt=workset;
      
      workmet=workphi=workset=0;
      pfcalcontainer.sig=algo3.significance(workmet,workphi,workset);
      pfcalcontainer.met=workmet;
      pfcalcontainer.phi=workphi;
      pfcalcontainer.sumEt=workset;
      EventVec1.clear();
      EventVec1.resize(0);
      EventVec2.clear();
      EventVec2.resize(0);
      EventVec3.clear();
      EventVec3.resize(0);
    }
    evcontainer.runnum=ev.id().run();
    evcontainer.evnum=ev.id().event();
    littletree->Fill();
    //    h_met->Fill(mets->et());
    //    h_metsignif->Fill(mets->significance());
    //    h_metsum->Fill(mets->sumEt());
    //    h_metsignfsumEt->Fill(mets->sumEt(),mets->significance());
    //    h_metsignfsqrtsumEt->Fill(mets->et()/sqrt(mets->sumEt()),mets->significance());
    //    std::cout << mets->significance() << " " << mets->sumEt() << std::endl;
  }// end of event loop.

//   TCanvas *workcanv = new TCanvas();
//   workcanv->cd();
//   h_vtx_z->Draw();
//   workcanv->Modified();
//   //  workcanv->Print("vertex_z.gif");
//   workcanv->SetLogy();
//   h_metsignif->Draw();
//   workcanv->Modified();
//   //  workcanv->Print("signif.gif");
//   h_met->Draw();
//   workcanv->Modified();
//   //  workcanv->Print("basicmet.gif");
//   workcanv->SetLogy(0);
//   h_metsignfsqrtsumEt->Draw("col");
//   workcanv->Modified();
//   //  workcanv->Print("sign_vs_stdsign.gif");
  for(ii=0; ii<badruns.size(); ++ii)
    std::cout << "additional bad run: " << badruns[ii] << std::endl;
  outputfile->Save();
  outputfile->Write();
  outputfile->Close();
}
