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
#include <string>
#include "DataFormats/Math/interface/deltaR.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
// INCLUDE ALL CMSSW FORMATS HERE!
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/MET.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"

#endif // end of INCLUDE ALL CMSSW FORMATS HERE!

void makeplot_fwlite(){
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
 
 // good runs:
  std::vector<uint64_t> goodruns;
  // got these from Artur
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
  ifstream filelist("filelist.txt");
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
      
      if(tempname.size()>0){
	gSystem->Exec(command.Data());
	fileNames.push_back(tempname);
      }
    }
  }

  std::cout << "reading " << fileNames.size() << " files." << std::endl;
  filelist.close();
 
  fwlite::ChainEvent ev(fileNames);  

  fwlite::Handle<L1GlobalTriggerReadoutRecord>  handle_techbits;
  fwlite::Handle<std::vector<reco::CaloMET> > handle_mets;
  fwlite::Handle<std::vector<reco::Vertex> > handle_vtxs;
  std::vector<reco::Vertex> vtxs;

  TFile *outputfile = new TFile("outputhistos.root","recreate");
  outputfile->cd();
  TH1D *h_vtx_z = new TH1D("h_vtx_z","",100,-20,20);
  h_vtx_z->SetXTitle("PV z (cm)");
  h_vtx_z->SetYTitle("events");
  
  TH1D *h_vtx_ndof = new TH1D("h_vtx_ndof","",100,0,100);
  h_vtx_ndof->SetXTitle("PV ndof");
  h_vtx_ndof->SetYTitle("events");

  TH1D *h_met = new TH1D("h_met","",100,0,10);
  h_met->SetXTitle("missing E_{T} (GeV)");
  h_met->SetYTitle("events");
  //  outputfile->Add(h_met);

  TH1D *h_metsignif = new TH1D("h_metsignif","",104,-2,50);
  h_metsignif->SetXTitle("missing E_{T} significance");
  h_metsignif->SetYTitle("events");

  TH1D *h_metsum = new TH1D("h_metsum","",100,0,100);
  h_metsum->SetXTitle("#Sigma E_{T}");
  h_metsum->SetYTitle("events");

  TH2D *h_metsignfsumEt = new TH2D("h_metsignfsumEt","",100,0,100,100,0,20);
  h_metsignfsumEt->SetXTitle("#Sigma E_{T}");
  h_metsignfsumEt->SetYTitle("missing E_{T} significance");
  TH2D *h_metsignfsqrtsumEt = new TH2D("h_metsignfsqrtsumEt","",100,0,5,100,0,20);
  
  h_metsignfsqrtsumEt->SetXTitle("missing E_{T}/sqrt(#Sigma E_{T})");
  h_metsignfsqrtsumEt->SetYTitle("missing E_{T} significance");

  TTree *littletree = new TTree("t","t");
  struct stuff{
    Float_t met;
    Float_t sig;
    Float_t phi;
    Float_t sumEt;
    Int_t runnum;
  };
  stuff container;

  littletree->Branch("met",&container,"met/F:sig/F:phi/F:sumEt/F:runnum/I");

  size_t iev=0;
  size_t nevmax=-1;
  size_t ii,jj,kk;
  std::pair<uint64_t,uint64_t> currentevent;

  for( ev.toBegin();
       ! ev.atEnd() && (iev<nevmax || nevmax<1);
       ++ev) {
    iev++;
    
    if(iev%100==0)
      std::cout << iev << " - run " <<ev.id().run()<< ", evt " << ev.id().event() <<  std::endl;

    if( std::find(goodruns.begin(),goodruns.end(),ev.id().run())==goodruns.end() && ev.id().run()>100){
     std::cout << "run " << ev.id().run() << " is not in good run list" << std::endl;
     continue;
    }
    // and check for bad events:
    currentevent.first=ev.id().run();
    currentevent.second=ev.id().event();
    if(std::find(badevents.begin(),badevents.end(),currentevent)!=badevents.end()){ 
      std::cout << "run " << ev.id().run() << ", event " << ev.id().event() << " is in bad event list. " << std::endl;
      continue;
    
    }
    handle_techbits.getByLabel(ev,"gtDigis");
    handle_mets.getByLabel(ev,"met");
    handle_vtxs.getByLabel(ev,"offlinePrimaryVertices");
     
    //    std::cout << handle_techbits.ptr()->technicalTriggerWord().at(0) << std::endl;
    
    // make sure bit 0 is on.
    if(handle_techbits.ptr()->technicalTriggerWord().at(0)==0)
      continue;

    //    for(int ttr=0; ttr<50; ttr++)
    //      std::cout << handle_techbits.ptr()->technicalTriggerWord().at(ttr) << ",";
    //    std::cout << std::endl;
    if(handle_mets.ptr()->size()<1)
      continue;
    if(handle_vtxs.ptr()->size()<1)
      continue;
    
    std::vector<reco::Vertex>::const_iterator vtxs = handle_vtxs.ptr()->begin();
    std::vector<reco::CaloMET>::const_iterator mets = handle_mets.ptr()->begin();
    
    //    std::cout << mets->et() << " " << vtxs->ndof() << " " << vtxs->isFake() << " " << vtxs->z() <<  std::endl;
    if(vtxs->isFake())
      continue;
    h_vtx_z->Fill(vtxs->z());
    
    if(fabs(vtxs->z())>15){
      std::cout << "PV(z) at " << vtxs->z() << ", continueing" << std::endl;
      continue;      
    }
    h_vtx_ndof->Fill(vtxs->ndof());
    if(vtxs->ndof()<5)
      continue;
    container.met=mets->et();
    container.sig=mets->significance();
    container.sumEt=mets->sumEt();
    container.phi=mets->phi();
    container.runnum=ev.id().run();
    littletree->Fill();
    h_met->Fill(mets->et());
    h_metsignif->Fill(mets->significance());
    h_metsum->Fill(mets->sumEt());
    h_metsignfsumEt->Fill(mets->sumEt(),mets->significance());
    h_metsignfsqrtsumEt->Fill(mets->et()/sqrt(mets->sumEt()),mets->significance());
    //    std::cout << mets->significance() << " " << mets->sumEt() << std::endl;
  }// end of event loop.

  TCanvas *workcanv = new TCanvas();
  workcanv->cd();
  h_vtx_z->Draw();
  workcanv->Modified();
  workcanv->Print("vertex_z.gif");
  workcanv->SetLogy();
  h_metsignif->Draw();
  workcanv->Modified();
  workcanv->Print("signif.gif");
  h_met->Draw();
  workcanv->Modified();
  workcanv->Print("basicmet.gif");
  workcanv->SetLogy(0);
  h_metsignfsqrtsumEt->Draw("col");
  workcanv->Modified();
  workcanv->Print("sign_vs_stdsign.gif");
  
  outputfile->Save();
  outputfile->Write();
  outputfile->Close();
}
