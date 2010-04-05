// -*- C++ -*-
//
// Package:    HCALDataIntegrityTools
// Class:      HCALDataIntegrityTools
// 
/**\class HCALDataIntegrityTools HCALDataIntegrityTools.cc JSturdy/HCALDataIntegrityTools/src/HCALDataIntegrityTools.cc

 Description: [Set of basic HCAL diagnostics tools]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jared Todd Sturdy,28 S-027,+41227676019,
//         Created:  Mon Apr  5 14:24:58 CEST 2010
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

#include "JSturdy/HCALDataIntegrityTools/interface/HCALDataIntegrityTools.h"

HCALDataIntegrityTools::HCALDataIntegrityTools(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


HCALDataIntegrityTools::~HCALDataIntegrityTools()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HCALDataIntegrityTools::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
HCALDataIntegrityTools::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HCALDataIntegrityTools::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HCALDataIntegrityTools);
