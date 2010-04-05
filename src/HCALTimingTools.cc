// -*- C++ -*-
//
// Package:    HCALTimingTools
// Class:      HCALTimingTools
// 
/**\class HCALTimingTools HCALTimingTools.cc JSturdy/HCALTimingTools/src/HCALTimingTools.cc

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

#include "JSturdy/HCALTimingTools/interface/HCALTimingTools.h"

HCALTimingTools::HCALTimingTools(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


HCALTimingTools::~HCALTimingTools()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HCALTimingTools::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

/* Find data that doesn't meet our expected BCN 
 * Determine if there are deposites in other sub detectors
 * 
 *
 */

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
HCALTimingTools::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HCALTimingTools::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HCALTimingTools);
