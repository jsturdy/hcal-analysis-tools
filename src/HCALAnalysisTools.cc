// -*- C++ -*-
//
// Package:    HCALAnalysisTools
// Class:      HCALAnalysisTools
// 
/**\class HCALAnalysisTools HCALAnalysisTools.cc JSturdy/HCALAnalysisTools/src/HCALAnalysisTools.cc

 Description: [one line class summary]

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
//
// class declaration
//

class HCALAnalysisTools : public edm::EDAnalyzer {
   public:
      explicit HCALAnalysisTools(const edm::ParameterSet&);
      ~HCALAnalysisTools();


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
HCALAnalysisTools::HCALAnalysisTools(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


HCALAnalysisTools::~HCALAnalysisTools()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HCALAnalysisTools::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
HCALAnalysisTools::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HCALAnalysisTools::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HCALAnalysisTools);
