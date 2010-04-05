// -*- C++ -*-
//
// Package:    HCALAnalysisTools
// Class:      HCALAnalysisTools
// 
/**\class HCALAnalysisTools HCALAnalysisTools.cc JSturdy/HCALAnalysisTools/src/HCALAnalysisTools.cc

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


#endif
