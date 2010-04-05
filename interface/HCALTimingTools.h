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
//
// class declaration
//

class HCALTimingTools : public edm::EDAnalyzer {
  public:
    explicit HCALTimingTools(const edm::ParameterSet&);
    ~HCALTimingTools();
  
  
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
