#ifndef HCALRBXANALYZER_H
#define HCALRBXANALYZER_H

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUnpacker.h"
// The following are needed for using pedestals in fC:
#include "CondFormats/HcalObjects/interface/HcalPedestal.h"
#include "CondFormats/HcalObjects/interface/HcalPedestalWidth.h"
#include "CondFormats/HcalObjects/interface/HcalLogicalMap.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalLogicalMapGenerator.h"

// Raw data stuff
#include "EventFilter/HcalRawToDigi/interface/HcalUnpacker.h"
#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"

// Use for stringstream
#include <iostream>
#include <iomanip>
#include <cmath>

#define NUMFEDS     32
#define NUMSPIGOTS  15
#define NUMCHANNELS 24
#define NUMFIBRES    8
#define NUMTS       10

/** \class HCALRBXAnalyzer
  *
  * $Date: 2010/04/30 08:05:12 $
  * $Revision: 1.1 $
  * \author J. Sturdy - Univ. of California
  */

class HCALRBXAnalyzer:  public edm::EDAnalyzer {
 public:
  HCALRBXAnalyzer();
  ~HCALRBXAnalyzer();
  
  void unpack(const FEDRawData& raw, const HcalElectronicsMap& emap);
  void setup(const edm::ParameterSet& ps, TTree* thetree);
  void reset();
  void clearME();
  void endOfJob();
 
  // processEvent routine -- specifies what inputs are looked at each event
  void processEvent(const FEDRawDataCollection& rawraw,
		    const HBHEDigiCollection&   hbhe,
		    const HODigiCollection&     ho,
		    const HFDigiCollection&     hf,
		    const HcalDbService&        cond,
		    const HcalUnpackerReport&   report,
		    const HcalElectronicsMap&   emap
		    );
  
  void process_RBX(const FEDRawDataCollection& rawraw,
		   const HBHEDigiCollection&   hbhe,
		   const HODigiCollection&     ho,
		   const HFDigiCollection&     hf,
		   const HcalUnpackerReport&   report,
		   const HcalElectronicsMap&   emap
		   );

  void process_raw(const FEDRawDataCollection& rawraw,
		   const HcalDbService&      cond,
		   const HcalUnpackerReport&   report,
		   const HcalElectronicsMap&   emap
		   );

  void process_digi(const HBHEDigiCollection& hbhe,
		    const HODigiCollection&   ho,
		    const HFDigiCollection&   hf,
		    const HcalDbService&      cond,
		    const HcalUnpackerReport& report,
		    const HcalElectronicsMap& emap
		    );

 private:
  
  //  int ievt_;
  //  MonitorElement* meEVT_;

  HcalLogicalMap theLogicalMap;  
  std::vector <int> fedUnpackList_;
  int firstFED_;
  // Add in histograms here (MonitorElements can handle TH1F, TH2F, TProfile plots)

}; // class HCALRBXAnalyzer

#endif  
