#ifndef GUARD_DQM_HCALMONITORTASKS_HCALLATENCYMONITOR_H
#define GUARD_DQM_HCALMONITORTASKS_HCALLATENCYMONITOR_H

#include "DQM/HcalMonitorTasks/interface/HcalBaseMonitor.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUnpacker.h"
// The following are needed for using pedestals in fC:
#include "CondFormats/HcalObjects/interface/HcalPedestal.h"
#include "CondFormats/HcalObjects/interface/HcalPedestalWidth.h"

// Raw data stuff
#include "DQM/HcalMonitorTasks/interface/HcalBaseMonitor.h"
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

/** \class HcalLatencyMonitor
  *
  * $Date: 2008/11/20 14:57:01 $
  * $Revision: 1.1 $
  * \author J. Sturdy - Univ. of California
  */

class HcalLatencyMonitor:  public HcalBaseMonitor {
 public:
  HcalLatencyMonitor();
  ~HcalLatencyMonitor();
  
  void unpack(const FEDRawData& raw, const HcalElectronicsMap& emap);
  void setup(const edm::ParameterSet& ps, DQMStore* dbe);
  void reset();
  void clearME();
  void endOfJob();
 
  // processEvent routine -- specifies what inputs are looked at each event
  void processEvent(const FEDRawDataCollection& rawraw,
		    const HcalUnpackerReport& report,
		    const HcalElectronicsMap& emap
		    );
  
  void process_LED_Latency(const FEDRawDataCollection& rawraw,
			   const HcalUnpackerReport& report,
			   const HcalElectronicsMap& emap);
  
 private:
  
  //  int ievt_;
  //  MonitorElement* meEVT_;
  
  std::vector <int> fedUnpackList_;
  int firstFED_;
  // Add in histograms here (MonitorElements can handle TH1F, TH2F, TProfile plots)

  MonitorElement* LEDLatency_;
  MonitorElement* LEDChargeRMS_;    //distribution of the RMS of average charge in each time bin
  MonitorElement* LEDTimingRMS_;    

  MonitorElement* LEDBimodal_;    
  MonitorElement* LEDBimodalBins_;    
  MonitorElement* LEDSteadyMinoritySpectrum_;
  MonitorElement* LEDFibreErrors_;

  MonitorElement* LEDCharge_[NUMFEDS];       //distribution of average charge in each time bin
  MonitorElement* LEDEventTiming_[NUMFEDS];  //distribution of charge weighted time in each event - binned per 0.1 TS
  MonitorElement* LEDBCNIdle_[NUMFEDS];      //idle-bcn

}; // class HcalLatencyMonitor

#endif  
