#include "JSturdy/HcalAnalysisTools/interface/HCALRBXAnalyzer.h"
// define sizes of ieta arrays for each subdetector

#define PI          3.1415926535897932

using namespace std;
using namespace edm;

/*  

v1.0
20 November 2008
by Jared Sturdy

*/

HCALRBXAnalyzer::HCALRBXAnalyzer()
{}

// destructor
HCALRBXAnalyzer::~HCALRBXAnalyzer() {}

void HCALRBXAnalyzer::setup(const edm::ParameterSet& ps, DQMStore* dbe)
{
  HcalLogicalMapGenerator mygen;
  theLogicalMap = mygen.createMap();

  //  ievt_=0; // event counter
  firstFED_ = FEDNumbering::getHcalFEDIds().first;
  //cout <<"FIRST FED = "<<firstFED_<<endl;

  for (int i=FEDNumbering::getHcalFEDIds().first; 
       i<=FEDNumbering::getHcalFEDIds().second; ++i) 
    {
      fedUnpackList_.push_back(i);
    } // for (int i=FEDNumbering::getHcalFEDIds().first;...



  char name[128];
  for (int ii = 0; ii < NUMFEDS; ii++) {
  }

  return;

} // void HCALRBXAnalyzer::setup()

void HCALRBXAnalyzer::analye(edm::Event const&e, edm::EventSetup const&s) 
{
  // try to get digis
  edm::Handle<HBHEDigiCollection> hbhe_digi;
  edm::Handle<HODigiCollection>   ho_digi;
  edm::Handle<HFDigiCollection>   hf_digi;

  if (!(e.getByLabel(digiLabel_,hbhe_digi)))
    {
      edm::LogWarning("HCALRBXAnalyzer")<< digiLabel_<<" hbhe_digi not available";
      return;
    }
  
  if (!(e.getByLabel(digiLabel_,hf_digi)))
    {
      edm::LogWarning("HCALRBXAnalyzer")<< digiLabel_<<" hf_digi not available";
      return;
    }
  if (!(e.getByLabel(digiLabel_,ho_digi)))
    {
      edm::LogWarning("HCALRBXAnalyzer")<< digiLabel_<<" ho_digi not available";
      return;
    }

  c.get<HcalDbRecord>().get(conditions_);
  const HcalQIEShape*  shape = conditions_->getHcalShape();
  HcalCalibrations     calibs=conditions_->getHcalCalibrations(*chan);
  const HcalQIECoder*  channelCoder = conditions_->getHcalCoder(*chan);

  edm::Handle<FEDRawDataCollection> rawraw;
  if (!(e.getByLabel(FEDRawDataCollection_,rawraw)))
    {
      edm::LogWarning("HcalRawDataMonitor")<<" raw data with label "<<FEDRawDataCollection_ <<" not available";
      return;
    }

  edm::Handle<HcalUnpackerReport> report;  
  if (!(e.getByLabel(digiLabel_,report)))
    {
      edm::LogWarning("HCALRBXAnalyzer")<< digiLabel_<<" unpacker report not available";
      return;
    }

  DetIdVector::const_iterator dummy = report.bad_quality_begin();
  DetIdVector::const_iterator dumm2 = report.bad_quality_end();
  
  //std::cout<<"unpacker digi at event: "<<eventNum<<std::endl;
  
  for ( DetIdVector::const_iterator baddigi_iter=report.bad_quality_begin(); 
	baddigi_iter != report.bad_quality_end();
	++baddigi_iter)
    {
      HcalDetId id(baddigi_iter->rawId());
      int rDepth = id.depth();
      int rPhi   = id.iphi();
      int rEta   = id.ieta();
      int binEta = CalcEtaBin(id.subdet(), rEta, rDepth); // why is this here?
      std::cout<<"Bad unpacker digi at event/LS: "<<eventNum<<"/"<<lumiBlock<<" in "<<theLogicalMap.getFrontEndId(id).rbx()<<std::endl;
      std::cout<<"iEta: "<<rEta<<"  iPhi"<<rPhi<<std::endl;
      if (id.subdet()==HcalBarrel) ++hbHists.count_bad;
      else if (id.subdet()==HcalEndcap) ++heHists.count_bad;
      else if (id.subdet()==HcalForward) 
      else if (id.subdet()==HcalOuter) 
      else 
	continue; // skip anything that isn't HB, HE, HO, HF
      // extra protection against nonsensical values -- prevents occasional crashes
    //if (binEta < 85 && binEta >= 0 
    //	  && (rPhi-1) >= 0 && (rPhi-1)<72 
    //	  && (rDepth-1) >= 0 && (rDepth-1)<4)
    //	{
    //	  ++badunpackerreport[binEta][rPhi-1][rDepth-1];
    //	  ++baddigis[binEta][rPhi-1][rDepth-1];  
    //	  //if (binEta < 30 && binEta > 15 && (rPhi-1) < 63 && (rPhi-1) > 58)
    //	  //if (rEta < 31 && rEta > 14 && (rPhi-1) < 64 && (rPhi-1) > 56)
    //
    //	}
      
    }
  // all objects grabbed; event is good
  if (debug_>1) std::cout <<"\t<HCALRBXAnalyzer::analyze>  Processing good event! event # = "<<ievt_<<std::endl;
  processEvent(*rawraw, *hbhe_digi, *ho_digi, *hf_digi, *cond, *report, *emap);
}
void HCALRBXAnalyzer::processEvent(const FEDRawDataCollection& rawraw,
				   const HBHEDigiCollection&   hbhe,
				   const HODigiCollection&     ho,
				   const HFDigiCollection&     hf,
				   const HcalDbService&        cond,
				   const HcalUnpackerReport&   report,
				   const HcalElectronicsMap&   emap
				   )
  
{
  process_RBX(rawraw, hbhe, ho, hf, cond, report, emap);
  return;
} // void HCALRBXAnalyzer::processEvent




void HCALRBXAnalyzer::process_RBX(const FEDRawDataCollection& rawraw,
				  const HBHEDigiCollection&   hbhe,
				  const HODigiCollection&     ho,
				  const HFDigiCollection&     hf,
				  const HcalDbService&        cond,
				  const HcalUnpackerReport&   report,
				  const HcalElectronicsMap&   emap
				  )
{

}

void HCALRBXAnalyzer::process_raw(const FEDRawDataCollection& rawraw,
				  const HcalDbService&        cond,
				  const HcalUnpackerReport&   report,
				  const HcalElectronicsMap&   emap
				  )
{
  /*
    This processes Raw Data.
    Additional info on working with Raw Data can be found in 
    HcalDataFormatMonitor.cc
  */
  


  // Loop over all FEDs reporting the event, unpacking if good.
  for (vector<int>::const_iterator i=fedUnpackList_.begin();i!=fedUnpackList_.end(); i++) 
    {
      const FEDRawData& fed = rawraw.FEDData(*i);
      if (fed.size()<12) continue; // Was 16.
      unpack(fed,emap);
    } // for (vector<int>::const_iterator i=fedUnpackList_.begin();...
  
  return;
  
} // void HCALRBXAnalyzer::process_RBX(const FEDRawDataCollection& rawraw,


// Process one FED's worth (one DCC's worth) of the event data.
void HCALRBXAnalyzer::unpack(const FEDRawData& raw, 
			     const HcalElectronicsMap& emap)
{
  /* 
   * This unpacks raw data info.  Additional info on working with the raw 
   * data can be found in the unpack method of HcalDataFormatMonitor.cc
  */
  

  // get the DCC header
  const HcalDCCHeader* dccHeader=(const HcalDCCHeader*)(raw.data());
  if(!dccHeader) return;
  
  // get the DCC trailer 
  unsigned char* trailer_ptr = (unsigned char*) (raw.data()+raw.size()-sizeof(uint64_t));
  FEDTrailer trailer = FEDTrailer(trailer_ptr);
  
  //DCC Event Fragment sizes distribution, in bytes.
  //  int rawsize=raw.size();
  
  int dccid=dccHeader->getSourceId();
  
  HcalHTRData htr;  
  for (int spigot=0; spigot<HcalDCCHeader::SPIGOT_COUNT; spigot++) {    
    if (!dccHeader->getSpigotPresent(spigot)) continue;
    activePlots[dccid-700][spigot]=true;
    
    // Load the given decoder with the pointer and length from this spigot.
    dccHeader->getSpigotData(spigot,htr, raw.size()); 
    const unsigned short* HTRraw = htr.getRawData();
    unsigned short HTRwdcount = HTRraw[htr.getRawLength() - 2];
    
    //fix me!
    HTRwdcount=htr.getRawLength();
    
    // Size checks for internal consistency
    // getNTP(), get NDD() seems to be mismatched with format. Manually:
    int NTP = ((htr.getExtHdr6() >> 8) & 0x00FF);
    int NDAQ = (HTRraw[htr.getRawLength() - 4] & 0x7FF);
    
    if ( !  ((HTRwdcount != 8)               ||
	     (HTRwdcount != 12 + NTP + NDAQ) ||
	     (HTRwdcount != 20 + NTP + NDAQ)    )) {
      //incompatible Sizes declared. Skip it.
      continue; }
    bool EE = ((dccHeader->getSpigotErrorBits(spigot) >> 2) & 0x01);
    if (EE) { 
      if (HTRwdcount != 8) {
	//incompatible Sizes declared. Skip it.
	continue; 
      }
    }
    else{ //For non-EE,
      if ((HTRwdcount-NDAQ-NTP) != 20) {
	//incompatible Sizes declared. Skip it.
	continue;
      }
    }
    
    if (htr.isHistogramEvent()) continue;
    
    
    // Fish out Front-End Errors from the precision channels
    const short unsigned int* daq_first, *daq_last, *tp_first, *tp_last;
    const HcalQIESample* qie_begin, *qie_end, *qie_work;
    
    // get pointers
    htr.dataPointers(&daq_first,&daq_last,&tp_first,&tp_last);
    
    qie_begin=(HcalQIESample*)daq_first;
    qie_end=(HcalQIESample*)(daq_last+1); // one beyond last..
    
    int lastcapid=-1;
    int lastfibchan =0, samplecounter=0;
    int channum=0; // Valid: [1,24]
    double weightedval=0.0;
    double charge=0.0;
    
    // Loop over DAQ words for this spigot
    for (qie_work=qie_begin; qie_work!=qie_end; ) {
      if (qie_work->raw()==0xFFFF) {
	qie_work++; // filler word
	continue;
      }
      // Beginning digitized hit of this Half-HTR?
      if (qie_work==qie_begin) {
	channum= (3*qie_work->fiber()) - 2 + qie_work->fiberChan();  
	lastcapid=qie_work->capid();
	samplecounter=1;
      }
      
      // or the first TS of a this channel's DAQ data?
      else if (qie_work->fiberAndChan() != lastfibchan) {
	channum= (3*qie_work->fiber()) - 2 + qie_work->fiberChan();
	if (samplecounter != htr.getNDD() ) {
	}
	samplecounter=1;
      }
      
      else { //precision samples not the first timeslice
	int hope = lastcapid +1;
	if (hope==4) hope = 0;
	if (qie_work->capid() != hope){
	}
	samplecounter++;
      }
      
      if (qie_work->dv() && !(qie_work->er())) {
	if (samplecounter==10) {
	}
	//For every precision data sample in Hcal:
      }
      // Prepare for the next round...
      lastcapid=qie_work->capid();
      lastfibchan=qie_work->fiberAndChan();
      qie_work++;
    }
    if (dccid==723 && spigot==3) { //the ZDC spigot
    }
  }
  
  uint64_t* lastDataWord = (uint64_t*) ( raw.data()+raw.size()-(2*sizeof(uint64_t)) );
  int EvFragLength = ((*lastDataWord>>32)*8);
  EvFragLength = raw.size();
  
  
  if (!dccHeader->thereIsASecondCDFHeaderWord()) 
    {
      cout <<"No second CDF header found!"<<endl;
    }
  
  return;
} // void HCALRBXAnalyzer::unpack(...)

void HCALRBXAnalyzer::endOfJob()
{

}
