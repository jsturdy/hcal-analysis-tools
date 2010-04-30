#include "DQM/HcalMonitorTasks/interface/HcalLatencyMonitor.h"
// define sizes of ieta arrays for each subdetector

#define PI          3.1415926535897932

using namespace std;
using namespace edm;

/*  

v1.0
20 November 2008
by Jared Sturdy

*/

double LEDLatencySum[NUMFEDS][NUMSPIGOTS][NUMCHANNELS]     = { 0.0 };
int    LEDLatencyCounter[NUMFEDS][NUMSPIGOTS][NUMCHANNELS] = { 0 };

double LEDTotalChargeSum[NUMFEDS][NUMSPIGOTS][NUMCHANNELS][NUMTS] = { 0.0 };
int     LEDChargeCounter[NUMFEDS][NUMSPIGOTS][NUMCHANNELS][NUMTS] = { 0 };

bool activePlots[NUMFEDS][NUMSPIGOTS] = {false};
// constructor
HcalLatencyMonitor::HcalLatencyMonitor()
{}

// destructor
HcalLatencyMonitor::~HcalLatencyMonitor() {}

void HcalLatencyMonitor::reset() {}

void HcalLatencyMonitor::clearME()
{

} // void HcalLatencyMonitor::clearME()


void HcalLatencyMonitor::setup(const edm::ParameterSet& ps, DQMStore* dbe)
{
  HcalBaseMonitor::setup(ps,dbe);  // perform setups of base class

  //  ievt_=0; // event counter
  baseFolder_ = rootFolder_ + "LatencyMonitor/"; // Will create an "LatencyMonitor" subfolder in .root output
  if (fVerbosity) cout <<"<HcalLatencyMonitor::setup> Setup in progress"<<endl;

  
  if(fVerbosity) cout << "About to pushback fedUnpackList_" << endl;
  firstFED_ = FEDNumbering::getHcalFEDIds().first;
  //cout <<"FIRST FED = "<<firstFED_<<endl;

  for (int i=FEDNumbering::getHcalFEDIds().first; 
       i<=FEDNumbering::getHcalFEDIds().second; ++i) 
    {
      if(fVerbosity) cout << "<HcalLatencyMonitor::setup>:Pushback for fedUnpackList_: " << i <<endl;
      fedUnpackList_.push_back(i);
    } // for (int i=FEDNumbering::getHcalFEDIds().first;...



  if (m_dbe)
    {
      m_dbe->setCurrentFolder(baseFolder_);
    
      char name[128];
      for (int ii = 0; ii < NUMFEDS; ii++) {
	sprintf(name,"LEDChargeTiming");
	m_dbe->setCurrentFolder(baseFolder_+name);//Per event
	sprintf(name,"LEDTimingDistributionFED%i",ii+700);
	LEDEventTiming_[ii] = m_dbe->book2D(name,name,(NUMSPIGOTS * (NUMCHANNELS+1)),0,(NUMSPIGOTS * (NUMCHANNELS+1)),502,-0.02,10.02);
	LEDEventTiming_[ii]->setAxisTitle("Mean LED Timing",2);
	
	sprintf(name,"LEDBCN_Idle");
	m_dbe->setCurrentFolder(baseFolder_+name);//Per event
	sprintf(name,"LEDBCNIdleFED%i",ii+700);
	LEDBCNIdle_[ii] = m_dbe->book2D(name,name,(NUMSPIGOTS * (NUMFIBRES+1)),0,(NUMSPIGOTS * (NUMFIBRES+1)),20,1089.5,1109.5);
	LEDBCNIdle_[ii]->setAxisTitle("Idle BCN",2);
      }

    } // if (m_dbe)

  return;

} // void HcalLatencyMonitor::setup()


void HcalLatencyMonitor::processEvent(const FEDRawDataCollection& rawraw,
				      const HcalUnpackerReport& report,
				      const HcalElectronicsMap& emap
				      )
  
{
  if (!m_dbe)
    {
      if (fVerbosity) cout <<"HcalLatencyMonitor::processEvent   DQMStore not instantiated!!!"<<endl;
      return;
    }

  process_LED_Latency(rawraw, report, emap);
  //printf("%12i  %12i  %7i\n",evNum,orbNum,EEcheck);
  return;
} // void HcalLatencyMonitor::processEvent




void HcalLatencyMonitor::process_LED_Latency(const FEDRawDataCollection& rawraw,
					     const HcalUnpackerReport& report,
					     const HcalElectronicsMap& emap)
{
  /*
    This processes Raw Data.
    Additional info on working with Raw Data can be found in 
    HcalDataFormatMonitor.cc
  */
  

  // Should not see this error
  if(!m_dbe) 
    {
      cout <<"HcalLatencyMonitor::process_LED_Latency:  DQMStore not instantiated!!!\n"<<endl;
      return;
    }

  // Loop over all FEDs reporting the event, unpacking if good.
  for (vector<int>::const_iterator i=fedUnpackList_.begin();i!=fedUnpackList_.end(); i++) 
    {
      const FEDRawData& fed = rawraw.FEDData(*i);
      if (fed.size()<12) continue; // Was 16.
      unpack(fed,emap);
    } // for (vector<int>::const_iterator i=fedUnpackList_.begin();...

  return;

} // void HcalLatencyMonitor::process_LED_Latency(const FEDRawDataCollection& rawraw,


// Process one FED's worth (one DCC's worth) of the event data.
void HcalLatencyMonitor::unpack(const FEDRawData& raw, 
				const HcalElectronicsMap& emap)
{
  /* 
     This unpacks raw data info.  Additional info on working with the raw data can be found in the unpack method of HcalDataFormatMonitor.cc
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
    LEDBCNIdle_[dccid-700]->Fill(((NUMFIBRES+1)*spigot)+1,htr.getFib1OrbMsgBCN());

    LEDBCNIdle_[dccid-700]->Fill(((NUMFIBRES+1)*spigot)+2,htr.getFib2OrbMsgBCN());

    LEDBCNIdle_[dccid-700]->Fill(((NUMFIBRES+1)*spigot)+3,htr.getFib3OrbMsgBCN());

    LEDBCNIdle_[dccid-700]->Fill(((NUMFIBRES+1)*spigot)+4,htr.getFib4OrbMsgBCN());

    LEDBCNIdle_[dccid-700]->Fill(((NUMFIBRES+1)*spigot)+5,htr.getFib5OrbMsgBCN());

    LEDBCNIdle_[dccid-700]->Fill(((NUMFIBRES+1)*spigot)+6,htr.getFib6OrbMsgBCN());

    LEDBCNIdle_[dccid-700]->Fill(((NUMFIBRES+1)*spigot)+7,htr.getFib7OrbMsgBCN());

    LEDBCNIdle_[dccid-700]->Fill(((NUMFIBRES+1)*spigot)+8,htr.getFib8OrbMsgBCN());

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

	weightedval+=samplecounter*qie_work->nominal_fC();
	charge+=qie_work->nominal_fC();
	LEDTotalChargeSum[dccid-700][spigot][channum-1][samplecounter-1]+=qie_work->nominal_fC();
	LEDChargeCounter[dccid-700][spigot][channum-1][samplecounter-1]++;
//	if (dccid==704&&channum==17&&spigot==11) {
//	  printf("FED %d, chan %.2d, spigot %.2d, TS %d, counter %d nominal2fC %2.2f    ",dccid,channum,spigot,samplecounter,LEDChargeCounter[dccid-700][spigot][channum-1][samplecounter-1],qie_work->nominal_fC());
//	  printf("weightedval %2.2f, total charge DAQ sample %2.2f, total charge on TS %2.2f\n",weightedval,charge,LEDTotalChargeSum[dccid-700][spigot][channum-1][samplecounter-1]);
//	}
	if (samplecounter==10) {
	  m_dbe->setCurrentFolder(baseFolder_+"LEDChargeTiming");//Per event
	  if (charge>100&&charge<6000) {
	    LEDLatencySum[dccid-700][spigot][channum-1]+=weightedval/charge;
	    LEDEventTiming_[dccid-700]->Fill(((NUMCHANNELS+1)*spigot)+channum,weightedval/charge);
	  }
//	  if (charge!=0) {
//	    LEDLatencySum[dccid-700][spigot][channum-1]+=weightedval/charge;
//	    LEDEventTiming_[dccid-700]->Fill(((NUMCHANNELS+1)*spigot)+channum,weightedval/charge);
//	  }
//	  else {
//	    LEDLatencySum[dccid-700][spigot][channum-1]+=0;
//	    LEDEventTiming_[dccid-700]->Fill(((NUMCHANNELS+1)*spigot)+channum,0.0);
//	  }
	  LEDLatencyCounter[dccid-700][spigot][channum-1]++;
	  weightedval=0.0;
	  charge=0.0;
	}
	
	//For every precision data sample in Hcal:
	
      }
      // Prepare for the next round...
      lastcapid=qie_work->capid();
      lastfibchan=qie_work->fiberAndChan();
      qie_work++;
    }
    if (dccid==723 && spigot==3) { //the ZDC spigot
      //const unsigned short* zdcRAW =  htr.getRawData();
      //std::cout  << "ZDC ===> " << zdcRAW[0] << std::endl;
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
} // void HcalLatencyMonitor::unpack(...)

void HcalLatencyMonitor::endOfJob()
{
  char name[128];
  int bimodproblemcounter = 0;
  int rmsproblemcounter   = 0;
  int meanproblemcounter  = 0;
  int meanproblemcounter2  = 0;
  int nomeancounter       = 0;

  TH2F* my2dhist;
  TH1F* my1dhist;

  m_dbe->setCurrentFolder(baseFolder_);//Per run
  sprintf(name,"LEDFibreErrors");
  LEDFibreErrors_ = m_dbe->book2D(name,name,NUMFEDS*2,0,NUMFEDS*2,NUMSPIGOTS*4,0,NUMSPIGOTS*4);

  sprintf(name,"LEDLatency");
  m_dbe->setCurrentFolder(baseFolder_+name);//Per run
  sprintf(name,"LEDChargeRMSDistribution");
  LEDChargeRMS_ = m_dbe->book2D(name,name,(NUMFEDS * (NUMSPIGOTS+1)),0,(NUMFEDS * (NUMSPIGOTS+1)),NUMCHANNELS,0,NUMCHANNELS);

  sprintf(name,"LEDTimingRMSDistribution");
  LEDTimingRMS_ = m_dbe->book1D(name,name,1000,-0.01,9.99);

  sprintf(name,"LEDBimodalDistribution");
  LEDBimodal_ = m_dbe->book1D(name,name,1000,-0.01,9.99);

  sprintf(name,"LEDBimodalBinDifference");
  LEDBimodalBins_ = m_dbe->book1D(name,name,100,-0.5,99.5);

  sprintf(name,"LEDLatencySpectrum");
  LEDSteadyMinoritySpectrum_ = m_dbe->book1D(name,name,100,-2,2);

  sprintf(name,"LED2fCLatency");
  LEDLatency_ = m_dbe->book2D(name,name,(NUMFEDS * (NUMSPIGOTS+1)),0,(NUMFEDS * (NUMSPIGOTS+1)),NUMCHANNELS,0,NUMCHANNELS);
  


  for (int ii = 0; ii < NUMFEDS; ii++) {
    sprintf(name,"LEDMeanCharge");
    m_dbe->setCurrentFolder(baseFolder_+name);//Per run
    sprintf(name,"LEDChargeDistributionFED7%.2d",ii);
    LEDCharge_[ii] = m_dbe->book2D(name,name,(NUMSPIGOTS * (NUMCHANNELS+1)),0,(NUMSPIGOTS * (NUMCHANNELS+1)),NUMTS,0,NUMTS);

    sprintf(name,"FED 7%.2d",ii);
    LEDChargeRMS_->setBinLabel(1+((NUMSPIGOTS+1)*ii),name,1);
    LEDLatency_->setBinLabel(1+((NUMSPIGOTS+1)*ii),name,1);
    LEDFibreErrors_->setBinLabel((ii*2)+1,name,1);

    for (int tt = 0; tt < NUMTS; tt++){
      sprintf(name,"TS %.2d",tt+1);
      LEDCharge_[ii]->setBinLabel(tt+1,name,2);}
    
    for (int jj = 0; jj < NUMSPIGOTS; jj++) {
      sprintf(name,"Spigot %.2d",jj);
      char fibnum[128];
      sprintf(fibnum,"%s\n\nFibre %.2d, %.2d",name,1,5);
      LEDFibreErrors_->setBinLabel((jj*4)+1,fibnum,2);
      sprintf(fibnum,"Fibre %.2d, %.2d",2,6);
      LEDFibreErrors_->setBinLabel((jj*4)+2,fibnum,2);
      sprintf(fibnum,"Fibre %.2d, %.2d",3,7);
      LEDFibreErrors_->setBinLabel((jj*4)+3,fibnum,2);
      sprintf(fibnum,"Fibre %.2d, %.2d",4,8);
      LEDFibreErrors_->setBinLabel((jj*4)+4,fibnum,2);

      LEDEventTiming_[ii]->setBinLabel(1+((NUMCHANNELS+1)*jj),name,1);
      LEDBCNIdle_[ii]->setBinLabel(1+((NUMFIBRES+1)*jj),name,1);
      LEDCharge_[ii]->setBinLabel(1+((NUMCHANNELS+1)*jj),name,1);

      LEDChargeRMS_->setBinLabel(1+((NUMSPIGOTS+1)*ii)+(jj+1),name,1);
      LEDLatency_->setBinLabel(1+((NUMSPIGOTS+1)*ii)+(jj+1),name,1);

      if (activePlots[ii][jj]==1) {
	for (int kk = 0; kk < NUMFIBRES; kk++) {
	  int binnum = 1+((NUMFIBRES+1)*jj)+(kk+1);
	  sprintf(name,"Fibre %.2d",kk+1);
	  LEDBCNIdle_[ii]->setBinLabel(binnum,name,1);}
	for (int kk = 0; kk < NUMCHANNELS; kk++) {
	  if (((kk+1)%3)==2){
	    int binnum = 1+((NUMCHANNELS+1)*jj)+(kk+1);
	    sprintf(name,"Fibre %.2d",(kk/3)+1);
	    LEDEventTiming_[ii]->setBinLabel(binnum,name,1);
	    LEDCharge_[ii]->setBinLabel(binnum,name,1);

	    LEDLatency_->setBinLabel(kk+1,name,2);
	    LEDChargeRMS_->setBinLabel(kk+1,name,2);}
	  
	  if (LEDLatencyCounter[ii][jj]!=0) {
	    double meanvalue = LEDLatencySum[ii][jj][kk]/LEDLatencyCounter[ii][jj][kk];
	    //	    printf("FED 7%.2d, spigot %.2d, channel %.2d, latency %2.2f\n",ii,jj,kk+1,meanvalue);
	    LEDLatency_->Fill(((NUMSPIGOTS+1)*ii)+(jj+1)+0.25,kk+0.25,meanvalue);
	    my2dhist = (TH2F*)LEDEventTiming_[ii]->getTH2F();
	    my1dhist = (TH1F*)my2dhist->ProjectionY("",1+((NUMCHANNELS+1)*jj)+(kk+1),1+((NUMCHANNELS+1)*jj)+(kk+1));
	    double chargerms = my1dhist->GetRMS();
	    if (chargerms!=0) LEDChargeRMS_->Fill(((NUMSPIGOTS+1)*ii)+(jj+1),kk+0.25,chargerms);
	    else              LEDChargeRMS_->Fill(((NUMSPIGOTS+1)*ii)+(jj+1),kk+0.25,-1);
	    LEDTimingRMS_->Fill(chargerms);

	    for (int ll = 0; ll < NUMTS; ll++) {
	      //	      printf("FED 7%.2d, Sp %.2d, Ch %.2d, TS %.2d - total charge %2.2f, counter %d\n",ii,jj,kk+1,ll,LEDTotalChargeSum[ii][jj][kk][ll],LEDChargeCounter[ii][jj][kk][ll]);
	      if (LEDChargeCounter[ii][jj][kk][ll]!=0){
		double meanchargeonts = LEDTotalChargeSum[ii][jj][kk][ll]/LEDChargeCounter[ii][jj][kk][ll];
		LEDCharge_[ii]->Fill(((NUMCHANNELS+1)*jj)+(kk+1),ll,meanchargeonts);}}}}
	
		
	double fibreLEDmean[8]  = {0.0};
	double fibreLEDrms[8]   = {0.0};
	double fibreLEDmode[8]  = {0.0};
	int    fibreLEDmax[8]   = { 0 };
	int    maxBinDiff[8]    = { 0 };
	int    countThresh      = 50;
	double ledRMStol        = 0.1;  //these values need tuning
	double MeanDiffTol      = 0.25;  //
	//	double MeanModeDiffTol1 = 0.1; //
	//	double MeanModeDiffTol2 = 0.12; //
	double RMSDiffTol       = 0.05; //
	for (int pp = 0;  pp < 8; pp++) {
	  int binnum = 1 + ((NUMCHANNELS+1)*jj) + (pp*3) + 1;
	  my2dhist = (TH2F*)LEDEventTiming_[ii]->getTH2F();
	  my1dhist = (TH1F*)my2dhist->ProjectionY("",binnum,binnum+2);
	  fibreLEDmean[pp] = my1dhist->GetMean();
	  fibreLEDrms[pp] = my1dhist->GetRMS();
	  // need to recalculate the mode since the number of bins has changed to 502 from 102, the step size is now 0.02 instead of 0.1
	  fibreLEDmode[pp] = (my1dhist->GetMaximumBin()-1)*0.02-0.02;
	  fibreLEDmax[pp] = my1dhist->GetMaximumBin();
	  int nbins = my1dhist->GetNbinsX();
	  for (int gg = 0; gg < nbins; gg++) {
	    double content = my1dhist->GetBinContent(gg+1);
	    if(content>countThresh) {
	      if (abs((gg+1)-fibreLEDmax[pp])>maxBinDiff[pp]) {
		maxBinDiff[pp] = abs((gg+1)-fibreLEDmax[pp]);}}}
	  LEDBimodal_->Fill(abs(fibreLEDmean[pp]-fibreLEDmode[pp]));
	  LEDBimodalBins_->Fill(abs(maxBinDiff[pp]));
	}
	
	int splitmean[2] = {0};
	for (int qq = 0; qq < 4; qq++) {
	  int fibIndicator[2] = {0};

	  if (fibreLEDrms[qq]>ledRMStol) {
	    if (maxBinDiff[qq]>35) {
	      fibIndicator[0]+=2;
	      if (jj<12) bimodproblemcounter++;}
	    else {
	      fibIndicator[0]+=1;
	      if (jj<12) rmsproblemcounter++;}}
//	  else if (fibreLEDrms[qq]>ledRMStol) {
//	    fibIndicator[0]+=1;
//	    if (jj<12) rmsproblemcounter++;}
	  else if (maxBinDiff[qq]>35) {
	    fibIndicator[0]+=2;
	    if (jj<12) bimodproblemcounter++;}
    
	  
	  if (fibreLEDrms[qq+4]>ledRMStol) {
	    if (maxBinDiff[qq+4]>35) {
	      fibIndicator[1]+=2;
	      if (jj<12) bimodproblemcounter++;}
	    else {
	      fibIndicator[1]+=1;
	      if (jj<12) rmsproblemcounter++;}}
//	  else if (fibreLEDrms[qq+4]>ledRMStol) {
//	    fibIndicator[1]+=1;
//	    if (jj<12) rmsproblemcounter++;}
	  else if (maxBinDiff[qq+4]>35) {
	    fibIndicator[1]+=2;
	    if (jj<12) bimodproblemcounter++;}
	  
	  int rmsErrCounter[2]  = {0};
	  int meanErrCounter[2] = {0};
	  
	  double latentmean[2]  = {0.0};

	  for (int oo = 0; oo < 4; oo++) {
	    if ((abs(fibreLEDmean[qq]-fibreLEDmean[oo]))>MeanDiffTol) meanErrCounter[0]++;
	    if ((abs(fibreLEDrms[qq]-fibreLEDrms[oo]))>RMSDiffTol) rmsErrCounter[0]++;
	    if ((abs(fibreLEDmean[qq+4]-fibreLEDmean[oo+4]))>MeanDiffTol) meanErrCounter[1]++;
	    if ((abs(fibreLEDrms[qq+4]-fibreLEDrms[oo+4]))>RMSDiffTol)  rmsErrCounter[1]++;}
	  
	  if (fibreLEDmean[qq]!=0) {
	    //	      printf("7%.2d, sp%.2d, fi%.2d, mean errors %d\n",ii,jj,qq+1,meanErrCounter[0]);
	    if (meanErrCounter[0]>2) {
	      //the rms and mean are different from the other 3 - split timing
	      if (!(rmsErrCounter[0]>2)) {
		fibIndicator[0]+=4;
		if (jj<12) {
		  meanproblemcounter++;
		  for (int oo = 0; oo < 4; oo++) {
		    if (oo==qq) continue;
		    else latentmean[0]+=fibreLEDmean[oo];}
		  LEDSteadyMinoritySpectrum_->Fill(fibreLEDmean[qq]-(latentmean[0]/3));
		}}}
	    else if (meanErrCounter[0]<2);
	    else {
	      splitmean[0] += 1;}}
	  //	      printf("first group 7%.2d, sp%.2d, fi%.2d\n",ii,jj,qq+1);}}

	  else {//mean is zero, indicating a not present fibre
	    fibIndicator[0]=-7;
	    nomeancounter++;}
	  
	  if (fibreLEDmean[qq+4]!=0) {
	    //	      printf("7%.2d, sp%.2d, fi%.2d, mean errors %d\n",ii,jj,qq+5,meanErrCounter[1]);
	    if (meanErrCounter[1]>2) {
	      if (!(rmsErrCounter[1]>2)) {
		fibIndicator[1]+=4;
		if (jj<12) {
		  meanproblemcounter++;
		  for (int oo = 0; oo < 4; oo++) {
		    if (oo+4==qq+4) continue;
		    else latentmean[1]+=fibreLEDmean[oo+4];}
		  LEDSteadyMinoritySpectrum_->Fill(fibreLEDmean[qq+4]-(latentmean[1]/3));
		}}}
	    else if (meanErrCounter[1]<2);
	    else {
	      splitmean[1] += 1;}}
	  //	      printf("second group 7%.2d, sp%.2d, fi%.2d\n",ii,jj,qq+5);}}

	  else {//mean is zero, indicating a not present fibre
	    fibIndicator[1]=-7;
	    nomeancounter++;}
	
	  LEDFibreErrors_->Fill((ii*2)+0.5,(jj*4)+qq+0.5,fibIndicator[0]);
	  LEDFibreErrors_->Fill((ii*2)+1.5,(jj*4)+qq+0.5,fibIndicator[1]);

	  if (qq==3) {
	    if (splitmean[0]>1) {	
	      if (jj<12) printf("7%.2d    Sp %.2d    first group    %d\n",ii,jj,splitmean[0]);
	      meanproblemcounter2++;}
	    if (splitmean[1]>1) {
	      if (jj<12) printf("7%.2d    Sp %.2d    second group    %d\n",ii,jj,splitmean[1]);
	      meanproblemcounter2++;}}
	  if ((fibIndicator[0]>1)&&jj<12) printf("7%.2d    Sp %.2d    %.2d    %d\n",ii,jj,qq+1,fibIndicator[0]);
	  if ((fibIndicator[1]>1)&&jj<12) printf("7%.2d    Sp %.2d    %.2d    %d\n",ii,jj,qq+5,fibIndicator[1]);}
      }
    }
  }
  printf("===========================Summary================================\n");
  printf("%d fibres indicated a bimodal timing distribution\n",bimodproblemcounter);
  printf("%d fibres indicated a wider than expected timing rms\n",rmsproblemcounter);
  printf("%d fibres indicated an error in the mean of the timing distribution, w/r to the other 3 fibres\n",meanproblemcounter);
  printf("%d groups of 4 fibres indicated a split mean of the timing distribution\n",meanproblemcounter2);
  printf("%d fibres indicated a mean of 0 in the timing distribution\n",nomeancounter);
}
