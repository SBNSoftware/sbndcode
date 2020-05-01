

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "sbndcode/CRT/CRTProducts/CRTData.hh"
#include "art_root_io/TFileService.h"
#include "TH1.h"
#include <bitset>


//  class CRTTrigFilter : public art::EDFilter(fhicl::ParameterSet const& p) {
 class CRTTrigFilter : public art::EDFilter {
   public:
     explicit CRTTrigFilter(fhicl::ParameterSet const& p);
     // FlashTimeFilter(FlashTimeFilter const&) = delete;
     // FlashTimeFilter(FlashTimeFilter&&) = delete;
     // FlashTimeFilter& operator=(FlashTimeFilter const&) = delete;
     // FlashTimeFilter& operator=(FlashTimeFilter&&) = delete;

     //
    virtual ~CRTTrigFilter() { }
    virtual bool filter(art::Event& e) override;
    void    reconfigure(fhicl::ParameterSet const& p);

  private:


   art::ServiceHandle<art::TFileService> tfs;

   std::string fCRTStripModuleLabel;
   std::vector<int>    fmodlistUTopL;       
   std::vector<int>    fmodlistUBotL;       
   std::vector<int>    fmodlistUTopR;       
   std::vector<int>    fmodlistUBotR;       
   std::vector<int>    fmodlistDTopL;       
   std::vector<int>    fmodlistDBotL;       
   std::vector<int>    fmodlistDTopR;       
   std::vector<int>    fmodlistDBotR;       
   bool fstripmatch;
   int fstripshift;
   int fedgecut;
   float    fTimeCoinc;
   float    fADCthresh;


   TH1F *hits;
   TH1F *trigt;
   TH1F *trig;
                    
  };

  CRTTrigFilter::CRTTrigFilter(fhicl::ParameterSet const& p): EDFilter{p} 
  {


    // Create histograms
    hits = tfs->make<TH1F>("hits","hits",14,39.5,63.5);
    //    hits->GetXaxis()->SetTitle("Track Length (cm)");
    trigt = tfs->make<TH1F>("trigt","trigt",200,-3000,3000);  // trig time in us
    trig = tfs->make<TH1F>("trig","trig",2,-0.5,1.5);  

    this->reconfigure(p);
  }

  void CRTTrigFilter::reconfigure(fhicl::ParameterSet const& p)
  {

    fCRTStripModuleLabel = p.get< std::string >("CRTStripModuleLabel","crt");
    fmodlistUTopL = p.get<std::vector<int>>("ModuleListUpstreamTopLeft");
    fmodlistUBotL = p.get<std::vector<int>>("ModuleListUpstreamBotLeft");
    fmodlistUTopR = p.get<std::vector<int>>("ModuleListUpstreamTopRight");
    fmodlistUBotR = p.get<std::vector<int>>("ModuleListUpstreamBotRight");
    fmodlistDTopL = p.get<std::vector<int>>("ModuleListDownstreamTopLeft");
    fmodlistDBotL = p.get<std::vector<int>>("ModuleListDownstreamBotLeft");
    fmodlistDTopR = p.get<std::vector<int>>("ModuleListDownstreamTopRight");
    fmodlistDBotR = p.get<std::vector<int>>("ModuleListDownstreamBotRight");
    fstripmatch = p.get<bool>("RequireStripMatch",true);
    fstripshift = p.get<int>("AllowedStripShift",4);
    fedgecut = p.get<int>("EdgeCutInStrips",12);
    if (fedgecut>31) fedgecut=31;
    fTimeCoinc = p.get<float>("StripTimeCoincidence",0.2);
    fADCthresh = p.get<float>("ADCthresh",500.0);
  }

  bool CRTTrigFilter::filter(art::Event& e)
  {

    bool KeepMe = false;
    int event = e.id().event();
    if (event%1000==0) std::cout << "event " << event << std::endl;
    //    unsigned long long planeUtop, planeUbot, planeDtop, planeDbot;
    uint64_t planeUtopL, planeUbotL, planeDtopL, planeDbotL;
    uint64_t planeUtopR, planeUbotR, planeDtopR, planeDbotR;
    uint64_t planeUpStL, planeDownStL;
    uint64_t planeUpStR, planeDownStR;

    //    uint64 planeUpSt,planeDownSt;

    int nstr=0;
    art::Handle<std::vector<sbnd::crt::CRTData> > crtStripListHandle;
    std::vector<art::Ptr<sbnd::crt::CRTData> > striplist;
    if (e.getByLabel(fCRTStripModuleLabel, crtStripListHandle))  {
      //    if (e.getByLabel("crt", crtStripListHandle))  {
      art::fill_ptr_vector(striplist, crtStripListHandle);
      nstr = striplist.size();
    }
    //    std::cout << "number of crt strips " << nstr << std::endl;
    
    bool trigKeep = false;   
    uint64_t edgecutL,edgecutR;
    edgecutR = ~((uint64_t)(pow(2,fedgecut)-1) | ((uint64_t)(pow(2,32)-1)<< 32));
    edgecutL = ~(((uint64_t)(pow(2,fedgecut)-1) << (32-fedgecut))| ((uint64_t)(pow(2,32)-1) <<32));
    //    std::cout << "left " << std::bitset<64>(edgecutL) << " right " << std::bitset<64>(edgecutR) << std::endl;
    //  left  0000000000000000000000000000000000000000000011111111111111111111 
    //  right 0000000000000000000000000000000011111111111111111111000000000000
    
    float stripwidth = 11.2; //cm
    float driftvel = 0.16;  //cm/us

    // set cut limits on drift window time expectation of track from CRT hit.
    float rwcut = fedgecut*stripwidth/driftvel;
    float rwcutlow = -200.0 - rwcut;
    float rwcuthigh = 1500.00 + rwcut;

    for (int i = 0; i<nstr-2; i+=2){
      if ((striplist[i]->ADC()+striplist[i+1]->ADC())>fADCthresh) { 
	uint32_t chan1 = striplist[i]->Channel();
	int strip1 = (chan1 >> 1) & 15;
	int module1 = (chan1>> 5);
	uint32_t ttime1 = striplist[i]->T0();
	float ctime1 = ttime1/16.;
	if (ttime1 > 2147483648) {
	  ctime1 = ((ttime1-4294967296)/16.);
	}

	//
	for (int j = i+2; j<nstr; j+=2){
	  if ((striplist[j]->ADC()+striplist[j+1]->ADC())>fADCthresh) { 	
	  bool match = false;
	  uint32_t chan2 = striplist[j]->Channel();
	  int strip2 = (chan2 >> 1) & 15;
	  int module2 = (chan2>> 5);
	  //
	  uint32_t ttime2 = striplist[j]->T0();
	  //  T0 in units of ticks, but clock frequency (16 ticks = 1 us) is wrong.
	  // ints were stored as uints, need to patch this up
	  float ctime2 = ttime2/16.;
	  if (ttime2 > 2147483648) {
	    ctime2 = ((ttime2-4294967296)/16.);
	  }
	  //
	  float diff = abs(ctime1-ctime2);
	  //	  if (diff<=fTimeCoinc && ctime1>-1300 && ctime1<1300) {
	  if (diff<=fTimeCoinc ) {
	    
	    
	    planeUtopL=0; planeUbotL=0; planeDtopL=0; planeDbotL=0;
	    planeUtopR=0; planeUbotR=0; planeDtopR=0; planeDbotR=0;
	    for (size_t im=0;im<fmodlistUTopL.size();++im) {
	      size_t iswap = fmodlistUTopL.size()-im-1;
	      if (module1==fmodlistUTopL[im]) planeUtopL+=(1 << ((15-strip1)+iswap*16));
	      if (module2==fmodlistUTopL[im]) planeUtopL+=(1 << ((15-strip2)+iswap*16));
	    }	
	    for (size_t im=0;im<fmodlistUBotL.size();++im) {
	      size_t iswap = fmodlistUBotL.size()-im-1;
	      if (module1==fmodlistUBotL[im]) planeUbotL+=(1 << ((15-strip1)+iswap*16));
	      if (module2==fmodlistUBotL[im]) planeUbotL+=(1 << ((15-strip2)+iswap*16));
	    }	
	    for (size_t im=0;im<fmodlistUTopR.size();++im) {
	      size_t iswap = fmodlistUTopR.size()-im-1;
	      if (module1==fmodlistUTopR[im]) planeUtopR+=(1 << ((15-strip1)+iswap*16));
	      if (module2==fmodlistUTopR[im]) planeUtopR+=(1 << ((15-strip2)+iswap*16));
	    }	
	    for (size_t im=0;im<fmodlistUBotR.size();++im) {
	      size_t iswap = fmodlistUBotR.size()-im-1;
	      if (module1==fmodlistUBotR[im]) planeUbotR+=(1 << ((15-strip1)+iswap*16));
	      if (module2==fmodlistUBotR[im]) planeUbotR+=(1 << ((15-strip2)+iswap*16));
	    }	
	    for (size_t im=0;im<fmodlistDTopL.size();++im) {
	      size_t iswap = fmodlistDTopL.size()-im-1;
	      if (module1==fmodlistDTopL[im]) planeDtopL+=(1 << ((15-strip1)+iswap*16));
	      if (module2==fmodlistDTopL[im]) planeDtopL+=(1 << ((15-strip2)+iswap*16));
	    }	
	    for (size_t im=0;im<fmodlistDBotL.size();++im) {
	      size_t iswap = fmodlistDBotL.size()-im-1;
	      if (module1==fmodlistDBotL[im]) planeDbotL+=(1 << ((15-strip1)+iswap*16));
	      if (module2==fmodlistDBotL[im]) planeDbotL+=(1 << ((15-strip2)+iswap*16));
	    }	
	    for (size_t im=0;im<fmodlistDTopR.size();++im) {
	      size_t iswap = fmodlistDTopR.size()-im-1;
	      if (module1==fmodlistDTopR[im]) planeDtopR+=(1 << ((15-strip1)+iswap*16));
	      if (module2==fmodlistDTopR[im]) planeDtopR+=(1 << ((15-strip2)+iswap*16));
	    }	
	    for (size_t im=0;im<fmodlistDBotR.size();++im) {
	      size_t iswap = fmodlistDBotR.size()-im-1;
	      if (module1==fmodlistDBotR[im]) planeDbotR+=(1 << ((15-strip1)+iswap*16));
	      if (module2==fmodlistDBotR[im]) planeDbotR+=(1 << ((15-strip2)+iswap*16));
	    }	

	    planeUpStL = planeUtopL | planeUbotL;
	    planeUpStR = planeUtopR | planeUbotR;
	    planeDownStL = planeDtopL | planeDbotL;
	    planeDownStR = planeDtopR | planeDbotR;
	    if ((planeUpStR & edgecutR) && (planeDownStR & edgecutR)) {

	      if ( planeUpStR & planeDownStR ) match=true;
	      for (int is=1;is<=fstripshift;++is) {
		uint64_t temp; temp=(planeUpStR << is);
		if (temp & planeDownStR) match=true;
		temp=(planeDownStR << is);
		if (temp & planeUpStR)  match=true;
	      }
	    }
	    if ((planeUpStL & edgecutL) && (planeDownStL & edgecutL)) {
	      if ( planeUpStL & planeDownStL ) match=true;
	      for (int is=1;is<=fstripshift;++is) {
		uint64_t temp; temp=(planeUpStL << is);
		if (temp & planeDownStL) match=true;
		temp=(planeDownStL << is);
		if (temp & planeUpStL)  match=true;
	      }
	    }
	    if (match) {
	      // std::cout << " before drift window cut " << std::endl;
	      // std::cout << " module and strip 1 " << module1 << " " << strip1 << std::endl;
	      // std::cout << " module and strip 2 " << module2 << " " << strip2 << std::endl;
	      // std::cout << "  upstream " << std::bitset<64>(planeUpStR) <<
	      // 	std::endl << "downstream " << std::bitset<64>(planeDownStR) << 
	      // 	std::endl;
	      // std::cout << " ----------------------------------" << std::endl;
	      float xpos;  // in cm
	      if (module1<module2) xpos=((int(module1/2)-20)*16+strip1+0.5)*stripwidth-358.4;
	      else xpos=((int(module2/2)-20)*16+strip2+0.5)*stripwidth-358.4;
	      float dtime;  // in us
	      if (xpos>0) dtime = ctime1+((200.0-xpos)/driftvel);
	      else dtime = ctime1+((200.0+xpos)/driftvel);
	      if (dtime<rwcutlow || dtime>rwcuthigh ) match=false;
	      // std::cout << "Event " << event << " " << module1 << " " << xpos << " " <<  ctime1 << " " << dtime << std::endl;
	    }
	  }	  
	  if (match) {
	   // std::cout << "Found One! Event " << event << 
	   //     "   m/s1 m/s2 " << module1 << " " << strip1 << " " << 
	   //     module2 << " " << strip2 << " " << ctime1 << " " <<ctime2 << std::endl;
	    KeepMe=true;
	  }
	  }} // end loop over second strip
	}} // end loop over first strip
    if (trigKeep) trig->Fill(1.0); else trig->Fill(0.0);
    return KeepMe;

  }

  // A macro required for a JobControl module.
  DEFINE_ART_MODULE(CRTTrigFilter)


