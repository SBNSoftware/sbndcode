
#include "DigitalNoiseChannelStatus.h"
#include "TMath.h"
#include "lardataobj/RawData/raw.h"
#include "art/Framework/Principal/Event.h"

sbnd::DigitalNoiseChannelStatus::DigitalNoiseChannelStatus(fhicl::ParameterSet const& p, art::ActivityRegistry& areg)
// :
// Initialize member data here.
{
  fRMSCutWire = p.get<float>("RMSCutWire",100.0);
  fRMSCutRawDigit = p.get<float>("RMSCutRawDigit",100.0);
  fNBADCutRawDigit = p.get<float>("NBADCutRawDigit",5);
  fRawDigitLabel = p.get<std::string>("RawDigitLabel","daq");
  fRecobWireLabel = p.get<std::string>("RecobWireLabel","caldata");
  fNAwayFromPedestalRawDigit = p.get<int>("NAwayFromPedestalRawDigit",100);
  fDistFromPedestalRawDigit = p.get<int>("DistFromPedestalRawDigit",100);
  fNAwayFromPedestalRecobWire = p.get<int>("NAwayFromPedestalRecobWire",100);
  fDistFromPedestalRecobWire = p.get<float>("DistFromPedestalRecobWire",100);
  
  areg.sPreProcessEvent.watch(this, &sbnd::DigitalNoiseChannelStatus::pub_PrepEvent);
}

bool sbnd::DigitalNoiseChannelStatus::IsBad(raw::ChannelID_t chan) const
{
  return (fDNChannels.find(chan) != fDNChannels.end());
}

const std::unordered_set<raw::ChannelID_t> & sbnd::DigitalNoiseChannelStatus::GetSetOfBadChannels() const
{
  return fDNChannels;
}

size_t sbnd::DigitalNoiseChannelStatus::NBadChannels() const
{
  return fDNChannels.size();
}

// prepare channel status data for the event

// put cuts on RMS, number of 0xBAD samples, and whether the differences between a sample and
// the first are all even (catches all the powers of two).

void sbnd::DigitalNoiseChannelStatus::pub_PrepEvent(const art::Event& evt, art::ScheduleContext)
{
  fDNChannels.clear();
  if (fRawDigitLabel != "")
    {
      auto const& rawdigits = evt.getProduct<std::vector<raw::RawDigit>>(fRawDigitLabel);
      for (const auto& rd : rawdigits)
	{
	  bool chanbad = false;
	  chanbad |= (rd.GetSigma() > fRMSCutRawDigit);
	  int nhexbad = 0;
	  std::vector<short> rawadc;
	  rawadc.resize(rd.Samples());
	  raw::Uncompress(rd.ADCs(), rawadc, rd.GetPedestal(), rd.Compression());
	  bool alleven = true;
	  int naway = 0;
	  const short adc0 = rawadc.at(0);
	  for (size_t i=0; i< rd.Samples(); ++i)
	    {
	      const short adc = rawadc.at(i);
	      if (adc == 0xBAD) ++nhexbad;
	      alleven &= ( ((adc - adc0) % 2) == 0 );
	      if (std::abs(adc - rd.GetPedestal()) > fDistFromPedestalRawDigit)
		{
		  ++naway;
		}
	    }
	  chanbad |= (nhexbad > fNBADCutRawDigit);
	  chanbad |= alleven;
	  chanbad |= (naway > fNAwayFromPedestalRawDigit);
	  if (chanbad) fDNChannels.emplace(rd.Channel());
	}
    }
  else if (fRecobWireLabel != "")
    {
      art::InputTag rwtag(fRecobWireLabel);
      auto const& rbwires = evt.getProduct<std::vector<recob::Wire>>(fRecobWireLabel);
      for (const auto& rw : rbwires)
	{
	  bool chanbad = false;
	  auto rms = TMath::RMS(rw.NSignal(),rw.Signal().data());
	  chanbad |= (rms > fRMSCutWire);
	  auto mns = TMath::Mean(rw.NSignal(),rw.Signal().data());
	  int naway = 0;
	  for (size_t i=0; i<rw.NSignal(); ++i)
	    {
	      if (std::abs(rw.Signal()[i] - mns) > fDistFromPedestalRecobWire)
		{
		  naway ++;
		}
	    }
	  chanbad |= (naway > fNAwayFromPedestalRecobWire);
	  if (chanbad) fDNChannels.emplace(rw.Channel());
	}
    }
  //std::cout << "sbnd::DigitalNoiseStatus_service: NBadChannels: " << fDNChannels.size() << std::endl;
  //for (const auto& c : fDNChannels)
  //  {
  //    std::cout << "noisy chan: " << c << std::endl;
  //  }
}

DEFINE_ART_SERVICE(sbnd::DigitalNoiseChannelStatus)


