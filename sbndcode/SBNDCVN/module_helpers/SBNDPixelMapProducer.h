#ifndef LCVN_SBNDPIXELMAPPRODUCER_H
#define LCVN_SBNDPIXELMAPPRODUCER_H

#include  <iostream>
#include  <ostream>
#include  <list>
#include  <algorithm>
#include <numeric>
#include <cstdlib>

#include "canvas/Persistency/Common/Ptr.h"
#include "larrecodnn/CVN/interfaces/PixelMapProducer.h"
#include "sbndcode/SBNDCVN/module_helpers/SBNDPixelMap.h"

namespace lcvn
{
  template <class T, class U> class SBNDPixelMapProducer : public PixelMapProducer<T,U>
  {
    public:
	SBNDPixelMapProducer(unsigned int nWire, unsigned int nTdc, double tRes, double threshold = 0.):PixelMapProducer<T,U>::PixelMapProducer(nWire, nTdc, tRes, threshold){std::cout << "============ Calling the function SBNDPixelMapProducer::SBNDPixelMapProducer() ==============\n";}
        SBNDPixelMapProducer():PixelMapProducer<T,U>::PixelMapProducer(){std::cout << "============ Calling the function SBNDPixelMapProducer::SBNDPixelMapProducer() ==============\n";}
	SBNDPixelMapProducer(const fhicl::ParameterSet& pset):PixelMapProducer<T,U>::PixelMapProducer(pset),fverbose(pset.get<bool>("verbose")),fChangeWireNo(pset.get<bool>("ChangeWireNo")),fReadoutSize(pset.get<double>("ReadoutSize")),fShiftT(pset.get<float>("ShiftT")),fInductionWires(pset.get<int>("InductionWires")),fFlipInductionView(pset.get<bool>("FlipInductionView")),fUseT(pset.get<bool>("UseT")) {std::cout << "============ Calling the function SBNDPixelMapProducer::SBNDPixelMapProducer() ==============\n";}
	Boundary DefineBoundary(detinfo::DetectorPropertiesData const& detProp,const std::vector< const T* >& cluster) override;
	void ConvertLocaltoGlobal(geo::WireID wireid, unsigned int &globalWire, unsigned int &globalPlane) const override; 
	void ConvertLocaltoGlobalTDC(geo::WireID wireid, double localTDC, unsigned int &globalWire, unsigned int &globalPlane, double &globalTDC) const override;
	PixelMap CreateMapGivenBoundary(detinfo::DetectorPropertiesData const& detProp,const std::vector< const T* >& cluster,const Boundary& bound) override;
	PixelMap CreateMap(detinfo::DetectorPropertiesData const& detProp,const std::vector<art::Ptr<T>>& cluster) override;
	PixelMap CreateMap(detinfo::DetectorPropertiesData const& detProp,const std::vector< const T* >& cluster) override;
        SBNDPixelMap SBNDCreateMapGivenBoundary(detinfo::DetectorPropertiesData const& detProp,const std::vector< const T* >& cluster,const Boundary& bound);
	SBNDPixelMap SBNDCreateMap(detinfo::DetectorPropertiesData const& detProp,const std::vector<art::Ptr<T>>& cluster);
	SBNDPixelMap SBNDCreateMap(detinfo::DetectorPropertiesData const& detProp,const std::vector< const T* >& cluster);
	void Set_fT0_value(float value);
    protected:
        bool fverbose;
        bool fChangeWireNo; 
        double fReadoutSize; // in time ticks
	float fShiftT; // size of the back/front porch
	int fInductionWires; // number of wires in the first induction plane
	bool fFlipInductionView; // should we flip the induction view
	bool fUseT; // 
	float fT0; // T0 coming from the PFP particles (in time ticks)
  };
  
  typedef SBNDPixelMapProducer<recob::Hit, lcvn::HitHelper> SBNDPixelMapHitProducer;
  typedef SBNDPixelMapProducer<recob::Wire, lcvn::WireHelper> SBNDPixelMapWireProducer;
  typedef SBNDPixelMapProducer<sim::SimChannel, lcvn::SimChannelHelper> SBNDPixelMapSimProducer;
}

#endif // LCVN_SBNDPIXELMAPPRODUCER_H
