#ifndef SPACECHARGE_SPACECHARGESBND_H
#define SPACECHARGE_SPACECHARGESBND_H

// LArSoft libraries
#include "larevt/SpaceCharge/SpaceCharge.h"

// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"

// Others
#include <string>
#include <vector>
#include <TGraph.h>
#include <TF1.h>
#include <TH3.h>
#include <TFile.h>

namespace spacecharge
{
    class SpaceChargeSBND : public SpaceCharge
    {

    public:
	explicit SpaceChargeSBND(fhicl::ParameterSet const& pset);
	SpaceChargeSBND(SpaceChargeSBND const&) = delete;
	virtual ~SpaceChargeSBND() = default;

	bool Configure(fhicl::ParameterSet const& pset);
	bool Update(uint64_t ts = 0);

	bool EnableSimSpatialSCE() const override;
	bool EnableSimEfieldSCE() const override;
	bool EnableCalSpatialSCE() const override;
    bool EnableCalEfieldSCE() const override;
	
	 bool EnableCorrSCE() const override {return (EnableCalSpatialSCE()||EnableCalEfieldSCE()) ;}
	 
	geo::Vector_t GetPosOffsets(geo::Point_t const& point) const override;
	geo::Vector_t GetEfieldOffsets(geo::Point_t const& point) const override;
	geo::Vector_t GetCalPosOffsets(geo::Point_t const& point, int const& TPCid = 1) const override;
	geo::Vector_t GetCalEfieldOffsets(geo::Point_t const& point, int const& TPCid = 1) const override { return {0.,0.,0.}; }

    private:
    protected:

	int initialSpatialFitPolN[3] = {3, 4, 3};
	int intermediateSpatialFitPolN[3] = {4, 4, 4};
	int initialEFieldFitPolN[3] = {3, 3, 3};
	int intermediateEFieldFitPolN[3] = {6, 4, 4};
	double DriftField = 500.0; // 500 V/cm

	bool fEnableSimSpatialSCE;
	bool fEnableSimEfieldSCE;
	bool fEnableCalSpatialSCE;
	bool fEnableCalEfieldSCE;
	bool fEnableCorrSCE;
	bool f_2D_drift_sim_hack;

	std::string fRepresentationType;
	std::string fInputFilename;

	std::vector<double> GetPosOffsetsParametric(double xVal, double yVal, double zVal) const;
	double GetOnePosOffsetParametric(double xVal, double yVal, double zVal, std::string axis) const;
	std::vector<double> GetEfieldOffsetsParametric(double xVal, double yVal, double zVal) const;
	double GetOneEfieldOffsetParametric(double xVal, double yVal, double zVal, std::string axis) const;
	double TransformX(double xVal) const;
	double TransformY(double yVal) const;
	double TransformZ(double zVal) const;
	bool IsInsideBoundaries(double xVal, double yVal, double zVal) const;

	//to store Voxelized_TH3 histograms
	std::vector<TH3F*> SCEhistograms = std::vector<TH3F*>(9);

	TGraph *gSpatialGraphX[99][99];
	TF1 *intermediateSpatialFitFunctionX[99];
	TF1 *initialSpatialFitFunctionX;
	TGraph *gSpatialGraphY[99][99];
	TF1 *intermediateSpatialFitFunctionY[99];
	TF1 *initialSpatialFitFunctionY;
	TGraph *gSpatialGraphZ[99][99];
	TF1 *intermediateSpatialFitFunctionZ[99];
	TF1 *initialSpatialFitFunctionZ;

	TGraph *gEFieldGraphX[99][99];
	TF1 *intermediateEFieldFitFunctionX[99];
	TF1 *initialEFieldFitFunctionX;
	TGraph *gEFieldGraphY[99][99];
	TF1 *intermediateEFieldFitFunctionY[99];
	TF1 *initialEFieldFitFunctionY;
	TGraph *gEFieldGraphZ[99][99];
	TF1 *intermediateEFieldFitFunctionZ[99];
	TF1 *initialEFieldFitFunctionZ;
}; // class SpaceChargeSBND
} //namespace spacecharge
#endif // SPACECHARGE_SPACECHARGESBND_H
