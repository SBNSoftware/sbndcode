///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SpaceChargeSBND.cxx; brief implementation of class for storing/accessing space charge distortions for SBND
// arbint@bnl.gov
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdio.h>

// LArSoft includes
#include "sbndcode/SpaceCharge/SpaceChargeSBND.h"

// Framework includes
#include "cetlib_except/exception.h"

spacecharge::SpaceChargeSBND::SpaceChargeSBND(fhicl::ParameterSet const& pset)
{
    Configure(pset);
}

bool spacecharge::SpaceChargeSBND::Configure(fhicl::ParameterSet const& pset)
{
    fEnableSimSpatialSCE = pset.get<bool>("EnableSimSpatialSCE");
    fEnableSimEfieldSCE = pset.get<bool>("EnableSimEfieldSCE");
    //fEnableCorrSCE = pset.get<bool>("EnableCorrSCE");
    fEnableCalSpatialSCE = pset.get<bool>("EnableCalSpatialSCE");
    fEnableCalEfieldSCE = pset.get<bool>("EnableCalEfieldSCE");
    f_2D_drift_sim_hack = pset.get<bool>("is2DdriftSimHack","false");

    if((fEnableSimSpatialSCE == true) || (fEnableSimEfieldSCE == true))
        {
            fRepresentationType = pset.get<std::string>("RepresentationType");
            fInputFilename = pset.get<std::string>("InputFilename");

            std::string fname;
            cet::search_path sp("FW_SEARCH_PATH");
            sp.find_file(fInputFilename, fname);

            std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
            if(!infile->IsOpen())
                {
                    throw cet::exception("SpaceChargeSBND") << "Could not find the space charge effect file '" << fname << "'!\n";
                }

            if(fRepresentationType == "Voxelized_TH3"){
      	      std::cout << "begin loading voxelized TH3s..." << std::endl;

      	      //Load in histograms
      	      TH3F* hTrueFwdX = (TH3F*) infile->Get("TrueFwd_Displacement_X");
      	      TH3F* hTrueFwdY = (TH3F*) infile->Get("TrueFwd_Displacement_Y");
      	      TH3F* hTrueFwdZ = (TH3F*) infile->Get("TrueFwd_Displacement_Z");
      	      TH3F* hTrueBkwdX = (TH3F*) infile->Get("TrueBkwd_Displacement_X");
      	      TH3F* hTrueBkwdY = (TH3F*) infile->Get("TrueBkwd_Displacement_Y");
      	      TH3F* hTrueBkwdZ = (TH3F*) infile->Get("TrueBkwd_Displacement_Z");
      	      TH3F* hTrueEFieldX = (TH3F*) infile->Get("True_ElecField_X");
      	      TH3F* hTrueEFieldY = (TH3F*) infile->Get("True_ElecField_Y");
      	      TH3F* hTrueEFieldZ = (TH3F*) infile->Get("True_ElecField_Z");

      	      //https://root.cern.ch/doc/master/classTH1.html#a0367fe04ae8709fd4b82795d0a5462c3
      	      //Set hist directories so they can be referenced elsewhere
      	      //This needs to be done because they were read in from ext file
      	      //Note this is not a property of the TH3F, so does't survive copying
      	      hTrueFwdX->SetDirectory(0);
      	      hTrueFwdY->SetDirectory(0);
      	      hTrueFwdZ->SetDirectory(0);
      	      hTrueBkwdX->SetDirectory(0);
      	      hTrueBkwdY->SetDirectory(0);
      	      hTrueBkwdZ->SetDirectory(0);
      	      hTrueEFieldX->SetDirectory(0);
      	      hTrueEFieldY->SetDirectory(0);
      	      hTrueEFieldZ->SetDirectory(0);

      	      //SCEhistograms can be accessed globally in this script
      	      SCEhistograms = {hTrueFwdX, hTrueFwdY, hTrueFwdZ,
      			       hTrueBkwdX, hTrueBkwdY, hTrueBkwdZ,
      			       hTrueEFieldX, hTrueEFieldY, hTrueEFieldZ};


      	      std::cout << "...finished loading TH3s" << std::endl;
      	    }else if(fRepresentationType == "Parametric")
                {
                    for(int i = 0; i < initialSpatialFitPolN[0] + 1; i++)
                        {
                            for(int j = 0; j < intermediateSpatialFitPolN[0] + 1; j++)
                                {
                                    gSpatialGraphX[i][j] = (TGraph*)infile->Get(Form("deltaX/g%i_%i", i, j));

                                }
                            intermediateSpatialFitFunctionX[i] = new TF1(Form("intermediateSpatialFitFunctionX_%i", i), Form("pol%i", intermediateSpatialFitPolN[0]));
                        }
                    for(int i = 0; i < initialSpatialFitPolN[1] + 1; i++)
                        {
                            for(int j = 0; j < intermediateSpatialFitPolN[1] + 1; j++)
                                {
                                    gSpatialGraphY[i][j] = (TGraph*)infile->Get(Form("deltaY/g%i_%i", i, j));

                                }
                            intermediateSpatialFitFunctionY[i] = new TF1(Form("intermediateSpatialFitFunctionY_%i", i), Form("pol%i", intermediateSpatialFitPolN[1]));
                        }
                    for(int i = 0; i < initialSpatialFitPolN[2] + 1; i++)
                        {
                            for(int j = 0; j < intermediateSpatialFitPolN[2] + 1; j++)
                                {
                                    gSpatialGraphZ[i][j] = (TGraph*)infile->Get(Form("deltaZ/g%i_%i", i, j));

                                }
                            intermediateSpatialFitFunctionZ[i] = new TF1(Form("intermediateSpatialFitFunctionZ_%i", i), Form("pol%i", intermediateSpatialFitPolN[2]));
                        }

                    initialSpatialFitFunctionX =  new TF1("initialSpatialFitFunctionX", Form("pol%i", initialSpatialFitPolN[0]));
                    initialSpatialFitFunctionY =  new TF1("initialSpatialFitFunctionY", Form("pol%i", initialSpatialFitPolN[1]));
                    initialSpatialFitFunctionZ =  new TF1("initialSpatialFitFunctionZ", Form("pol%i", initialSpatialFitPolN[2]));

                    for(int i = 0; i < initialEFieldFitPolN[0] + 1; i++)
                        {
                            for(int j = 0; j < intermediateEFieldFitPolN[0] + 1; j++)
                                {
                                    gEFieldGraphX[i][j] = (TGraph*)infile->Get(Form("deltaEx/g%i_%i", i, j));

                                }
                            intermediateEFieldFitFunctionX[i] = new TF1(Form("intermediateEFieldFitFunctionX_%i", i), Form("pol%i", intermediateEFieldFitPolN[0]));
                        }
                    for(int i = 0; i < initialEFieldFitPolN[1] + 1; i++)
                        {
                            for(int j = 0; j < intermediateEFieldFitPolN[1] + 1; j++)
                                {
                                    gEFieldGraphY[i][j] = (TGraph*)infile->Get(Form("deltaEy/g%i_%i", i, j));

                                }
                            intermediateEFieldFitFunctionY[i] = new TF1(Form("intermediateEFieldFitFunctionY_%i", i), Form("pol%i", intermediateEFieldFitPolN[1]));
                        }
                    for(int i = 0; i < initialEFieldFitPolN[2] + 1; i++)
                        {
                            for(int j = 0; j < intermediateEFieldFitPolN[2] + 1; j++)
                                {
                                    gEFieldGraphZ[i][j] = (TGraph*)infile->Get(Form("deltaEz/g%i_%i", i, j));

                                }
                            intermediateEFieldFitFunctionZ[i] = new TF1(Form("intermediateEFieldFitFunctionZ_%i", i), Form("pol%i", intermediateEFieldFitPolN[2]));
                        }

                    initialEFieldFitFunctionX =  new TF1("initialEFieldFitFunctionX", Form("pol%i", initialEFieldFitPolN[0]));
                    initialEFieldFitFunctionY =  new TF1("initialEFieldFitFunctionY", Form("pol%i", initialEFieldFitPolN[1]));
                    initialEFieldFitFunctionZ =  new TF1("initialEFieldFitFunctionZ", Form("pol%i", initialEFieldFitPolN[2]));
                }else{
                  std::cout << "fRepresentationType not known!!!" << std::endl;
                }
            infile->Close();
        }

    if(fEnableCorrSCE == true)
        {
            // Grab other parameters from pset
        }
    return true;
}

bool spacecharge::SpaceChargeSBND::Update(uint64_t ts)
{
    if (ts == 0)
        {
            return false;
        }

    return true;
}

// Whether or not to turn simulation of SCE on for spatial distortions
bool spacecharge::SpaceChargeSBND::EnableSimSpatialSCE() const
{
    return fEnableSimSpatialSCE;
}

// Whether or not to turn simulation of SCE on for E-field distortions
bool spacecharge::SpaceChargeSBND::EnableSimEfieldSCE() const
{
    return fEnableSimEfieldSCE;
}

// Whether or not to apply SCE corrections
//bool spacecharge::SpaceChargeSBND::EnableCorrSCE() const
//{
//    return fEnableCorrSCE;
//}

// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceChargeSBND::EnableCalSpatialSCE() const
{
  return fEnableCalSpatialSCE;
}

// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceChargeSBND::EnableCalEfieldSCE() const
{
  return fEnableCalEfieldSCE;
}

// Primary working method of service that provides position offsets
geo::Vector_t spacecharge::SpaceChargeSBND::GetPosOffsets(geo::Point_t const& point) const
{
    std::vector<double> thePosOffsets;
    double xx=point.X(), yy=point.Y(), zz=point.Z();

    if(fRepresentationType == "Voxelized_TH3"){
      //handle OOAV by projecting edge cases
      if(xx<-199.999){xx=-199.999;}
      else if(xx>199.999){xx=199.999;}
      if(yy<-199.999){yy=-199.999;}
      else if(yy>199.999){yy=199.999;}
      if(zz<0.001){zz=0.001;}
      else if(zz>499.999){zz=499.999;}
      //larsim requires negative sign in TPC 0
      int corr = 1;

      // ========================================================
      // This hack is to account for a known issue with the 
      // space charge implementation for the 2D simulation.
      // See https://cdcvs.fnal.gov/redmine/issues/28099
      // This should be removed once the appropriate upgrades
      // have been implemented.
      // ========================================================
      if(f_2D_drift_sim_hack == true)
	corr = -1; 

      if (xx < 0) {
	corr = -1; 
      }
            
      double offset_x=0., offset_y=0., offset_z=0.;
      offset_x = corr*SCEhistograms.at(0)->Interpolate(xx,yy,zz);
      offset_y = SCEhistograms.at(1)->Interpolate(xx,yy,zz);
      offset_z = SCEhistograms.at(2)->Interpolate(xx,yy,zz);
      thePosOffsets = {offset_x, offset_y, offset_z};

    }else if(fRepresentationType == "Parametric"){
      if(IsInsideBoundaries(point.X(), point.Y(), point.Z()) == false){
        thePosOffsets.resize(3, 0.0);
      }else{
        // GetPosOffsetsParametric returns m; the PosOffsets should be in cm
        thePosOffsets = GetPosOffsetsParametric(xx, yy, zz);
        for(int i=0; i<3; i++){thePosOffsets[i]=100.*thePosOffsets[i];}
      }
    }
    else{
      thePosOffsets.resize(3, 0.0);
    }

    return { thePosOffsets[0], thePosOffsets[1], thePosOffsets[2] };
}

// Provides backward position offset for analyzers (TH3)
geo::Vector_t spacecharge::SpaceChargeSBND::GetCalPosOffsets(geo::Point_t const& point, int const& TPCid ) const
{
  std::vector<double> theCalPosOffsets;
  double xx=point.X(), yy=point.Y(), zz=point.Z();

  if(fRepresentationType == "Voxelized_TH3"){
    //handle OOAV by projecting edge cases
    if(xx<-199.999){xx=-199.999;}
    else if(xx>199.999){xx=199.999;}
    if(yy<-199.999){yy=-199.999;}
    else if(yy>199.999){yy=199.999;}
    if(zz<0.001){zz=0.001;}
    else if(zz>499.999){zz=499.999;}
    //correct for charge drifted across cathode
    if ((TPCid == 0) and (xx > -2.5)) { xx = -2.5; }
    if ((TPCid == 1) and (xx < 2.5)) { xx = 2.5; }
    double offset_x=0., offset_y=0., offset_z=0.;
    offset_x = SCEhistograms.at(3)->Interpolate(xx,yy,zz);
    offset_y = SCEhistograms.at(4)->Interpolate(xx,yy,zz);
    offset_z = SCEhistograms.at(5)->Interpolate(xx,yy,zz);
    theCalPosOffsets = {offset_x, offset_y, offset_z};
    
  }else if(fRepresentationType == "Parametric"){     
    //this is not supported for parametric
    std::cout << "Change Representation Type to Voxelized TH3 if you want to use the backward offset function" << std::endl;
    theCalPosOffsets.resize(3, 0.0);
  }
  
  return { theCalPosOffsets[0], theCalPosOffsets[1], theCalPosOffsets[2] };
}



// Provides position offsets using a parametric representation
std::vector<double> spacecharge::SpaceChargeSBND::GetPosOffsetsParametric(double xVal, double yVal, double zVal) const
{
    std::vector<double> thePosOffsetsParametric;

    double xValNew = TransformX(xVal);
    double yValNew = TransformY(yVal);
    double zValNew = TransformZ(zVal);

    thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew, yValNew, zValNew, "X"));
    thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew, yValNew, zValNew, "Y"));
    thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew, yValNew, zValNew, "Z"));

    return thePosOffsetsParametric;
}

// Provides one position offset using a parametric representation, for a given axis
double spacecharge::SpaceChargeSBND::GetOnePosOffsetParametric(double xValNew, double yValNew, double zValNew, std::string axis) const
{
    double parA[99][99];
    double parB[99];

    for(int i = 0; i < 99; i++)
        {
            for(int j = 0; j < 99; j++)
                {
                    parA[i][j] = 0.0;
                }
            parB[i] = 0.0;
        }
    if(axis == "X")
        {
            for(int i = 0; i < initialSpatialFitPolN[0] + 1; i++)
                {
                    for(int j = 0; j < intermediateSpatialFitPolN[0] + 1; j++)
                        {
                            parA[i][j] = gSpatialGraphX[i][j]->Eval(zValNew);
                        }
                    intermediateSpatialFitFunctionX[i]->SetParameters(parA[i]);
                }
        }
    else if(axis == "Y")
        {
            for(int i = 0; i < initialSpatialFitPolN[1] + 1; i++)
                {
                    for(int j = 0; j < intermediateSpatialFitPolN[1] + 1; j++)
                        {
                            parA[i][j] = gSpatialGraphY[i][j]->Eval(zValNew);
                        }
                    intermediateSpatialFitFunctionY[i]->SetParameters(parA[i]);
                }
        }
    else if(axis == "Z")
        {
            for(int i = 0; i < initialSpatialFitPolN[2] + 1; i++)
                {
                    for(int j = 0; j < intermediateSpatialFitPolN[2] + 1; j++)
                        {
                            parA[i][j] = gSpatialGraphZ[i][j]->Eval(zValNew);
                        }
                    intermediateSpatialFitFunctionZ[i]->SetParameters(parA[i]);
                }
        }

    double aValNew;
    double bValNew;

    if(axis == "Y")
        {
            aValNew = xValNew;
            bValNew = yValNew;
        }
    else
        {
            aValNew = yValNew;
            bValNew = xValNew;
        }
    double offsetValNew = 0.0;
    if(axis == "X")
        {
            for(int i = 0; i < initialSpatialFitPolN[0] + 1; i++)
                {
                    parB[i] = intermediateSpatialFitFunctionX[i]->Eval(aValNew);
                }
            initialSpatialFitFunctionX->SetParameters(parB);
            offsetValNew = initialSpatialFitFunctionX->Eval(bValNew);
        }
    else if(axis == "Y")
        {
            for(int i = 0; i < initialSpatialFitPolN[1] + 1; i++)
                {
                    parB[i] = intermediateSpatialFitFunctionY[i]->Eval(aValNew);
                }
            initialSpatialFitFunctionY->SetParameters(parB);
            offsetValNew = initialSpatialFitFunctionY->Eval(bValNew);
        }
    else if(axis == "Z")
        {
            for(int i = 0; i < initialSpatialFitPolN[2] + 1; i++)
                {
                    parB[i] = intermediateSpatialFitFunctionZ[i]->Eval(aValNew);
                }
            initialSpatialFitFunctionZ->SetParameters(parB);
            offsetValNew = initialSpatialFitFunctionZ->Eval(bValNew);
        }

    return offsetValNew;
}

// Primary working method of service that provides E field offsets
geo::Vector_t spacecharge::SpaceChargeSBND::GetEfieldOffsets(geo::Point_t const& point) const
{
    std::vector<double> theEfieldOffsets;
    double xx=point.X(), yy=point.Y(), zz=point.Z();
    double offset_x=0., offset_y=0., offset_z=0.;

    if(fRepresentationType == "Voxelized_TH3"){
      //handle OOAV by projecting edge cases
      if(xx<-199.999){xx=-199.999;}
      else if(xx>199.999){xx=199.999;}
      if(yy<-199.999){yy=-199.999;}
      else if(yy>199.999){yy=199.999;}
      if(zz<0.001){zz=0.001;}
      else if(zz>499.999){zz=499.999;}
      offset_x = SCEhistograms.at(6)->Interpolate(xx, yy, zz);
      offset_y = SCEhistograms.at(7)->Interpolate(xx, yy, zz);
      offset_z = SCEhistograms.at(8)->Interpolate(xx, yy, zz);

      theEfieldOffsets = {offset_x, offset_y, offset_z};
      
    }else if(fRepresentationType == "Parametric"){

      if(IsInsideBoundaries(point.X(), point.Y(), point.Z()) == false){
	theEfieldOffsets.resize(3, 0.0);
      }
      else
        {
	  theEfieldOffsets = GetEfieldOffsetsParametric(point.X(), point.Y(), point.Z());

	  // GetOneEfieldOffsetParametric returns V/m
	  // The E-field offsets are returned as -dEx/|E_nominal|, -dEy/|E_nominal|, and -dEz/|E_nominal| where |E_nominal| is DriftField
	  theEfieldOffsets[0] = -1.0 * theEfieldOffsets[0] / (100.0 * DriftField);
	  theEfieldOffsets[1] = -1.0 * theEfieldOffsets[1] / (100.0 * DriftField);
	  theEfieldOffsets[2] = -1.0 * theEfieldOffsets[2] / (100.0 * DriftField);
	}
    }

    return { theEfieldOffsets[0], theEfieldOffsets[1], theEfieldOffsets[2] };
}

// Provides E-field offsets using a parametric representation
std::vector<double> spacecharge::SpaceChargeSBND::GetEfieldOffsetsParametric(double xVal, double yVal, double zVal) const
{
    std::vector<double> theEfieldOffsetsParametric;

    double xValNew = TransformX(xVal);
    double yValNew = TransformY(yVal);
    double zValNew = TransformZ(zVal);

    theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew, yValNew, zValNew, "X"));
    theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew, yValNew, zValNew, "Y"));
    theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew, yValNew, zValNew, "Z"));

    return theEfieldOffsetsParametric;
}

// Provides one E-field offset using a parametric representation, for a given axis
double spacecharge::SpaceChargeSBND::GetOneEfieldOffsetParametric(double xValNew, double yValNew, double zValNew, std::string axis) const
{
    double parA[99][99];
    double parB[99];

    for(int i = 0; i < 99; i++)
        {
            for(int j = 0; j < 99; j++)
                {
                    parA[i][j] = 0.0;
                }
            parB[i] = 0.0;
        }
    if(axis == "X")
        {
            for(int i = 0; i < initialEFieldFitPolN[0] + 1; i++)
                {
                    for(int j = 0; j < intermediateEFieldFitPolN[0] + 1; j++)
                        {
                            parA[i][j] = gEFieldGraphX[i][j]->Eval(zValNew);
                        }
                    intermediateEFieldFitFunctionX[i]->SetParameters(parA[i]);
                }
        }
    else if(axis == "Y")
        {
            for(int i = 0; i < initialEFieldFitPolN[1] + 1; i++)
                {
                    for(int j = 0; j < intermediateEFieldFitPolN[1] + 1; j++)
                        {
                            parA[i][j] = gEFieldGraphY[i][j]->Eval(zValNew);
                        }
                    intermediateEFieldFitFunctionY[i]->SetParameters(parA[i]);
                }
        }
    else if(axis == "Z")
        {
            for(int i = 0; i < initialEFieldFitPolN[2] + 1; i++)
                {
                    for(int j = 0; j < intermediateEFieldFitPolN[2] + 1; j++)
                        {
                            parA[i][j] = gEFieldGraphZ[i][j]->Eval(zValNew);
                        }
                    intermediateEFieldFitFunctionZ[i]->SetParameters(parA[i]);
                }
        }

    double aValNew;
    double bValNew;

    if(axis == "Y")
        {
            aValNew = xValNew;
            bValNew = yValNew;
        }
    else
        {
            aValNew = yValNew;
            bValNew = xValNew;
        }

    double offsetValNew = 0.0;
    if(axis == "X")
        {
            for(int i = 0; i < initialEFieldFitPolN[0] + 1; i++)
                {
                    parB[i] = intermediateEFieldFitFunctionX[i]->Eval(aValNew);
                }
            initialEFieldFitFunctionX->SetParameters(parB);
            offsetValNew = initialEFieldFitFunctionX->Eval(bValNew);
        }
    else if(axis == "Y")
        {
            for(int i = 0; i < initialEFieldFitPolN[1] + 1; i++)
                {
                    parB[i] = intermediateEFieldFitFunctionY[i]->Eval(aValNew);
                }
            initialEFieldFitFunctionY->SetParameters(parB);
            offsetValNew = initialEFieldFitFunctionY->Eval(bValNew);
        }
    else if(axis == "Z")
        {
            for(int i = 0; i < initialEFieldFitPolN[2] + 1; i++)
                {
                    parB[i] = intermediateEFieldFitFunctionZ[i]->Eval(aValNew);
                }
            initialEFieldFitFunctionZ->SetParameters(parB);
            offsetValNew = initialEFieldFitFunctionZ->Eval(bValNew);
        }

    return offsetValNew;
}

// Transform LarSoft-X (cm) to SCE-X (m) coordinate
// [-196.5, 196.5] to [0, 2.0]
double spacecharge::SpaceChargeSBND::TransformX(double xVal) const
{
    xVal = xVal / 100.0;

    // We use the same map twice; 0 is cathod
    // Map [0, 196.5] to [0, 2.0]
    double xValNew = (2.0 / 1.965) * fabs(xVal);

    return xValNew;
}

// Transform LarSoft-Y (cm) to SCE-Y (m) coordinate
// [-200.0, 200.0] to [0, 4.0]
double spacecharge::SpaceChargeSBND::TransformY(double yVal) const
{
    yVal = yVal / 100.0;

    return (yVal + 2.0);
}

// Transform LarSoft-Z (cm) to SCE-Z (m) coordinate
// [0, 500.0] to [0, 5.0]
double spacecharge::SpaceChargeSBND::TransformZ(double zVal) const
{
    return (zVal / 100.0);
}

// Check to see if point is inside boundaries of map
// x = [-196.5, 196.5], y = [-200, 200], and z = [0, 500]  is the active volume boundary
// For now, don't allow to go (slightly) out of range; will change if necessary
bool spacecharge::SpaceChargeSBND::IsInsideBoundaries(double xVal, double yVal, double zVal) const
{
    bool isInside = true;

    if((xVal < -196.5) || (xVal > 196.5) || (yVal < -200.0) || (yVal > 200.0) || (zVal < 0) || (zVal > 500.0))
        {
            isInside = false;
        }

    return isInside;
}
