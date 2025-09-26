# Position Reconstruction using Convolutional Visual Networks (PosRecoCVN)

This directory contains tools for neutrino interaction position reconstruction in SBND using Convolutional Neural Networks (CNNs) trained on photomultiplier tube (PMT) optical response patterns.

## Overview

The PosRecoCVN framework reconstructs the 3D position (x, y, z) of neutrino interactions by analyzing the spatial distribution of light detected by the SBND photomultiplier tube array. The system converts PMT hit patterns into 2D images that are fed to a trained CNN model.

## Directory Structure

```
PosRecoCVN/
├── inference/
│   ├── module/               # LArSoft modules for CNN inference
│   │   ├── PosRecoCVNProducer_module.cc/.hh    # Main inference module
│   │   └── PixelMapVars.h    # Data structure definitions
│   ├── tf/                   # TensorFlow interface
│   └── ReadPixelMapVars.ipynb # Analysis notebook
├── fcls/                     # FHiCL configuration files
│   └── run_pos_inference.fcl  # Main production configuration
├── pmt_maps/                 # PMT channel mapping files
├── training/                 # CNN training script (jupyter notebook)
└── README.md                # This file
```

## Key Components

### PosRecoCVNProducer Module

The main LArSoft module supports **three processing modes** for different data types and use cases:

#### Core Functionality

The module:
1. **Applies mode-specific filters** based on data type and processing requirements
2. **Processes optical data** from PMT hits and flashes with proper time handling
3. **Creates PE matrices** by mapping PMT channels to spatial positions
4. **Generates 2D images** suitable for CNN input
5. **Runs CNN inference** using TensorFlow to predict neutrino interaction positions
6. **Outputs results** in two complementary formats

#### Coordinate System

- **X-axis**: Drift direction, [0, 200] cm
- **Y-axis**: Vertical direction, [-200, 200] cm
- **Z-axis**: Beam direction, [0, 500] cm

#### Processing Modes

1. **MC_testing**: Full simulation analysis with ground truth
   - Processes Monte Carlo (MC) simulation files
   - Calculates ground truth energy barycenter from MCTruth and MCParticles
   - Runs CNN inference and compares with ground truth
   - Applies energy-based filters (dEtpc > 50 MeV + position cuts)
   - Outputs performance metrics and differences

2. **MC_inference**: Simulation-only inference mode
   - Processes Monte Carlo (MC) simulation files for inference only
   - Skips ground truth calculation to save processing time
   - Uses PE-based quality filters instead of energy filters
   - Suitable for large-scale MC processing where truth comparison isn't needed

3. **DATA_inference**: Real data processing mode
   - Processes DATA files for inference
   - No filtering based on ground truth variables nor comparison with ground truth energy barycenter
   - Uses PE-based quality filters instead of energy filters
   - Designed for production data processing

#### Mode-Specific Filtering

Each processing mode applies different event filters:

**MC_testing filters:**
- Neutrino filter (exactly 1 neutrino, configurable)
- Optical data filter (flash with OpHits required)
- Beam window filter (0.367 - 1.9 μs)
- Energy-based filter (dEtpc > 50 MeV + position cuts)

**MC_inference filters:**
- Neutrino filter (exactly 1 neutrino, configurable)
- Optical data filter (flash with OpHits required)
- Beam window filter (0.367 - 1.9 μs)
- PE-based filter (100 PE minimum, 3 channels minimum)

**DATA_inference filters:**
- Neutrino filter (automatically skipped - no MCTruth)
- Optical data filter (flash with OpHits required)
- Beam window filter (0.367 - 1.9 μs)
- PE-based filter (100 PE minimum, 3 channels minimum)

### Output Data Formats

The module produces **two complementary output formats** in the same ROOT file:

#### 1. TTree (`inference_tree`) - For Easy Analysis
Simple, flat structure with one row per event:

**Branch Structure (Mode-Dependent):**
```cpp
// Event identification (all modes)
int run, subrun, event;
bool passed_filters;

// CNN predictions (all modes)
double pred_x, pred_y, pred_z;

// Ground truth coordinates (MC_testing mode only)
double true_x, true_y, true_z;

// Performance metrics (MC_testing mode only)
double diff_x, diff_y, diff_z;     // pred - true
double error_3d;                   // 3D Euclidean distance error

// Physics variables (MC_testing mode only)
double nuv_t, nuv_z;              // Neutrino vertex time and Z
double deposited_energy;          // Total energy deposited in TPC
```

**Note:** Ground truth and performance branches are only filled in MC_testing mode. In MC_inference and DATA_inference modes, these branches are empty by definition and contain default values (-999.0).


#### 2. PixelMapVars Object - Dataproduct for debugging and complete analysis
Complex structure containing:
- **PE matrices** `[event][312_channels]` - Raw PMT response data
- **PE images** `[event][59][70][2]` - 2D images fed to CNN (coated/uncoated PMT maps)
- **PMT maps** - Channel-to-position mapping tables
- **Flash details** - Nested vectors with detailed optical hit information
- **Complete physics data** - All Monte Carlo truth and reconstruction variables


### Model Input Format, Requirements & Performance Metrics

#### Input Format
The CNN expects 2D images with shape `(batch, 59, 70, 2)`:
- **59×70 pixels**: Spatial PMT array mapping
- **2 channels**: Coated and uncoated PMT responses
- **Preprocessing**: PE values normalized to [0,1] range

#### Model Requirements
- **TensorFlow SavedModel format** compatible with LArSoft TensorFlow integration
- **Input shape**: `(batch_size, 59, 70, 2)`
- **Output**: 3 values representing normalized coordinates [0,1] or [-1,1]
- **Scaling**: Model outputs are inverse-scaled to physical coordinates

#### Performance metrics

The system calculates several performance metrics:

- **Individual coordinate errors**: `|pred - true|` for x, y, z
- **3D distance error**: `sqrt(Δx² + Δy² + Δz²)`
- **Coordinate differences**: `pred - true` (signed)


## Configuration

### Processing Mode Selection

Choose the appropriate mode in `run_pos_inference.fcl`:

```fcl
# Processing Mode Configuration
# "MC_testing" = calculate ground truth + prediction (for MC simulation)
# "MC_inference" = only run inference on MC, skip ground truth calculation
# "DATA_inference" = optimized for real data, simplified filters (PE + channels only)
ProcessingMode: "DATA_inference"  # Change as needed
```

### Key Parameters

```fcl
# Processing Mode (REQUIRED)
ProcessingMode: "MC_testing"               # MC_testing, MC_inference, or DATA_inference

# CNN Model Configuration
RunInference: true
ModelPath: "/path/to/saved_model"          # TensorFlow SavedModel directory (auto-discovered)
CustomNormFactor: XXX.XX                   # Image normalization factor

# PMT Mapping Files (auto-discovered)
CoatedPMTMapPath: "coated_pmt_map.csv"
UncoatedPMTMapPath: "uncoated_pmt_map.csv"

# Output Configuration
SavePixelMapVars: true                     # Save detailed PixelMapVars object
```

### Debug Information

Use `Verbosity: 2` or higher to get detailed processing information including:
- Event-by-event filter results
- Energy deposition summaries
- Flash selection logic
- TensorFlow inference details

## Usage

Interactive in console:

```bash
lar -c run_pos_inference.fcl file.root
```

## Contact

- Contact: Sergio Dominguez-Vidales
- Last updated: 26th September 2025