# Position Reconstruction using Convolutional Visual Networks (PosRecoCVN)

This directory contains tools for neutrino interaction position reconstruction in SBND using Convolutional Neural Networks (CNNs) trained on photomultiplier tube (PMT) response patterns.

## Overview

The PosRecoCVN framework reconstructs the 3D position (x, y, z) of neutrino interactions by analyzing the spatial distribution of light detected by the SBND photomultiplier tube array. The system converts PMT hit patterns into 2D images that are fed to a trained CNN model.

## Directory Structure

```
PosRecoCVN/
├── inference/
│   ├── module/               # LArSoft modules for CNN inference
│   │   ├── SBNDPDSProducer_module.cc/.hh    # Main inference module
│   │   └── PixelMapVars.h    # Data structure definitions
│   ├── tf/                   # TensorFlow interface
│   └── ReadPixelMapVars.ipynb # Analysis notebook
├── fcls/                     # FHiCL configuration files
│   └── run_sbndpds_prod.fcl  # Main production configuration
├── pmt_maps/                 # PMT channel mapping files
├── training/                 # CNN training scripts (if applicable)
└── README.md                # This file
```

## Key Components

### SBNDPDSProducer Module

The main LArSoft module that:
1. **Filters events** based on neutrino truth and flash quality criteria
2. **Processes optical data** from PMT hits and flashes
3. **Creates PE matrices** by mapping photomultiplier channels to spatial positions
4. **Generates 2D images** suitable for CNN input
5. **Runs CNN inference** using TensorFlow to predict neutrino interaction positions
6. **Outputs results** in two complementary formats

### Output Data Formats

The module produces **two complementary output formats** in the same ROOT file:

#### 1. TTree (`inference_tree`) - For Easy Analysis
Simple, flat structure with one row per event:

**Branch Structure:**
```cpp
// Event identification
int run, subrun, event;
bool passed_filters;

// Ground truth coordinates (cm)
double true_x, true_y, true_z;

// CNN predictions (cm)  
double pred_x, pred_y, pred_z;

// Performance metrics
double diff_x, diff_y, diff_z;     // pred - true
double error_3d;                   // 3D Euclidean distance error

// Physics variables
double nuv_t, nuv_z;              // Neutrino vertex time and Z
double deposited_energy;          // Total energy deposited in TPC
```

**Analysis Examples:**
```cpp
// Quick performance analysis
tree->Draw("error_3d");                    // Error distribution
tree->Draw("pred_x:true_x");               // X reconstruction quality
tree->Draw("error_3d", "deposited_energy>100");  // Filter by energy

// Python analysis with uproot
import uproot
data = uproot.open("file.root:inference_tree").arrays()
```

#### 2. PixelMapVars Object - For Detailed Analysis
Complex structure containing:
- **PE matrices** `[event][312_channels]` - Raw PMT response data
- **PE images** `[event][59][70][2]` - 2D images fed to CNN (coated/uncoated PMT maps)
- **PMT maps** - Channel-to-position mapping tables
- **Flash details** - Nested vectors with detailed optical hit information
- **Complete physics data** - All Monte Carlo truth and reconstruction variables

### Model Input Format

The CNN expects 2D images with shape `(batch, 59, 70, 2)`:
- **59×70 pixels**: Spatial PMT array mapping
- **2 channels**: Coated and uncoated PMT responses
- **Preprocessing**: PE values normalized to [0,1] range

### Coordinate System

- **X-axis**: Drift direction, [0, 200] cm (absolute values used for training consistency)
- **Y-axis**: Vertical direction, [-200, 200] cm  
- **Z-axis**: Beam direction, [0, 500] cm

## Configuration

### Verbosity Levels (`run_sbndpds_prod.fcl`)

Configure output detail level for different use cases:

```fcl
Verbosity: 0  # Minimal output for large-scale processing (>1000 events)
              # - Events that pass: "Run=X Event=Y True(x,y,z) Pred(x,y,z) Error3D=Z"
              # - Events that fail: "Run=X Subrun=Y Event=Z FAILED FILTER"

Verbosity: 1  # Basic processing information
              # - Level 0 output plus basic filter status

Verbosity: 2  # Detailed processing information  
              # - Level 1 output plus MCTruth details and energy depositions

Verbosity: >2 # Debug mode (development only)
              # - Full detailed output including TensorFlow debug info
```

### Key Parameters

```fcl
# CNN Model Configuration
RunInference: true
ModelPath: "/path/to/saved_model"          # TensorFlow SavedModel directory
CustomNormFactor: 11311.191                # Image normalization factor

# Event Selection Criteria
MCTruthOrigin: [1, 4]                     # BNB neutrinos and single particles
G4BufferBoxX: [-300, 300]                 # Fiducial volume cuts (cm)
G4BufferBoxY: [-400, 400]
G4BufferBoxZ: [-100, 600]
G4BeamWindow: [-10000, 12000]             # Time window (ns)

# PMT Mapping Files
CoatedPMTMapPath: "path/to/coated_pmt_map.csv"
UncoatedPMTMapPath: "path/to/uncoated_pmt_map.csv"
```

## Usage

### Basic Inference Job

```bash
lar -c run_sbndpds_prod.fcl input_file.root
```

### Large-Scale Processing

```bash
# Set verbosity to 0 for minimal output
# Redirect output to log file
lar -c run_sbndpds_prod.fcl input.root > inference_results.log 2>&1
```

### Analysis

```cpp
// ROOT analysis
TFile f("pixelmap_variables.root");

// Simple analysis with TTree
TTree* tree = f.Get<TTree>("inference_tree");
tree->Draw("error_3d");

// Detailed analysis with PixelMapVars
PixelMapVars* pv;
f.GetObject("PixelMapVars_opanatree__SBNDPDSProd", pv);
// Access PE images: pv->pe_images[event][y][z][map_type]
```

## Performance Metrics

The system calculates several performance metrics:

- **Individual coordinate errors**: `|pred - true|` for x, y, z
- **3D distance error**: `sqrt(Δx² + Δy² + Δz²)`
- **Coordinate differences**: `pred - true` (signed)

## Event Selection Pipeline

1. **MCTruth filtering**: Select BNB neutrinos with specified PDG codes
2. **Flash quality**: Require at least one optical flash with hits
3. **Flash selection**: Choose even/odd PMT groups based on PE content
4. **Energy deposition**: Require >50 MeV deposited energy in TPC
5. **Position cuts**: Ensure true position within detector bounds
6. **Image creation**: Generate normalized PE images for CNN input

## Model Requirements

- **TensorFlow SavedModel format** compatible with LArSoft TensorFlow integration
- **Input shape**: `(batch_size, 59, 70, 2)`
- **Output**: 3 values representing normalized coordinates [0,1] or [-1,1]
- **Scaling**: Model outputs are inverse-scaled to physical coordinates

## Troubleshooting

### Common Issues

1. **Empty output**: Check that input files contain required data products
   - MCTruth with label "generator"
   - OpFlashes with labels "opflashtpc0", "opflashtpc1"
   - OpHits with labels "ophitpmt", "ophitxarapuca"

2. **TensorFlow errors**: Verify model path and compatibility
   - Model must be TensorFlow SavedModel format
   - Check input/output tensor shapes match expectations

3. **Filter failures**: Events failing filters will have `passed_filters = false`
   - Common causes: no neutrino truth, no optical flashes, insufficient energy

### Debug Information

Use `Verbosity: 2` or higher to get detailed processing information including:
- Event-by-event filter results
- Energy deposition summaries
- Flash selection logic
- TensorFlow inference details

## Contact

- Contact: Sergio Dominguez-Vidales
- Last updated: 09/03/2025