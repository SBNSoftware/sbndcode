# PosRecoCVN — Position & Time Reconstruction using Neural Networks

Tools for neutrino interaction **3D position** and **interaction time (nuvT)** reconstruction in SBND using neural networks trained on PMT optical response patterns.

## Directory Structure

```
PosRecoCVN/
├── 1-training-data-preparation/   # art::EDAnalyzer to extract training data from MC
│   ├── module/PosRecoCVNDataPrep_module.cc
│   └── fcls/run_pos_dataprep.fcl
├── 2-cnn-training-notebooks/      # Python training scripts and notebooks
│   ├── nuvT_Reconstruction.ipynb  # Transformer model for nuvT reconstruction
│   ├── preprocess_temporal.py     # Preprocessing pipeline for nuvT training
│   ├── PosRecoCVN_Training.ipynb  # ResNet training for 3D position
│   ├── PosRecoCVN_Inference.ipynb # ResNet inference notebook
│   ├── preprocess_images.py       # Preprocessing for position images
│   ├── pmt_maps/                  # PMT channel-to-position mapping files
│   └── utils.py / config.py
├── 3-inference-larsoft-module/    # art modules for inference in LArSoft
│   ├── module/
│   │   ├── PosRecoCVNProducer_module.cc/.hh  # ResNet position inference
│   │   ├── NuSliceFilter_module.cc           # FV + NuScore filter
│   │   ├── NuSliceAnalyzer_module.cc         # Vertex + CNN predictions TTree
│   │   └── PixelMapVars.h
│   ├── fcls/
│   │   ├── run_pos_dir_inference.fcl
│   │   └── run_pos_inference_data_nuslice.fcl
│   └── tf/                        # TensorFlow interface
└── README.md
```

---

## 1. Training Data Preparation (`1-training-data-preparation/`)

`PosRecoCVNDataPrep_module.cc` is an `art::EDAnalyzer` that reads MC reco1 files and writes a ROOT TTree with:

- **MC truth**: neutrino vertex (x, y, z, t), PDG, energy
- **Labels**: energy-weighted barycenter, PCA direction
- **PE images**: raw (unnormalized) uncoated/coated PMT images `[ny × nz]`
- **OpHit/Flash data**: per-channel PE, time, channel ID

**Usage:**
```bash
lar -c run_pos_dataprep.fcl -s input_mc_reco1.root -n 100
```

Output: `training_data.root` with tree `dataprep/training_tree`.

---

## 2. CNN Training Notebooks (`2-cnn-training-notebooks/`)

### 2a. nuvT Reconstruction — Transformer model

`nuvT_Reconstruction.ipynb` + `preprocess_temporal.py`

Reconstructs the neutrino interaction time **nuvT** from PMT optical hit sequences using a Transformer with exponential time-weighted pooling.

#### Preprocessing pipeline (`preprocess_temporal.py`)

Reads ROOT files and produces NPZ arrays with 7 features per ophit:

| Feature | Description |
|---------|-------------|
| `t_us` | Ophit time [μs], corrected for (1) global ν ToF (`nuvZ/c`) and (2) per-channel photon ToF (`d(vertex,PMT)/v_VUV`) |
| `logPE_norm` | `log1p(PE)/log1p(pe_max)` |
| `det_type` | PMT type: coated (0) or uncoated (1) |
| `x_norm`, `y_norm`, `z_norm` | PMT position normalized to [−1, 1] |
| `delta_t_us` | `t_i − t_first_ophit` (flash shape feature) |

Key constants:
- `v_VUV = 100/7.48 ≈ 13.37 cm/ns` (LAr VUV group velocity, SBND measurement)
- `c_light = 29.98 cm/ns`

**Run preprocessing:**
```python
FORCE_REPROCESS = True   # set True to regenerate NPZ files
```

#### Transformer model

- Architecture: 4× TransformerBlock (d_model=128, heads=8, ff_dim=256) + exponential time-weighted pooling
- Loss: Huber(δ=0.05) — focuses on core Gaussian, down-weights 6% outlier events
- Label: `nuvT_tof = nuvT − nuvZ/c` [ns], normalized to [0,1] with MinMaxScaler

**Performance (test set, MC fall production):**

| Correction | σ [ns] |
|-----------|--------|
| None (raw nuvT) | ~10.9 |
| + global ν ToF (nuvZ/c) | ~6.7 |
| + per-channel photon ToF (d/v_VUV) | ~3.9 |

The BNB RF micro-bunch structure (period ≈ 18.936 ns) becomes visible in the predicted distribution at σ ≈ 4 ns.

**Remaining limitations:** Rayleigh scattering of VUV photons (λ_R ≈ 55 cm in LAr) and PMT transit time spread (~1–2 ns).

### 2b. 3D Position Reconstruction — ResNet

`PosRecoCVN_Training.ipynb` trains a ResNet18 on 2D PMT images `(59×70×2)` to predict the neutrino vertex (x, y, z).

---

## 3. Inference LArSoft Modules (`3-inference-larsoft-module/`)

### PosRecoCVNProducer

`art::EDProducer` that runs the ResNet position model inside LArSoft. Supports three modes:

| Mode | Description |
|------|-------------|
| `MC_testing` | Full MC: ground truth + CNN prediction + residuals |
| `MC_inference` | MC inference only (no ground truth comparison) |
| `DATA_inference` | Real data: PE-based filters, no MC truth |

Output: TTree `inference_tree` + `PixelMapVars` data product.

### NuSliceFilter

`art::EDFilter` — selects the best neutrino slice from Pandora:
- Fiducial volume cut
- `IsClearCosmic == 0`
- `NuScore > threshold` (default 0.6)
- Outputs `bestSliceID` as int data product

### NuSliceAnalyzer

`art::EDAnalyzer` — reads Pandora vertex + SpacePoints + CNN predictions and fills a TTree with vertex, barycenter, PCA direction, and CNN position/direction predictions per slice.

**Usage (data):**
```bash
lar -c run_pos_inference_data_nuslice.fcl -s input_data.root
```

---

## Contact

- Author: Sergio Dominguez-Vidales
- Last updated: March 2026
