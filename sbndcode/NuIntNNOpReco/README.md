# NuIntNNOpReco — Neutrino Interaction Reconstruction with Neural Networks

Tools for neutrino **3D position**, **direction**, and **interaction time** reconstruction
in SBND using neural networks trained on PMT optical response patterns.

---

## Structure

```
NuIntNNOpReco/
├── 1-training-data-preparation/   ← art::EDAnalyzer to extract training data from MC
├── 2-cnn-training-notebooks/      ← Python notebooks to train the networks
├── 3-inference-larsoft-module/    ← art modules to run inference in LArSoft
└── pmt_maps/                      ← PMT channel-to-position mapping files
```

Each subdirectory has its own README with details.

---

## Workflow

```
MC reco1 files
     │
     ▼
1-training-data-preparation    →   training_data.root
     │
     ▼
2-cnn-training-notebooks       →   trained models (SavedModel format)
     │
     ▼
3-inference-larsoft-module     →   reconstruction output (TTree / data products)
```

---

## 1. Training data preparation

`NuIntNNDataPrep_module.cc` reads MC reco1 files and writes `training_data.root`
with four TTrees: MC truth labels, 2D PE images, ophit sequences, and per-flash detail.

```bash
lar -c fcls/run_nuint_nn_dataprep.fcl -s input_mc_reco1.root -n 1000
```

→ See [`1-training-data-preparation/README.md`](1-training-data-preparation/README.md)

---

## 2. CNN training

Two notebooks train the three networks:

- **`PosDirCoords_Training.ipynb`** — ResNet-18 for 3D position (X, Y, Z) and direction (Y, Z projected)
- **`TimeCoord_Training.ipynb`** — Transformer and LSTM for interaction time (nuvT)

→ See [`2-cnn-training-notebooks/README.md`](2-cnn-training-notebooks/README.md)

---

## 3. Inference modules

Three art modules run the trained networks inside LArSoft:

- **`NuIntNNProducer_posdir`** — runs ResNet position and direction models, produces `PixelMapVars` data product
- **`NuSliceFilter`** — selects the best neutrino slice from Pandora (FV + NuScore cut)
- **`NuSliceAnalyzer`** — fills a TTree with Pandora vertex, barycenter, PCA direction, and CNN predictions

Models are loaded via the TensorFlow C++ API from the `tf/` interface and the
`SavedModel` directories bundled under `module/`.

```bash
# MC
lar -c fcls/run_pos_dir_inference.fcl -s input_mc.root

# Data (with Pandora NuScore slice selection)
lar -c fcls/run_pos_inference_data_nuslice.fcl -s input_data.root
```

---

## Contact

Author: Sergio Dominguez-Vidales
