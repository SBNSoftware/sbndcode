# CNN Training Notebooks

Python notebooks for training the neural networks that reconstruct neutrino interaction
properties from SBND optical detector data.

---

## Notebooks

### `PosDirCoords_Training.ipynb` — Position and direction models
Trains two ResNet-18 networks:
- **Position model** (3 outputs: X, Y, Z) — predicts the energy-weighted centroid of the
  neutrino interaction in physical coordinates (cm).
- **Direction model** (2 outputs: dirY, dirZ) — predicts the shower direction projected
  onto the YZ plane using a sign-invariant angular loss.

Both models take as input the 2D PE maps from `images/tree` (shape `[ny, nz, 2]`,
two channels: coated and uncoated PMTs).

### `TimeCoord_Training.ipynb` — Timing model
Trains two sequence models and compares them:
- **Transformer** — self-attention over the full ophit sequence.
- **LSTM** — bidirectional LSTM with exponential-weighted pooling favouring early hits.

Both models take as input the ophit sequences from `temporal/tree`
(shape `[seq_len, 7]`) and predict the neutrino interaction time after
Time-of-Flight correction.

---

## Input

Both notebooks read from the same ROOT file produced by `../1-training-data-preparation`:

```
training_data.root
├── labels/tree    ← MC truth labels
├── images/tree    ← 2D PE maps  (used by PosDirCoords)
└── temporal/tree  ← ophit sequences  (used by TimeCoord)
```

Both notebooks split the dataset **70 / 15 / 15 %** (train / val / test).
Normalisation is always computed on the training set only to avoid data leakage.

---

## Output

Each notebook creates its own timestamped folder:

```
# PosDirCoords_Training.ipynb
position/runs/run_YYYYMMDD_HHMM/
├── models/
│   ├── v<date>_trained_w_<N>_position/
│   │   ├── saved_model/         ← TF SavedModel (loaded by LArSoft C++ module)
│   │   └── inference_config.json
│   └── v<date>_trained_w_<N>_direction_2d/
│       ├── saved_model/
│       └── inference_config.json
└── plots/                       ← training curves and residual distributions

# TimeCoord_Training.ipynb
temporal/runs/run_YYYYMMDD_HHMM/
├── lstm_tof/
│   ├── saved_model/             ← TF SavedModel
│   └── inference_config.json
├── lstm_tof.keras               ← Keras format (for retraining / fine-tuning)
├── transformer_tof/
│   ├── saved_model/
│   └── inference_config.json
├── transformer_tof.keras
├── sc_t.pkl / sc_pe.pkl / sc_tof.pkl   ← MinMaxScalers (required at inference)
└── hist_lstm_tof.json / hist_transformer_tof.json
```

The LArSoft inference module (`../3-inference-larsoft-module`) loads models in
**SavedModel format** via the TensorFlow C++ API. The Keras format is kept separately
to allow retraining without losing the model architecture and custom layers.

Each `inference_config.json` records all parameters needed to reproduce the
preprocessing at inference time (PE normalisation factor, coordinate ranges,
scaler min/max values) together with model architecture info and validation performance.

---

## Current models

The latest production models are kept under `current_models_trained/` in this directory
so that the LArSoft module can reference them with a stable relative path:

```
current_models_trained/
├── v0410_trained_w_388k_position/
│   ├── saved_model/
│   └── inference_config.json   ← pe_normalization_factor, coord ranges, performance
├── v0411_trained_w_388k_direction_2d/
│   ├── saved_model/
│   └── inference_config.json   ← pe_normalization_factor, unit-vector output info
├── v0415_trained_w_388k_time_lstm/
│   ├── saved_model/
│   ├── sc_t.pkl / sc_pe.pkl / sc_tof.pkl
│   └── inference_config.json   ← scaler ranges, feature layout, performance
└── v0415_trained_w_388k_time_transformer/
    ├── saved_model/
    ├── sc_t.pkl / sc_pe.pkl / sc_tof.pkl
    └── inference_config.json
```

When a new training run produces better models, copy the relevant `saved_model/`,
scalers, and `inference_config.json` here and update the FCL model paths in
`../3-inference-larsoft-module/fcls/`.

---

## Dependencies

Tested with the following versions:

| Package | Version |
|---|---|
| tensorflow | 2.16.2 |
| keras | 3.11.2 |
| uproot | 5.3.11 |
| awkward | 2.9.0 |
| numpy | 1.26.4 |
| pandas | 2.3.3 |
| matplotlib | 3.10.8 |
| scipy | 1.15.2 |
| scikit-learn | 1.7.2 |
| joblib | 1.5.3 |

---

## Helper module

`plot_style.py` — shared matplotlib style (colours, fonts, figure sizes) used by
both notebooks.
