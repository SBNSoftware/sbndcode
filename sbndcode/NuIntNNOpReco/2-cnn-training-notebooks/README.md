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
│   ├── *_position/        ← position model (SavedModel format)
│   └── *_direction_2d/    ← direction model (SavedModel format)
└── plots/                 ← training curves and residual distributions

# TimeCoord_Training.ipynb
temporal/runs/run_YYYYMMDD_HHMM/
├── transformer_tof/       ← transformer model (SavedModel — loaded by LArSoft module)
├── transformer_tof.keras  ← transformer model (Keras — for retraining or fine-tuning)
├── lstm_tof/              ← LSTM model (SavedModel — loaded by LArSoft module)
└── lstm_tof.keras         ← LSTM model (Keras — for retraining or fine-tuning)
```

The LArSoft inference module (`../3-inference-larsoft-module`) loads models in
**SavedModel format** via the TensorFlow C++ API. The Keras format is kept separately
to allow retraining without losing the model architecture and custom layers.

Alongside each SavedModel a JSON metadata file stores the normalisation parameters
needed by the LArSoft C++ inference module.

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
