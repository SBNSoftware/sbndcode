#!/usr/bin/env python3
"""
Preprocess ROOT files for nuvT (neutrino time) temporal reconstruction.

Saves padded ophit sequences as NPZ files ready for LSTM/Transformer training.
Mirrors the structure of preprocess_images.py for positional reconstruction.

Output per event: (num_flashes * max_ophits_per_flash, 6) feature tensor
  Features: [t_norm, PE_norm, det_type, x_norm, y_norm, z_norm]
  Padding marker: -1000.0 (use MaskNegative1000 layer in model)

Usage:
    python preprocess_temporal.py -i <root_file> -o <output_dir> [options]

Example:
    python preprocess_temporal.py \\
        -i /exp/sbnd/data/.../training.root \\
        -o /exp/sbnd/data/.../preprocessed_temporal/training \\
        -v
"""

import argparse
import json
import os
import shutil
import sys
import time
import numpy as np
import uproot
import awkward as ak
import concurrent.futures

# Ensure local imports (config, utils) come from this directory
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
if parent_dir in sys.path:
    sys.path.remove(parent_dir)
if script_dir in sys.path:
    sys.path.remove(script_dir)
sys.path.insert(0, script_dir)

from config import DATA_CONFIG, FILTER_CONFIG
from utils import process_events

# Base temporal filter config — max_flashes_for_keep_all set dynamically in main()
# based on --num-flashes: if 1 → always select max-PE TPC; if >1 → keep up to N TPCs
TEMPORAL_FILTER_CONFIG_BASE = {**FILTER_CONFIG}

# ── Geometry normalization constants (from PDSMapTree) ────────────────────────
GEOM = {
    'x_max':  213.75,   # cm  → normalize to [-1, 1] by / x_max
    'y_max':  175.00,   # cm  → normalize to [-1, 1] by / y_max
    'z_min':   16.05,   # cm  \
    'z_max':  484.95,   # cm  /  → normalize to [0, 1]
}

# Padding marker (must match MaskNegative1000 layer in model)
PAD_VALUE = -1000.0

# Speed of light in vacuum (cm/ns) — neutrino ToF correction: nuvZ / C_LIGHT_CM_NS
C_LIGHT_CM_NS = 29.9792458  # cm/ns

# VUV scintillation light group velocity in LAr (SBND measurement: 1/vg = 7.48 ± 0.14 ns/m)
V_VUV_LIGHT_CM_NS = 100.0 / 7.48   # ≈ 13.37 cm/ns

# Network storage prefixes that benefit from local caching
NETWORK_PREFIXES = ('/exp/', '/pnfs/', '/cvmfs/')
MAX_COPY_SIZE_GB = 10.0


# ── I/O helpers ───────────────────────────────────────────────────────────────

def copy_to_local(file_path, verbose=True):
    """Copy file to /tmp for faster I/O if on network storage."""
    is_network = any(file_path.startswith(p) for p in NETWORK_PREFIXES)
    if not is_network:
        return file_path, False

    file_size_gb = os.path.getsize(file_path) / (1024**3)
    if file_size_gb > MAX_COPY_SIZE_GB:
        if verbose:
            print(f">> File too large to copy ({file_size_gb:.1f} GB), reading directly")
        return file_path, False

    basename = os.path.basename(file_path)
    local_path = os.path.join('/tmp', basename)

    if os.path.exists(local_path) and os.path.getsize(local_path) == os.path.getsize(file_path):
        if verbose:
            print(f">> Using cached local copy: {local_path}")
        return local_path, False

    if verbose:
        print(f">> Copying to /tmp for faster I/O ({file_size_gb:.2f} GB)...")
    t0 = time.time()
    shutil.copy2(file_path, local_path)
    if verbose:
        speed = file_size_gb * 1024 / (time.time() - t0)
        print(f"   * Done in {time.time()-t0:.1f}s ({speed:.0f} MB/s)")
    return local_path, True


def parse_args():
    parser = argparse.ArgumentParser(
        description='Preprocess ROOT files for nuvT temporal reconstruction'
    )
    parser.add_argument('--input',  '-i', required=True,  help='Input ROOT file')
    parser.add_argument('--output', '-o', required=True,  help='Output directory')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--max-events', '-n', type=int, default=None,
                        help='Max events to process (for testing)')
    parser.add_argument('--max-ophits', type=int, default=100,
                        help='Max ophits per flash (default: 100)')
    parser.add_argument('--num-flashes', type=int, default=1,
                        help='Number of flashes per event after padding (default: 1, TPC with max PE)')
    return parser.parse_args()


# ── ROOT loading ──────────────────────────────────────────────────────────────

TEMPORAL_KEYS = [
    'flash_ophit_pe', 'flash_ophit_ch', 'flash_ophit_time',
    'nuvT', 'nuvX', 'nuvY', 'nuvZ',
    'dEpromx', 'dEpromy', 'dEpromz',
    'dEdirx', 'dEdiry', 'dEdirz',
    'dEspreadx', 'dEspready', 'dEspreadz',
    'dEtpc',
]


def load_root_data(file_path, verbose=True, max_events=None):
    """Load arrays and channel geometry from ROOT file."""
    if verbose:
        print(f">> Loading ROOT file: {os.path.basename(file_path)}")
        t0 = time.time()

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=16)
    f = uproot.open(file_path)
    tree = f['opanatree']['OpAnaTree']

    n_entries = tree.num_entries
    n_load = min(max_events, n_entries) if max_events is not None else n_entries
    if verbose:
        print(f"   * Events in file: {n_entries:,}  |  Loading: {n_load:,}")

    arrays = tree.arrays(
        TEMPORAL_KEYS, library="ak",
        decompression_executor=executor,
        entry_stop=n_load
    )
    arrays = {k: arrays[k] for k in TEMPORAL_KEYS}

    # Build channel geometry dict: ch_id → (type, x, y, z)
    pds = f['opanatree']['PDSMapTree']
    ids   = np.array(pds['OpDetID'].array()[0],   dtype=int)
    types = np.array(pds['OpDetType'].array()[0],  dtype=int)
    xs    = np.array(pds['OpDetX'].array()[0],     dtype=np.float32)
    ys    = np.array(pds['OpDetY'].array()[0],     dtype=np.float32)
    zs    = np.array(pds['OpDetZ'].array()[0],     dtype=np.float32)
    channel_geom = {
        int(i): (int(t), float(x), float(y), float(z))
        for i, t, x, y, z in zip(ids, types, xs, ys, zs)
    }

    # Minimal channel_dict needed by process_events (type only)
    channel_dict = {ch: t for ch, (t, x, y, z) in channel_geom.items()}

    executor.shutdown(wait=False)

    if verbose:
        print(f"   * Loaded in {time.time()-t0:.1f}s  |  {len(channel_geom)} PDS channels")

    return f, arrays, channel_dict, channel_geom


# ── Temporal preprocessing ────────────────────────────────────────────────────

# Beam gate window (μs) — ophits outside this window are noise/cosmics
BEAM_GATE_MIN_US = -0.1   # μs  (small negative margin for pre-trigger)
BEAM_GATE_MAX_US =  1.8   # μs  (slightly beyond 1.6 μs beam window)


def filter_and_sort_ophits(f_ophit_t, f_ophit_PE, f_ophit_ch):
    """
    Keep only ophits within the beam gate [BEAM_GATE_MIN_US, BEAM_GATE_MAX_US],
    then sort by time ascending within each flash.
    """
    time_mask = (f_ophit_t >= BEAM_GATE_MIN_US) & (f_ophit_t <= BEAM_GATE_MAX_US)
    f_t_gated  = f_ophit_t[time_mask]
    f_PE_gated = f_ophit_PE[time_mask]
    f_ch_gated = f_ophit_ch[time_mask]

    idx = ak.argsort(f_t_gated, axis=-1, ascending=True)
    return f_t_gated[idx], f_PE_gated[idx], f_ch_gated[idx]


def pad_sequence(arr_1d, target_len, pad_value=PAD_VALUE):
    """Pad or truncate a 1-D array to target_len."""
    n = len(arr_1d)
    if n >= target_len:
        return np.array(arr_1d[:target_len], dtype=np.float32)
    out = np.full(target_len, pad_value, dtype=np.float32)
    out[:n] = arr_1d
    return out


def build_padded_arrays(f_ophit_t, f_ophit_PE, f_ophit_ch,
                        max_ophits, num_flashes):
    """
    Build fixed-shape arrays (N, num_flashes, max_ophits) for t, PE, ch.
    Padding entries filled with PAD_VALUE.
    """
    n_events = len(f_ophit_t)
    seq_len = num_flashes * max_ophits

    t_out  = np.full((n_events, seq_len), PAD_VALUE, dtype=np.float32)
    pe_out = np.full((n_events, seq_len), PAD_VALUE, dtype=np.float32)
    ch_out = np.full((n_events, seq_len), PAD_VALUE, dtype=np.float32)

    for ev in range(n_events):
        flashes_t  = f_ophit_t[ev]
        flashes_pe = f_ophit_PE[ev]
        flashes_ch = f_ophit_ch[ev]
        n_fl = min(len(flashes_t), num_flashes)

        for fl in range(n_fl):
            offset = fl * max_ophits
            t_out[ev, offset:offset+max_ophits]  = pad_sequence(flashes_t[fl],  max_ophits)
            pe_out[ev, offset:offset+max_ophits] = pad_sequence(flashes_pe[fl], max_ophits)
            ch_out[ev, offset:offset+max_ophits] = pad_sequence(flashes_ch[fl], max_ophits)

    return t_out, pe_out, ch_out


def build_feature_tensor(t_arr, pe_arr, ch_arr, channel_geom, pe_max):
    """
    Build (N, seq_len, 7) feature tensor.

    Features per ophit:
        0: t_us         raw μs  (no normalization — preserves linear correlation with nuvT)
        1: logPE_norm   [0, 1]  log1p(PE) / log1p(pe_max) — compresses skewed PE distribution
        2: det_type     {0,1,2,3}  integer
        3: x_norm       [-1, 1]
        4: y_norm       [-1, 1]
        5: z_norm       [0, 1]
        6: delta_t_us   time relative to first real ophit in event (μs) — encodes flash shape

    Padding entries (t == PAD_VALUE) keep PAD_VALUE across all 7 features.
    """
    n_events, seq_len = t_arr.shape
    features = np.full((n_events, seq_len, 7), PAD_VALUE, dtype=np.float32)

    z_span = GEOM['z_max'] - GEOM['z_min']

    # Precompute lookup arrays indexed by channel id (0..311)
    max_ch = max(channel_geom.keys()) + 1
    lut_type = np.full(max_ch, PAD_VALUE, dtype=np.float32)
    lut_x    = np.full(max_ch, PAD_VALUE, dtype=np.float32)
    lut_y    = np.full(max_ch, PAD_VALUE, dtype=np.float32)
    lut_z    = np.full(max_ch, PAD_VALUE, dtype=np.float32)

    for ch, (det_type, x, y, z) in channel_geom.items():
        lut_type[ch] = float(det_type)
        lut_x[ch]    = x / GEOM['x_max']
        lut_y[ch]    = y / GEOM['y_max']
        lut_z[ch]    = (z - GEOM['z_min']) / z_span

    # Build mask of real (non-padding) entries
    mask = (t_arr != PAD_VALUE)  # (N, seq_len)

    # Feature 0: raw time in μs — NOT normalized
    t_norm = np.where(mask, t_arr, PAD_VALUE)

    # Feature 1: log-PE normalized to [0, 1] — log1p compresses the skewed PE distribution
    log_pe_max = np.log1p(pe_max)
    pe_norm = np.where(mask, np.log1p(np.maximum(pe_arr, 0.0)) / log_pe_max, PAD_VALUE)
    pe_norm = np.clip(pe_norm, 0.0, 1.0)

    # Channel lookup (only for real entries)
    ch_int = ch_arr.astype(int)
    ch_int_safe = np.clip(ch_int, 0, max_ch - 1)

    det_type_arr = np.where(mask, lut_type[ch_int_safe], PAD_VALUE)
    x_arr        = np.where(mask, lut_x[ch_int_safe],    PAD_VALUE)
    y_arr        = np.where(mask, lut_y[ch_int_safe],    PAD_VALUE)
    z_arr        = np.where(mask, lut_z[ch_int_safe],    PAD_VALUE)

    # Feature 6: Δt = t_i - t_first_ophit (μs)
    # t_first = global minimum time across all real ophits per event (earliest photon arrival)
    t_safe  = np.where(mask, t_arr, np.inf)
    first_t = np.min(t_safe, axis=1, keepdims=True)         # (N, 1)
    first_t = np.where(first_t == np.inf, 0.0, first_t)    # fallback for all-pad events
    delta_t = np.where(mask, t_arr - first_t, PAD_VALUE)

    features[:, :, 0] = t_norm
    features[:, :, 1] = pe_norm
    features[:, :, 2] = det_type_arr
    features[:, :, 3] = x_arr
    features[:, :, 4] = y_arr
    features[:, :, 5] = z_arr
    features[:, :, 6] = delta_t

    return features, mask


# ── Save / info ───────────────────────────────────────────────────────────────

def save_data(output_dir, data, filename, verbose=True):
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, f"{filename}.npz")
    np.savez_compressed(out_file, **data)
    if verbose:
        size_mb = os.path.getsize(out_file) / 1024**2
        print(f"   * Saved: {out_file}  ({data['nuvT'].shape[0]:,} events, {size_mb:.1f} MB)")
    return out_file


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: Input file not found: {args.input}")
        sys.exit(1)

    os.makedirs(args.output, exist_ok=True)
    print("=" * 60)
    print("  PosRecoCVN Temporal Preprocessor (nuvT)")
    print("=" * 60)
    print(f"Input:       {args.input}")
    print(f"Output:      {args.output}")
    print(f"Max ophits:  {args.max_ophits} per flash")
    print(f"Num flashes: {args.num_flashes} per event")
    print(f"Seq length:  {args.num_flashes * args.max_ophits} total ophits")
    print("=" * 60)

    t_total = time.time()

    # Step 0: local copy for faster I/O
    local_input, cleanup = copy_to_local(args.input, verbose=args.verbose)

    # Step 1: load ROOT
    _, arrays, channel_dict, channel_geom = load_root_data(
        local_input, verbose=args.verbose, max_events=args.max_events
    )

    # Step 2: filter (same as position reconstruction)
    # max_flashes_for_keep_all: 0 → always pick max-PE flash (1 TPC)
    #                           N → keep up to N flashes (both TPCs when N≥2)
    temporal_filter_config = {
        **TEMPORAL_FILTER_CONFIG_BASE,
        'max_flashes_for_keep_all': args.num_flashes,
    }
    if args.verbose:
        print(f"\n>> Applying filters (max_flashes_for_keep_all={args.num_flashes})...")

    results = process_events(
        arrays['nuvT'],
        arrays['flash_ophit_pe'],
        arrays['flash_ophit_ch'],
        arrays['flash_ophit_time'],
        arrays['dEpromx'], arrays['dEpromy'], arrays['dEpromz'],
        arrays['dEtpc'],
        arrays['nuvZ'],
        channel_dict, temporal_filter_config, verbose=args.verbose,
        nuvX=arrays['nuvX'], nuvY=arrays['nuvY'],
        dEdirx=arrays['dEdirx'], dEdiry=arrays['dEdiry'], dEdirz=arrays['dEdirz'],
        dEspreadx=arrays['dEspreadx'], dEspready=arrays['dEspready'],
        dEspreadz=arrays['dEspreadz'],
    )

    (nuvT_f, f_pe_f, f_ch_f, f_t_f,
     dEpromx_f, dEpromy_f, dEpromz_f, dEtpc_f,
     nuvX_f, nuvY_f, nuvZ_f,
     dEdirx_f, dEdiry_f, dEdirz_f,
     dEspreadx_f, dEspready_f, dEspreadz_f,
     selected_tpc_f, _stats) = results

    n_events = len(nuvT_f)
    print(f"\n>> Events after filtering: {n_events:,}")
    if n_events == 0:
        print("ERROR: No events passed filters!")
        sys.exit(1)

    # Step 3: filter ophits to beam gate + sort by time
    if args.verbose:
        print(f"\n>> Filtering ophits to beam gate [{BEAM_GATE_MIN_US}, {BEAM_GATE_MAX_US}] μs and sorting...")
    f_t_sorted, f_pe_sorted, f_ch_sorted = filter_and_sort_ophits(f_t_f, f_pe_f, f_ch_f)

    # Step 4: pad to fixed shape (N, num_flashes * max_ophits)
    if args.verbose:
        print(f"\n>> Padding to ({n_events:,}, {args.num_flashes * args.max_ophits}) ...")
    t0 = time.time()
    t_arr, pe_arr, ch_arr = build_padded_arrays(
        f_t_sorted, f_pe_sorted, f_ch_sorted,
        args.max_ophits, args.num_flashes
    )
    if args.verbose:
        print(f"   * Done in {time.time()-t0:.1f}s")

    # Step 4b: ToF-inside correction — subtract nuvZ/c from ophit times
    # Same correction applied to nuvT: removes neutrino flight time inside detector (z=0 → vertex)
    # nuvZ in cm, C_LIGHT_CM_NS in cm/ns → correction in ns → /1000 to convert to μs
    nuvZ_np = np.array(ak.flatten(nuvZ_f), dtype=np.float32)
    tof_correction_us = (nuvZ_np / C_LIGHT_CM_NS / 1000.0)[:, np.newaxis]  # (N, 1) μs
    real_mask = (t_arr != PAD_VALUE)
    t_arr = np.where(real_mask, t_arr - tof_correction_us, PAD_VALUE)
    if args.verbose:
        t_corr_min = float(t_arr[real_mask].min())
        t_corr_max = float(t_arr[real_mask].max())
        print(f"\n>> ToF-inside correction applied (nuvZ/c, c={C_LIGHT_CM_NS} cm/ns)")
        print(f"   * t_ophit corrected range: [{t_corr_min:.4f}, {t_corr_max:.4f}] μs")

    # Step 4c: Per-channel photon ToF correction (ΔTγ per channel)
    # For each ophit j at channel ch_j: t_j -= d(vertex, PMT_j) / v_VUV
    # This removes the geometry-dependent photon propagation delay, leaving only
    #   nuvT_tof + scintillation_emission_delay  (spread ~τ_fast ≈ 6 ns)
    # Training: uses MC true vertex (nuvX, nuvY, nuvZ).
    # Inference: replace with PosRecoCVN predicted vertex (position uncertainty ~5-10 cm → δt < 1 ns).
    if args.verbose:
        print(f"\n>> Per-channel photon ToF correction "
              f"(ΔTγ_j = d(vtx,PMT_j)/v_VUV, v_VUV={V_VUV_LIGHT_CM_NS:.2f} cm/ns)...")

    nuvX_np = np.array(ak.flatten(nuvX_f), dtype=np.float32)  # (N,) cm — vertex X
    nuvY_np = np.array(ak.flatten(nuvY_f), dtype=np.float32)  # (N,) cm — vertex Y
    # nuvZ_np already extracted in Step 4b

    # PMT position lookup table indexed by channel id
    max_ch_id = max(channel_geom.keys()) + 1
    lut_pmt_x = np.zeros(max_ch_id, dtype=np.float32)
    lut_pmt_y = np.zeros(max_ch_id, dtype=np.float32)
    lut_pmt_z = np.zeros(max_ch_id, dtype=np.float32)
    for ch_id, (_, px, py, pz) in channel_geom.items():
        lut_pmt_x[ch_id] = px
        lut_pmt_y[ch_id] = py
        lut_pmt_z[ch_id] = pz

    real_mask_tof = (t_arr != PAD_VALUE)                               # (N, seq_len)
    ch_safe       = np.clip(ch_arr.astype(int), 0, max_ch_id - 1)     # safe channel indices

    # Vertex → PMT distance per ophit slot: (N, seq_len)
    dx = nuvX_np[:, np.newaxis] - lut_pmt_x[ch_safe]
    dy = nuvY_np[:, np.newaxis] - lut_pmt_y[ch_safe]
    dz = nuvZ_np[:, np.newaxis] - lut_pmt_z[ch_safe]
    d_vtx_pmt = np.sqrt(dx**2 + dy**2 + dz**2)                        # cm

    perchan_tof_us = d_vtx_pmt / V_VUV_LIGHT_CM_NS / 1000.0           # μs
    t_arr = np.where(real_mask_tof, t_arr - perchan_tof_us, PAD_VALUE)

    if args.verbose:
        real_after = (t_arr != PAD_VALUE)
        t_min = float(t_arr[real_after].min())
        t_max = float(t_arr[real_after].max())
        print(f"   * t_ophit after ΔTγ correction range: [{t_min:.4f}, {t_max:.4f}] μs")
        print(f"   * Spread (expect << 1.6 μs; only scintillation delay remains): "
              f"{(t_max - t_min)*1000:.1f} ns")

    # Step 5: normalization ranges (time is kept raw in μs, only PE is normalized)
    real_mask  = (t_arr != PAD_VALUE)
    pe_max_val = float(pe_arr[real_mask].max())

    if args.verbose:
        t_data_min = float(t_arr[real_mask].min())
        t_data_max = float(t_arr[real_mask].max())
        print(f"\n>> Normalization ranges:")
        print(f"   * t   : raw μs  (data range: [{t_data_min:.4f}, {t_data_max:.4f}] μs)")
        print(f"   * PE  : [0, {pe_max_val:.2f}]")

    # Step 6: build 6-feature tensor (N, seq_len, 6)
    if args.verbose:
        print(f"\n>> Building feature tensor (N, {args.num_flashes * args.max_ophits}, 6)...")
    t0 = time.time()
    features, mask = build_feature_tensor(
        t_arr, pe_arr, ch_arr, channel_geom,
        pe_max=pe_max_val
    )
    if args.verbose:
        print(f"   * Done in {time.time()-t0:.1f}s  |  shape: {features.shape}")

    # Step 7: flatten nuvT (single neutrino per event after filters)
    nuvT_np   = np.array(ak.flatten(nuvT_f),  dtype=np.float32)
    # nuvX_np, nuvY_np already extracted in Step 4c; nuvZ_np in Step 4b
    tpc_np    = np.array(selected_tpc_f,       dtype=np.int8)

    # nuvT with ToF-inside correction (same nuvZ/c subtracted from ophit times)
    nuvT_tof_np = nuvT_np - nuvZ_np / C_LIGHT_CM_NS  # ns

    # Step 8: save
    if args.verbose:
        print(f"\n>> Saving {n_events:,} events...")

    input_name = os.path.splitext(os.path.basename(args.input))[0]

    data = {
        'ophit_features': features,     # (N, seq_len, 6)  float32
        'ophit_mask':     mask,          # (N, seq_len)     bool  — True = real ophit
        'nuvT':           nuvT_np,       # (N,)  raw ns
        'nuvT_tof':       nuvT_tof_np,   # (N,)  ns  — ToF-inside corrected (nuvT - nuvZ/c)
        'nuvX':           nuvX_np,       # (N,)  cm
        'nuvY':           nuvY_np,       # (N,)  cm
        'nuvZ':           nuvZ_np,       # (N,)  cm
        'selected_tpc':   tpc_np,        # (N,)  0 or 1
    }
    out_file = save_data(args.output, data, input_name, verbose=args.verbose)

    # Step 9: save preprocessing metadata
    info = {
        'input_file':    os.path.abspath(args.input),
        'output_file':   os.path.basename(out_file),
        'total_events':  int(n_events),
        'seq_len':       int(args.num_flashes * args.max_ophits),
        'num_flashes':   int(args.num_flashes),
        'max_ophits':    int(args.max_ophits),
        'n_features':    7,
        'features':      ['t_us', 'logPE_norm', 'det_type', 'x_norm', 'y_norm', 'z_norm', 'delta_t_us'],
        'pad_value':     PAD_VALUE,
        'normalization': {
            't_note':      'μs, corrected for (1) nuvZ/c global ν-ToF [step 4b] and '
                           '(2) per-channel d(vtx,PMT)/v_VUV photon ToF [step 4c, v_VUV='
                           f'{V_VUV_LIGHT_CM_NS:.2f} cm/ns]; nuvT scaling at training time',
            'pe_max':      pe_max_val,
            'x_range_cm':  [-GEOM['x_max'], GEOM['x_max']],
            'y_range_cm':  [-GEOM['y_max'], GEOM['y_max']],
            'z_range_cm':  [GEOM['z_min'],  GEOM['z_max']],
        },
        'nuvT_range_ns': [float(nuvT_np.min()), float(nuvT_np.max())],
        'processing_time_seconds': time.time() - t_total,
    }
    info_file = os.path.join(args.output, f"{input_name}_info.json")
    with open(info_file, 'w') as f_json:
        json.dump(info, f_json, indent=2)
    if args.verbose:
        print(f">> Saved preprocessing info: {info_file}")

    # Cleanup
    if cleanup and os.path.exists(local_input):
        os.remove(local_input)

    print("\n" + "=" * 60)
    print("  Temporal Preprocessing Complete!")
    print("=" * 60)
    print(f"Events:      {n_events:,}")
    print(f"Output:      {out_file}")
    print(f"Total time:  {time.time()-t_total:.1f}s")
    print("=" * 60)


if __name__ == '__main__':
    main()
