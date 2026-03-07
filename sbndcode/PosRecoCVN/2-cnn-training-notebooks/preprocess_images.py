#!/usr/bin/env python3
"""
Preprocess ROOT files to create PE images for CNN training.
Saves images in chunks to avoid memory issues in Jupyter.

Usage:
    python preprocess_images.py --input <root_file> --output <output_dir> [options]

Example:
    python preprocess_images.py -i training.root -o ./preprocessed/ -v
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

# Ensure we can import local modules from THIS directory (copia_de_seguridad)
# This is critical to avoid importing from parent directory
script_dir = os.path.dirname(os.path.abspath(__file__))

# Remove parent directory from sys.path if present to avoid import conflicts
parent_dir = os.path.dirname(script_dir)
if parent_dir in sys.path:
    sys.path.remove(parent_dir)

# Add current directory at the beginning to ensure local imports
if script_dir in sys.path:
    sys.path.remove(script_dir)
sys.path.insert(0, script_dir)

from config import DATA_CONFIG, FILTER_CONFIG, IMAGE_CONFIG
from utils import process_events, create_pe_matrix, create_pe_images


# Network storage prefixes that benefit from local caching
NETWORK_PREFIXES = ('/exp/', '/pnfs/', '/cvmfs/')
# Max file size to copy to local (10 GB) - larger files are read directly
MAX_COPY_SIZE_GB = 10.0


def copy_to_local(file_path, verbose=True):
    """
    Copy file to /tmp for faster I/O if it's on network storage.
    Only copies files smaller than MAX_COPY_SIZE_GB.

    Returns:
        tuple: (local_path, should_cleanup) - local_path to use, whether to delete after
    """
    # Check if file is on network storage
    is_network = any(file_path.startswith(prefix) for prefix in NETWORK_PREFIXES)

    if not is_network:
        return file_path, False

    # Check file size - skip copy for very large files
    file_size_gb = os.path.getsize(file_path) / (1024**3)
    if file_size_gb > MAX_COPY_SIZE_GB:
        if verbose:
            print(f">> File too large to copy ({file_size_gb:.1f} GB > {MAX_COPY_SIZE_GB} GB)", flush=True)
            print(f"   * Reading directly from network storage", flush=True)
        return file_path, False

    # Create local copy in /tmp
    basename = os.path.basename(file_path)
    local_path = os.path.join('/tmp', basename)

    # Check if already cached
    if os.path.exists(local_path):
        local_size = os.path.getsize(local_path)
        remote_size = os.path.getsize(file_path)
        if local_size == remote_size:
            if verbose:
                print(f">> Using cached local copy: {local_path}")
            return local_path, False  # Don't delete cached file

    # Copy to local
    file_size_gb = os.path.getsize(file_path) / (1024**3)
    if verbose:
        print(f">> Copying to local storage for faster I/O...", flush=True)
        print(f"   * Source: {file_path}", flush=True)
        print(f"   * Destination: {local_path}", flush=True)
        print(f"   * Size: {file_size_gb:.2f} GB", flush=True)

    start_time = time.time()
    shutil.copy2(file_path, local_path)
    copy_time = time.time() - start_time

    if verbose:
        speed_mbs = (file_size_gb * 1024) / copy_time
        print(f">> Copy completed in {copy_time:.1f}s ({speed_mbs:.1f} MB/s)", flush=True)

    return local_path, True  # Cleanup after processing


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Preprocess ROOT files to PE images for CNN training',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process training file
  python preprocess_images.py -i training.root -o ./preprocessed/ -v

  # Process only first 10000 events (for testing)
  python preprocess_images.py -i training.root -o ./preprocessed/ -v --max-events 10000
        """
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Input ROOT file path')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory for preprocessed data')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable verbose output')
    parser.add_argument('--max-events', '-n', type=int, default=None,
                        help='Maximum number of events to process (for testing). If not specified, processes all events.')
    return parser.parse_args()


def load_root_data(file_path, verbose=True, max_events=None):
    """
    Load arrays from ROOT file.

    Args:
        file_path: Path to ROOT file
        verbose: Print progress information
        max_events: Maximum number of events to load (None = all events)

    Returns:
        tuple: (file_handle, arrays_dict, channel_dict)
    """
    import concurrent.futures

    if verbose:
        print(f">> Loading ROOT file: {os.path.basename(file_path)}")
        start_time = time.time()

    # Use thread pool for faster decompression (significant speedup for network storage)
    # More workers can help with network I/O bound operations
    n_workers = 16
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=n_workers)

    file = uproot.open(file_path)
    optree = file['opanatree']['OpAnaTree']

    n_entries = optree.num_entries

    # Determine how many events to load
    events_to_load = min(max_events, n_entries) if max_events is not None else n_entries

    if verbose:
        print(f"   * Total events in file: {n_entries:,}")
        if max_events is not None:
            print(f"   * Loading first {events_to_load:,} events (max_events={max_events:,})")

    # Keys to load (must match order expected by process_events)
    keys = DATA_CONFIG['keys_to_load']

    if verbose:
        file_size_gb = os.path.getsize(file_path) / (1024**3)
        print(f"   * File size: {file_size_gb:.1f} GB")
        print(f"   * Loading {len(keys)} arrays with {n_workers} threads...")
        if max_events is None:
            print(f"   * This may take several minutes for large files on network storage...")
        else:
            print(f"   * Loading subset of events - should be faster...")

    # Load arrays at once with threading (faster decompression)
    # Use entry_stop to limit number of events if max_events is specified
    arrays = optree.arrays(keys, library="ak", decompression_executor=executor, entry_stop=events_to_load)

    # Convert to dict format
    arrays = {key: arrays[key] for key in keys}

    # Create channel dictionary from PDSMapTree
    PDSMap = file['opanatree']['PDSMapTree']
    ID = PDSMap['OpDetID'].array()
    Type = PDSMap['OpDetType'].array()
    channel_dict = {int(id_val): int(type_val) for id_val, type_val in zip(ID[0], Type[0])}

    # Cleanup
    executor.shutdown(wait=False)

    if verbose:
        loading_time = time.time() - start_time
        print(f">> Data loaded in {loading_time:.1f} seconds")
        print(f"   * Channel dictionary: {len(channel_dict)} channels")

    return file, arrays, channel_dict


def save_data(output_dir, data, filename, verbose=True):
    """
    Save preprocessed data as compressed NPZ.

    Args:
        output_dir: Directory to save data
        data: Dictionary with arrays to save
        filename: Output filename (without extension)
        verbose: Print progress
    """
    output_file = os.path.join(output_dir, f"{filename}.npz")

    np.savez_compressed(
        output_file,
        images=data['images'],
        dEpromx_abs=data['dEpromx_abs'],
        dEpromy=data['dEpromy'],
        dEpromz=data['dEpromz'],
        dEdirx=data['dEdirx'],
        dEdiry=data['dEdiry'],
        dEdirz=data['dEdirz'],
        dEspreadx=data['dEspreadx'],
        dEspready=data['dEspready'],
        dEspreadz=data['dEspreadz'],
        nuvX=data['nuvX'],
        nuvY=data['nuvY'],
        nuvZ=data['nuvZ'],
        selected_tpc=data['selected_tpc'],
        max_pe_values=data['max_pe_values']  # For later normalization
    )

    if verbose:
        file_size_mb = os.path.getsize(output_file) / (1024 * 1024)
        print(f"   * Saved: {output_file} ({len(data['images']):,} events, {file_size_mb:.1f} MB)")

    return output_file


def save_preprocessing_info(output_dir, info_dict, prefix="", verbose=True):
    """Save preprocessing metadata as JSON."""
    info_file = os.path.join(output_dir, f"{prefix}info.json")

    with open(info_file, 'w') as f:
        json.dump(info_dict, f, indent=2)

    if verbose:
        print(f">> Saved preprocessing info: {info_file}")

    return info_file


def main():
    args = parse_args()

    # Validate input file exists
    if not os.path.exists(args.input):
        print(f"ERROR: Input file not found: {args.input}")
        sys.exit(1)

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    print("=" * 60)
    print("  PosRecoCVN Image Preprocessor")
    print("=" * 60)
    print(f"Input:      {args.input}")
    print(f"Output:     {args.output}")
    print("=" * 60)

    total_start_time = time.time()

    # Step 0: Copy to local storage if on network (much faster I/O)
    local_input, cleanup_local = copy_to_local(args.input, verbose=args.verbose)

    # Step 1: Load ROOT data
    file, arrays, channel_dict = load_root_data(local_input, verbose=args.verbose, max_events=args.max_events)

    # Unpack arrays in the order expected by process_events
    f_ophit_PE = arrays['flash_ophit_pe']
    f_ophit_ch = arrays['flash_ophit_ch']
    f_ophit_t = arrays['flash_ophit_time']
    nuvT = arrays['nuvT']
    dEpromx = arrays['dEpromx']
    dEpromy = arrays['dEpromy']
    dEpromz = arrays['dEpromz']
    dEdirx = arrays['dEdirx']
    dEdiry = arrays['dEdiry']
    dEdirz = arrays['dEdirz']
    dEspreadx = arrays['dEspreadx']
    dEspready = arrays['dEspready']
    dEspreadz = arrays['dEspreadz']
    dEtpc = arrays['dEtpc']
    nuvX = arrays['nuvX']
    nuvY = arrays['nuvY']
    nuvZ = arrays['nuvZ']

    # Step 2: Apply filters
    if args.verbose:
        print("\n>> Applying filters...")

    results = process_events(
        nuvT, f_ophit_PE, f_ophit_ch, f_ophit_t,
        dEpromx, dEpromy, dEpromz, dEtpc, nuvZ,
        channel_dict, FILTER_CONFIG, verbose=args.verbose,
        nuvX=nuvX, nuvY=nuvY,
        dEdirx=dEdirx, dEdiry=dEdiry, dEdirz=dEdirz,
        dEspreadx=dEspreadx, dEspready=dEspready, dEspreadz=dEspreadz
    )

    # Unpack filtered results
    (nuvT_final, f_ophit_PE_final, f_ophit_ch_final, f_ophit_t_final,
     dEpromx_final, dEpromy_final, dEpromz_final, dEtpc_final,
     nuvX_final, nuvY_final, nuvZ_final,
     dEdirx_final, dEdiry_final, dEdirz_final,
     dEspreadx_final, dEspready_final, dEspreadz_final,
     selected_tpc_final, stats) = results

    n_events = len(nuvT_final)
    print(f"\n>> Events after filtering: {n_events:,}")

    if n_events == 0:
        print("ERROR: No events passed filters!")
        sys.exit(1)

    # Step 3: Create PE matrix
    if args.verbose:
        print("\n>> Creating PE matrix...")

    pe_matrix = create_pe_matrix(
        f_ophit_PE_final,
        f_ophit_ch_final,
        IMAGE_CONFIG['max_channels']
    )

    # Step 4: Load PMT maps
    if args.verbose:
        print("\n>> Loading PMT maps...")

    uncoated_map = np.loadtxt(DATA_CONFIG['pmt_maps']['uncoated'], delimiter=",", dtype=int)
    coated_map = np.loadtxt(DATA_CONFIG['pmt_maps']['coated'], delimiter=",", dtype=int)

    if args.verbose:
        print(f"   * Uncoated map shape: {uncoated_map.shape}")
        print(f"   * Coated map shape: {coated_map.shape}")

    # Step 5: Create PE images (raw, unnormalized)
    if args.verbose:
        print("\n>> Creating PE images (raw, unnormalized)...")

    images, max_pe_values = create_pe_images(
        pe_matrix, uncoated_map, coated_map,
        method=IMAGE_CONFIG['selection_method'],
        normalize=False  # Don't normalize - save raw PE values
    )

    print(f">> Images created: {images.shape}")
    print(f"   * Max PE values (for later normalization): {max_pe_values}")

    # Step 6: Prepare coordinate arrays (in cm, not scaled)
    # Convert awkward arrays to numpy - arrays are already flat after process_events
    def to_numpy_flat(arr):
        """Convert awkward array to flat numpy array, flattening only if nested."""
        if arr is None:
            return None
        # Check if array is nested (has variable-length sublists)
        if hasattr(arr, 'ndim') and arr.ndim > 1:
            return np.array(ak.flatten(arr), dtype=np.float32)
        # Already flat or regular array
        return np.array(arr, dtype=np.float32)

    dEpromx_abs = np.abs(to_numpy_flat(dEpromx_final))
    dEpromy_arr = to_numpy_flat(dEpromy_final)
    dEpromz_arr = to_numpy_flat(dEpromz_final)

    # PCA direction variables
    if dEdirx_final is not None:
        dEdirx_arr = to_numpy_flat(dEdirx_final)
        dEdiry_arr = to_numpy_flat(dEdiry_final)
        dEdirz_arr = to_numpy_flat(dEdirz_final)
    else:
        # Fill with NaN if not available
        dEdirx_arr = np.full(n_events, np.nan, dtype=np.float32)
        dEdiry_arr = np.full(n_events, np.nan, dtype=np.float32)
        dEdirz_arr = np.full(n_events, np.nan, dtype=np.float32)

    # PCA spread variables
    if dEspreadx_final is not None:
        dEspreadx_arr = to_numpy_flat(dEspreadx_final)
        dEspready_arr = to_numpy_flat(dEspready_final)
        dEspreadz_arr = to_numpy_flat(dEspreadz_final)
    else:
        # Fill with NaN if not available
        dEspreadx_arr = np.full(n_events, np.nan, dtype=np.float32)
        dEspready_arr = np.full(n_events, np.nan, dtype=np.float32)
        dEspreadz_arr = np.full(n_events, np.nan, dtype=np.float32)

    # Neutrino positions (for debugging)
    if nuvX_final is not None:
        nuvX_arr = to_numpy_flat(nuvX_final)
        nuvY_arr = to_numpy_flat(nuvY_final)
        nuvZ_arr = to_numpy_flat(nuvZ_final)
    else:
        # Fill with NaN if not available
        nuvX_arr = np.full(n_events, np.nan, dtype=np.float32)
        nuvY_arr = np.full(n_events, np.nan, dtype=np.float32)
        nuvZ_arr = np.full(n_events, np.nan, dtype=np.float32)

    # Selected TPC
    selected_tpc_arr = np.array(selected_tpc_final, dtype=np.int8)

    # Step 7: Save data to single file
    input_basename = os.path.basename(args.input)
    input_name = os.path.splitext(input_basename)[0]  # Remove .root

    if args.verbose:
        print(f"\n>> Saving {n_events:,} events...")

    data = {
        'images': images,
        'dEpromx_abs': dEpromx_abs,
        'dEpromy': dEpromy_arr,
        'dEpromz': dEpromz_arr,
        'dEdirx': dEdirx_arr,
        'dEdiry': dEdiry_arr,
        'dEdirz': dEdirz_arr,
        'dEspreadx': dEspreadx_arr,
        'dEspready': dEspready_arr,
        'dEspreadz': dEspreadz_arr,
        'nuvX': nuvX_arr,
        'nuvY': nuvY_arr,
        'nuvZ': nuvZ_arr,
        'selected_tpc': selected_tpc_arr,
        'max_pe_values': np.array(max_pe_values, dtype=np.float32)
    }

    output_file = save_data(args.output, data, input_name, verbose=args.verbose)

    # Step 8: Save preprocessing info
    total_time = time.time() - total_start_time

    preprocessing_info = {
        'input_file': os.path.abspath(args.input),
        'output_file': os.path.basename(output_file),
        'total_events': int(n_events),
        'image_shape': list(images.shape[1:]),
        'max_pe_values': [float(f) for f in max_pe_values],  # Max PE for later normalization
        'filter_config': {
            'min_energy_tpc': FILTER_CONFIG['min_energy_tpc'],
            'depromx_range': FILTER_CONFIG['depromx_range'],
            'depromy_range': FILTER_CONFIG['depromy_range'],
            'depromz_range': FILTER_CONFIG['depromz_range'],
            'active_volume': FILTER_CONFIG.get('active_volume', None)
        },
        'image_config': {
            'max_channels': IMAGE_CONFIG['max_channels'],
            'selection_method': IMAGE_CONFIG['selection_method']
        },
        'coordinate_ranges': {
            'dEpromx_abs': [float(np.min(dEpromx_abs)), float(np.max(dEpromx_abs))],
            'dEpromy': [float(np.min(dEpromy_arr)), float(np.max(dEpromy_arr))],
            'dEpromz': [float(np.min(dEpromz_arr)), float(np.max(dEpromz_arr))],
            'dEdirx': [float(np.nanmin(dEdirx_arr)), float(np.nanmax(dEdirx_arr))],
            'dEdiry': [float(np.nanmin(dEdiry_arr)), float(np.nanmax(dEdiry_arr))],
            'dEdirz': [float(np.nanmin(dEdirz_arr)), float(np.nanmax(dEdirz_arr))],
            'dEspreadx': [float(np.nanmin(dEspreadx_arr)), float(np.nanmax(dEspreadx_arr))],
            'dEspready': [float(np.nanmin(dEspready_arr)), float(np.nanmax(dEspready_arr))],
            'dEspreadz': [float(np.nanmin(dEspreadz_arr)), float(np.nanmax(dEspreadz_arr))]
        },
        'processing_time_seconds': total_time
    }

    save_preprocessing_info(args.output, preprocessing_info, prefix=f"{input_name}_", verbose=args.verbose)

    # Cleanup local copy if needed
    if cleanup_local and os.path.exists(local_input):
        os.remove(local_input)
        if args.verbose:
            print(f">> Cleaned up local copy: {local_input}")

    # Summary
    print("\n" + "=" * 60)
    print("  Preprocessing Complete!")
    print("=" * 60)
    print(f"Total events:  {n_events:,}")
    print(f"Total time:    {total_time:.1f} seconds")
    print(f"Output:        {output_file}")
    print("=" * 60)


# Utility function for loading in Jupyter
def load_preprocessed_data(npz_file_or_dir):
    """
    Load preprocessed data from NPZ file.

    Args:
        npz_file_or_dir: Path to .npz file or directory containing .npz files

    Returns:
        tuple: (images, coords, pca_dir, pca_spread, info)
            - images: numpy array (N, H, W, C)
            - coords: numpy array (N, 3) with [x_abs, y, z] in cm
            - pca_dir: numpy array (N, 3) with [dirx, diry, dirz] PCA direction
            - pca_spread: numpy array (N, 3) with [spreadx, spready, spreadz] PCA spread
            - info: dict with preprocessing metadata (or None if no info file)
    """
    import glob

    # Handle both file and directory inputs
    if os.path.isdir(npz_file_or_dir):
        npz_files = sorted(glob.glob(os.path.join(npz_file_or_dir, "*.npz")))
        if not npz_files:
            raise FileNotFoundError(f"No .npz files found in {npz_file_or_dir}")
        npz_file = npz_files[0]
        info_files = glob.glob(os.path.join(npz_file_or_dir, "*_info.json"))
        info_file = info_files[0] if info_files else None
    else:
        npz_file = npz_file_or_dir
        info_file = npz_file.replace('.npz', '_info.json')
        if not os.path.exists(info_file):
            info_file = None

    # Load data
    data = np.load(npz_file)
    images = data['images']
    coords = np.column_stack([
        data['dEpromx_abs'],
        data['dEpromy'],
        data['dEpromz']
    ])

    # Load PCA variables
    pca_dir = np.column_stack([
        data['dEdirx'],
        data['dEdiry'],
        data['dEdirz']
    ])

    pca_spread = np.column_stack([
        data['dEspreadx'],
        data['dEspready'],
        data['dEspreadz']
    ])

    # Load metadata if available
    info = None
    if info_file and os.path.exists(info_file):
        with open(info_file, 'r') as f:
            info = json.load(f)

    return images, coords, pca_dir, pca_spread, info


if __name__ == '__main__':
    main()
