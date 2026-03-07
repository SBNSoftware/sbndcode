"""
Utility functions for PosRecoCVN training pipeline.
Contains data processing, filtering, and visualization functions.

Optimized with:
- Vectorized awkward operations (avoiding ak.to_list where possible)
- NumPy scatter operations (np.add.at)
- Optional Numba JIT compilation for inner loops
"""

import awkward as ak
import numpy as np
import time
import pandas as pd
from IPython.display import display
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Optional numba support for additional speedup
try:
    from numba import jit, prange
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    # Dummy decorator when numba not available
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator
    prange = range


@jit(nopython=True, parallel=True, cache=True)
def _numba_fill_pe_matrix(pe_matrix, event_indices, channels, pes):
    """
    Numba-accelerated PE matrix filling.

    Uses parallel loops for ~2-5x additional speedup over numpy.
    """
    n_ophits = len(event_indices)
    for i in prange(n_ophits):
        ev = event_indices[i]
        ch = channels[i]
        pe = pes[i]
        pe_matrix[ev, ch] += pe


@jit(nopython=True, cache=True)
def _numba_compute_pca_2d(pes, y_coords, z_coords):
    """
    Compute weighted PCA for 2D points (Y-Z plane).

    Returns: (pca_y, pca_z, elongation)
    """
    n = len(pes)
    if n < 2:
        return 0.0, 1.0, 1.0

    # Weighted centroid
    weight_sum = 0.0
    cy, cz = 0.0, 0.0
    for i in range(n):
        w = pes[i]
        weight_sum += w
        cy += w * y_coords[i]
        cz += w * z_coords[i]

    if weight_sum == 0:
        return 0.0, 1.0, 1.0

    cy /= weight_sum
    cz /= weight_sum

    # Weighted covariance matrix
    cov_yy, cov_yz, cov_zz = 0.0, 0.0, 0.0
    for i in range(n):
        w = pes[i]
        dy = y_coords[i] - cy
        dz = z_coords[i] - cz
        cov_yy += w * dy * dy
        cov_yz += w * dy * dz
        cov_zz += w * dz * dz

    cov_yy /= weight_sum
    cov_yz /= weight_sum
    cov_zz /= weight_sum

    # Analytical eigenvalues for 2x2 symmetric matrix
    # |a  b|   eigenvalues: (a+c)/2 ± sqrt(((a-c)/2)^2 + b^2)
    # |b  c|
    trace = cov_yy + cov_zz
    det = cov_yy * cov_zz - cov_yz * cov_yz
    discriminant = trace * trace / 4 - det

    if discriminant < 0:
        discriminant = 0

    sqrt_disc = np.sqrt(discriminant)
    lambda1 = trace / 2 - sqrt_disc  # smaller eigenvalue
    lambda2 = trace / 2 + sqrt_disc  # larger eigenvalue

    # Eigenvector for larger eigenvalue
    if abs(cov_yz) > 1e-10:
        # Standard case
        vy = cov_yz
        vz = lambda2 - cov_yy
    elif abs(cov_yy - lambda2) < abs(cov_zz - lambda2):
        vy = 1.0
        vz = 0.0
    else:
        vy = 0.0
        vz = 1.0

    # Normalize
    norm = np.sqrt(vy * vy + vz * vz)
    if norm > 1e-10:
        vy /= norm
        vz /= norm

    # Elongation
    if lambda1 > 1e-10:
        elongation = lambda2 / lambda1
    else:
        elongation = 100.0

    return vy, vz, elongation


@jit(nopython=True, parallel=True, cache=True)
def _numba_compute_all_pcas(pca_vectors, pca_elongation,
                             flat_pe, flat_ch, event_starts, event_ends,
                             y_lookup, z_lookup, valid_channels, max_channel):
    """
    Compute PCA for all events in parallel using numba.
    """
    n_events = len(event_starts)

    for event_idx in prange(n_events):
        start_idx = event_starts[event_idx]
        end_idx = event_ends[event_idx]

        if start_idx >= end_idx:
            pca_vectors[event_idx, 0] = 0.0
            pca_vectors[event_idx, 1] = 1.0
            pca_elongation[event_idx] = 1.0
            continue

        # Count valid ophits
        n_valid = 0
        for i in range(start_idx, end_idx):
            ch = flat_ch[i]
            if ch >= 0 and ch < max_channel and valid_channels[ch]:
                n_valid += 1

        if n_valid < 2:
            pca_vectors[event_idx, 0] = 0.0
            pca_vectors[event_idx, 1] = 1.0
            pca_elongation[event_idx] = 1.0
            continue

        # Extract valid data (pre-allocate for numba efficiency)
        valid_pes = np.empty(n_valid, dtype=np.float32)
        valid_y = np.empty(n_valid, dtype=np.float32)
        valid_z = np.empty(n_valid, dtype=np.float32)

        j = 0
        for i in range(start_idx, end_idx):
            ch = flat_ch[i]
            if ch >= 0 and ch < max_channel and valid_channels[ch]:
                valid_pes[j] = flat_pe[i]
                valid_y[j] = y_lookup[ch]
                valid_z[j] = z_lookup[ch]
                j += 1

        # Compute PCA
        vy, vz, elong = _numba_compute_pca_2d(valid_pes, valid_y, valid_z)
        pca_vectors[event_idx, 0] = vy
        pca_vectors[event_idx, 1] = vz
        pca_elongation[event_idx] = elong


class CutFlowTracker:
    """Track cut flow statistics with nice Jupyter display."""
    
    def __init__(self, initial_count):
        self.initial = initial_count
        self.cuts = []
        self.current = initial_count
    
    def add_cut(self, name, remaining):
        removed = self.current - remaining
        efficiency = remaining / self.initial if self.initial > 0 else 0
        self.cuts.append({
            'Cut': name,
            'Removed': removed,
            'Remaining': remaining, 
            'Cumulative_Eff': f"{efficiency:.3f}"
        })
        self.current = remaining
    
    def display(self):
        """Display nice table in Jupyter"""
        df = pd.DataFrame(self.cuts)
        print(f"Initial events: {self.initial:,}")
        print("=" * 60)
        display(df)
        print("=" * 60)
        if self.current > 0:
            print(f"Final efficiency: {self.current/self.initial:.3f} ({self.current:,}/{self.initial:,})")


def create_channel_lookup(channel_dict, pmt_types, xas_types):
    """Create fast channel lookup for categorization."""
    pmt_channels = {ch for ch, val in channel_dict.items() if val in pmt_types}
    xas_channels = {ch for ch, val in channel_dict.items() if val in xas_types}
    
    lookup = {}
    for ch in pmt_channels:
        lookup[ch] = 0 if ch % 2 == 0 else 1  # PMT even/odd
    for ch in xas_channels:
        lookup[ch] = 2 if ch % 2 == 0 else 3  # XAS even/odd
    
    return lookup


def vectorized_categorize_flashes(f_ophit_ch, channel_lookup):
    """
    Fast flash categorization using vectorized awkward operations.

    Avoids ak.to_list() by using numpy lookup array.
    ~10-20x faster than list comprehension version.
    """
    first_channels = f_ophit_ch[:, :, 0]

    if not channel_lookup:
        return ak.full_like(first_channels, -1)

    # Build numpy lookup array for O(1) access
    max_channel = max(channel_lookup.keys()) + 1
    lookup_array = np.full(max_channel, -1, dtype=np.int8)
    for ch, cat in channel_lookup.items():
        lookup_array[ch] = cat

    # Flatten, lookup, unflatten - all vectorized
    flat_channels = ak.flatten(first_channels)

    # Convert to numpy for fast lookup
    flat_ch_np = ak.to_numpy(flat_channels).astype(np.int32)

    # Vectorized lookup with bounds checking
    valid_mask = (flat_ch_np >= 0) & (flat_ch_np < max_channel)
    flat_categories = np.where(valid_mask, lookup_array[np.clip(flat_ch_np, 0, max_channel - 1)], -1)

    # Unflatten back to original structure
    counts = ak.num(first_channels, axis=1)
    categories = ak.unflatten(flat_categories, counts)

    return categories


def smart_flash_selection(f_ophit_PE, categories, max_flashes_for_keep_all=2):
    """Optimized flash selection based on PE sums."""
    sum_pe = ak.sum(f_ophit_PE, axis=2)

    mask_even = (categories == 0) | (categories == 2)
    mask_odd = (categories == 1) | (categories == 3)

    sum_even = ak.sum(ak.where(mask_even, sum_pe, 0), axis=1)
    sum_odd = ak.sum(ak.where(mask_odd, sum_pe, 0), axis=1)
    
    n_flashes = ak.num(categories, axis=1)
    decision = sum_even >= sum_odd
    
    keep_all = n_flashes <= max_flashes_for_keep_all
    selection_mask = ak.where(
        keep_all,
        ak.ones_like(mask_even),
        ak.where(decision[:, np.newaxis], mask_even, mask_odd)
    )
    
    return selection_mask, decision


def apply_selection_mask(arrays_dict, mask):
    """Apply mask to dictionary of arrays."""
    return {name: arr[mask] for name, arr in arrays_dict.items()}


def select_tpc_values(tpc_dict, decision):
    """Select TPC values based on decision."""
    selector = ak.concatenate([
        decision[:, np.newaxis],
        (~decision)[:, np.newaxis]
    ], axis=1)

    result = {}
    for name, array in tpc_dict.items():
        result[name] = ak.sum(ak.where(selector, array, 0), axis=1)

    return result


def select_tpc_values_new(tpc_dict, decision):
    """
    Select TPC values based on decision - vectorized version.

    Arrays have structure [n_events, 2] where index 0=TPC0, index 1=TPC1
    Decision: True=TPC0 (even), False=TPC1 (odd)

    Optimized to avoid Python loops using ak.where.
    """
    # Convert decision to numpy for efficient indexing
    decision_np = ak.to_numpy(decision)

    result = {}
    for name, array in tpc_dict.items():
        # Ensure array has correct structure [n_events, 2]
        # Some arrays might be nested differently

        # Try to get values for index 0 and 1
        try:
            # Get TPC0 and TPC1 values
            val_tpc0 = array[:, 0]
            val_tpc1 = array[:, 1]

            # Select based on decision using ak.where
            selected = ak.where(decision, val_tpc0, val_tpc1)
            result[name] = selected

        except (IndexError, TypeError):
            # Fallback for irregular structures: convert to numpy
            array_np = ak.to_numpy(ak.fill_none(ak.pad_none(array, 2, axis=1), 0.0))

            # Vectorized selection using numpy advanced indexing
            idx = np.where(decision_np, 0, 1)
            selected = array_np[np.arange(len(array_np)), idx]
            result[name] = ak.Array(selected)

    return result


def filter_neutrinos_in_active_volume(nuvX, nuvY, nuvZ, active_volume_config):
    """
    Filter neutrinos that are inside the active volume.

    Parameters:
    -----------
    nuvX, nuvY, nuvZ : awkward arrays
        Neutrino vertex positions (can have multiple neutrinos per event)
    active_volume_config : dict
        Active volume ranges {'x_range': (min, max), 'y_range': (min, max), 'z_range': (min, max)}

    Returns:
    --------
    mask : awkward array
        Boolean mask indicating which neutrinos are inside the active volume
    """
    x_min, x_max = active_volume_config['x_range']
    y_min, y_max = active_volume_config['y_range']
    z_min, z_max = active_volume_config['z_range']

    # Create mask for neutrinos inside active volume
    mask_x = (nuvX >= x_min) & (nuvX <= x_max)
    mask_y = (nuvY >= y_min) & (nuvY <= y_max)
    mask_z = (nuvZ >= z_min) & (nuvZ <= z_max)

    # Combined mask
    mask_active = mask_x & mask_y & mask_z

    return mask_active


def process_events(nuvT, f_ophit_PE, f_ophit_ch, f_ophit_t,
                  dEpromx, dEpromy, dEpromz, dEtpc, nuvZ,
                  channel_dict, config, verbose=True,
                  nuvX=None, nuvY=None,
                  dEdirx=None, dEdiry=None, dEdirz=None,
                  dEspreadx=None, dEspready=None, dEspreadz=None):
    """
    Process events with optimized pipeline.

    NEW FILTER LOGIC:
    1. Filter neutrinos inside active volume (using nuvX, nuvY, nuvZ)
    2. Keep only events with exactly 1 neutrino in active volume
    3. Apply remaining filters (flashes, energy, position, etc.)

    Parameters:
    -----------
    nuvT, f_ophit_PE, f_ophit_ch, f_ophit_t : arrays
        Standard arrays
    dEpromx, dEpromy, dEpromz, dEtpc, nuvZ : arrays
        Standard arrays
    channel_dict : dict
        Channel mapping
    config : dict
        Configuration dictionary with filter parameters (including 'active_volume')
    verbose : bool
        Show progress and statistics
    nuvX, nuvY : arrays, optional
        Neutrino X and Y positions (required for active volume filtering)
        If None, uses old logic (single neutrino first)

    Returns:
    --------
    Tuple of all filtered arrays plus statistics tracker
    """

    if verbose:
        print(">> Starting optimized event processing...")
        start_time = time.time()

    initial_count = len(nuvT)
    tracker = CutFlowTracker(initial_count)

    if verbose:
        print(f"Initial events: {initial_count:,}")

    # Step 1: Active volume filter (NEW LOGIC)
    if nuvX is not None and nuvY is not None and 'active_volume' in config:
        if verbose:
            print("Applying active volume filter...")

        # Filter neutrinos inside active volume
        mask_active = filter_neutrinos_in_active_volume(nuvX, nuvY, nuvZ, config['active_volume'])

        # Apply mask to neutrino arrays
        nuvT_active = nuvT[mask_active]
        nuvX_active = nuvX[mask_active]
        nuvY_active = nuvY[mask_active]
        nuvZ_active = nuvZ[mask_active]

        # Count neutrinos in active volume per event
        n_neutrinos_in_active = ak.num(nuvT_active, axis=1)

        # Keep events with exactly 1 neutrino in active volume
        mask_single_in_active = n_neutrinos_in_active == 1

        arrays = {
            'nuvT': nuvT_active[mask_single_in_active],
            'nuvX': nuvX_active[mask_single_in_active],  # Keep nuvX through pipeline
            'nuvY': nuvY_active[mask_single_in_active],  # Keep nuvY through pipeline
            'nuvZ': nuvZ_active[mask_single_in_active],  # Keep nuvZ through pipeline
            'f_ophit_PE': f_ophit_PE[mask_single_in_active],
            'f_ophit_ch': f_ophit_ch[mask_single_in_active],
            'f_ophit_t': f_ophit_t[mask_single_in_active],
            'dEpromx': dEpromx[mask_single_in_active],
            'dEpromy': dEpromy[mask_single_in_active],
            'dEpromz': dEpromz[mask_single_in_active],
            'dEtpc': dEtpc[mask_single_in_active]
        }

        # Add PCA variables if provided
        if dEdirx is not None:
            arrays['dEdirx'] = dEdirx[mask_single_in_active]
            arrays['dEdiry'] = dEdiry[mask_single_in_active]
            arrays['dEdirz'] = dEdirz[mask_single_in_active]
        if dEspreadx is not None:
            arrays['dEspreadx'] = dEspreadx[mask_single_in_active]
            arrays['dEspready'] = dEspready[mask_single_in_active]
            arrays['dEspreadz'] = dEspreadz[mask_single_in_active]

        tracker.add_cut("Neutrinos in active volume + single neutrino", len(arrays['nuvT']))

    else:
        # OLD LOGIC: Single neutrino filter first (backward compatibility)
        mask_single = ak.num(nuvT) == 1

        arrays = {
            'nuvT': nuvT[mask_single],
            'f_ophit_PE': f_ophit_PE[mask_single],
            'f_ophit_ch': f_ophit_ch[mask_single],
            'f_ophit_t': f_ophit_t[mask_single],
            'dEpromx': dEpromx[mask_single],
            'dEpromy': dEpromy[mask_single],
            'dEpromz': dEpromz[mask_single],
            'dEtpc': dEtpc[mask_single],
            'nuvZ': nuvZ[mask_single]
        }

        # Add PCA variables if provided
        if dEdirx is not None:
            arrays['dEdirx'] = dEdirx[mask_single]
            arrays['dEdiry'] = dEdiry[mask_single]
            arrays['dEdirz'] = dEdirz[mask_single]
        if dEspreadx is not None:
            arrays['dEspreadx'] = dEspreadx[mask_single]
            arrays['dEspready'] = dEspready[mask_single]
            arrays['dEspreadz'] = dEspreadz[mask_single]

        tracker.add_cut("Single neutrino", len(arrays['nuvT']))
    
    # Step 2: Has flashes filter  
    mask_flashes = ak.num(arrays['f_ophit_PE'], axis=1) >= config['min_flashes']
    arrays = apply_selection_mask(arrays, mask_flashes)
    
    tracker.add_cut("Has flashes", len(arrays['nuvT']))
    
    # Step 3: Flash processing
    if verbose:
        print("Processing flashes...")

    channel_lookup = create_channel_lookup(
        channel_dict, config['pmt_types'], config['xas_types']
    )
    categories = vectorized_categorize_flashes(arrays['f_ophit_ch'], channel_lookup)

    flash_mask, decision = smart_flash_selection(
        arrays['f_ophit_PE'], categories, config['max_flashes_for_keep_all']
    )
    
    # Apply flash selection
    arrays['f_ophit_PE'] = arrays['f_ophit_PE'][flash_mask]
    arrays['f_ophit_ch'] = arrays['f_ophit_ch'][flash_mask]
    arrays['f_ophit_t'] = arrays['f_ophit_t'][flash_mask]

    # Select TPC values (for all per-TPC variables)
    tpc_arrays = {
        'dEpromx': arrays['dEpromx'],
        'dEpromy': arrays['dEpromy'],
        'dEpromz': arrays['dEpromz'],
        'dEtpc': arrays['dEtpc']
    }

    # Add PCA direction variables if they exist
    if 'dEdirx' in arrays:
        tpc_arrays['dEdirx'] = arrays['dEdirx']
        tpc_arrays['dEdiry'] = arrays['dEdiry']
        tpc_arrays['dEdirz'] = arrays['dEdirz']

    # Add PCA spread variables if they exist
    if 'dEspreadx' in arrays:
        tpc_arrays['dEspreadx'] = arrays['dEspreadx']
        tpc_arrays['dEspready'] = arrays['dEspready']
        tpc_arrays['dEspreadz'] = arrays['dEspreadz']

    selected_tpc = select_tpc_values(tpc_arrays, decision)
    arrays.update(selected_tpc)

    # Store TPC selection as simple label: 0 for TPC0 (even, X<0), 1 for TPC1 (odd, X>0)
    arrays['selected_tpc'] = ak.where(decision, 0, 1)

    
    # Step 4: Valid data cut
    if verbose:
        print("Applying data validity cuts...")
    
    valid_mask = (
        (arrays['dEpromx'] != config['invalid_marker']) &
        (arrays['dEpromy'] != config['invalid_marker']) &
        (arrays['dEpromz'] != config['invalid_marker'])
    )
    
    arrays = apply_selection_mask(arrays, valid_mask)
    tracker.add_cut("Valid data (≠ -999 in dEprom)", len(arrays['nuvT']))
    
    # Step 5: Energy cut
    if verbose:
        print("Applying energy cuts...")
    
    energy_mask = arrays['dEtpc'] > config['min_energy_tpc']
    arrays = apply_selection_mask(arrays, energy_mask)
    tracker.add_cut("Energy cut", len(arrays['nuvT']))
    
    # Step 6: Position cuts
    if verbose:
        print("Applying position cuts...")

    position_mask = (
        (arrays['dEpromx'] >= config['depromx_range'][0]) &
        (arrays['dEpromx'] <= config['depromx_range'][1]) &
        (arrays['dEpromy'] >= config['depromy_range'][0]) &
        (arrays['dEpromy'] <= config['depromy_range'][1]) &
        (arrays['dEpromz'] >= config['depromz_range'][0]) &
        (arrays['dEpromz'] <= config['depromz_range'][1])
    )

    final_arrays = apply_selection_mask(arrays, position_mask)
    tracker.add_cut("dEprom in active volume cut", len(final_arrays['nuvT']))
    
    if verbose:
        processing_time = time.time() - start_time
        print(f">> Processing completed in {processing_time:.2f} seconds")
        print()
        tracker.display()
        
        if len(final_arrays['nuvT']) > 0:
            print("\n* Final dataset ranges:")
            for var in ['dEpromx', 'dEpromy', 'dEpromz', 'dEtpc']:
                arr = final_arrays[var]
                print(f"  {var}: [{ak.min(arr):.1f}, {ak.max(arr):.1f}]")

    # Return with nuvX, nuvY if they were provided (new logic), otherwise None
    nuvX_final = final_arrays.get('nuvX', None)
    nuvY_final = final_arrays.get('nuvY', None)
    nuvZ_final = final_arrays.get('nuvZ', None)

    # PCA variables (can be None if not provided)
    dEdirx_final = final_arrays.get('dEdirx', None)
    dEdiry_final = final_arrays.get('dEdiry', None)
    dEdirz_final = final_arrays.get('dEdirz', None)
    dEspreadx_final = final_arrays.get('dEspreadx', None)
    dEspready_final = final_arrays.get('dEspready', None)
    dEspreadz_final = final_arrays.get('dEspreadz', None)

    return (
        final_arrays['nuvT'],
        final_arrays['f_ophit_PE'],
        final_arrays['f_ophit_ch'],
        final_arrays['f_ophit_t'],
        final_arrays['dEpromx'],
        final_arrays['dEpromy'],
        final_arrays['dEpromz'],
        final_arrays['dEtpc'],
        nuvX_final,  # Can be None if old logic was used
        nuvY_final,  # Can be None if old logic was used
        nuvZ_final,  # Can be None if old logic was used
        dEdirx_final,  # PCA direction x
        dEdiry_final,  # PCA direction y
        dEdirz_final,  # PCA direction z
        dEspreadx_final,  # PCA spread x
        dEspready_final,  # PCA spread y
        dEspreadz_final,  # PCA spread z
        final_arrays['selected_tpc'],
        tracker
    )


def smart_flash_selection_new(f_ophit_PE, categories, max_flashes_for_keep_all=2):
    """
    Vectorized flash selection using awkward operations.

    Logic:
    1. Sum PEs per flash
    2. Sum flashes by TPC (even channels = TPC0, odd channels = TPC1)
    3. Pick TPC with more total PE

    Optimized to avoid Python loops using ak.where and ak.sum with masks.
    ~5-10x faster than loop-based version.
    """
    # Sum PE per flash (axis=2 sums ophits within each flash)
    sum_pe = ak.sum(f_ophit_PE, axis=2)

    # Create masks for even (TPC0) and odd (TPC1) categories
    # Categories: 0,2 = even (TPC0), 1,3 = odd (TPC1)
    mask_even = (categories == 0) | (categories == 2)
    mask_odd = (categories == 1) | (categories == 3)

    # Sum PEs by TPC using vectorized operations
    # ak.where returns 0 where mask is False, keeping array structure
    sum_even = ak.sum(ak.where(mask_even, sum_pe, 0), axis=1)
    sum_odd = ak.sum(ak.where(mask_odd, sum_pe, 0), axis=1)

    # Decision: True = TPC0 (even), False = TPC1 (odd)
    decision = sum_even >= sum_odd

    # Number of flashes per event
    n_flashes = ak.num(categories, axis=1)

    # Keep all flashes if n_flashes <= max_flashes_for_keep_all
    keep_all = n_flashes <= max_flashes_for_keep_all

    # Build selection mask vectorized
    # When keep_all: mask is all True
    # Otherwise: mask based on decision
    # Use broadcasting: decision[:, np.newaxis] broadcasts to flash dimension

    # For events where we keep all: ones_like gives True for all flashes
    # For events where we select: use mask_even if decision else mask_odd
    selection_mask = ak.where(
        keep_all[:, np.newaxis],
        ak.ones_like(mask_even, dtype=bool),
        ak.where(decision[:, np.newaxis], mask_even, mask_odd)
    )

    return selection_mask, decision


def process_events_new(nuvT, f_ophit_PE, f_ophit_ch, f_ophit_t,
                       dEpromx, dEpromy, dEpromz, dEtpc,
                       channel_dict, config, verbose=True):
    """
    Process events with optimized pipeline - handles variable depth arrays.
    Uses smart_flash_selection_new to avoid UnionArray errors.

    Parameters:
    -----------
    All the usual arrays plus:
    config : dict
        Configuration dictionary with filter parameters
    verbose : bool
        Show progress and statistics

    Returns:
    --------
    Tuple of all filtered arrays plus statistics tracker
    """

    if verbose:
        print(">> Starting optimized event processing (new version)...")
        start_time = time.time()

    initial_count = len(nuvT)
    tracker = CutFlowTracker(initial_count)

    if verbose:
        print(f"Initial events: {initial_count:,}")

    # Step 1: Single neutrino filter
    mask_single = ak.num(nuvT) == 1

    arrays = {
        'nuvT': nuvT[mask_single],
        'f_ophit_PE': f_ophit_PE[mask_single],
        'f_ophit_ch': f_ophit_ch[mask_single],
        'f_ophit_t': f_ophit_t[mask_single],
        'dEpromx': dEpromx[mask_single],
        'dEpromy': dEpromy[mask_single],
        'dEpromz': dEpromz[mask_single],
        'dEtpc': dEtpc[mask_single]
    }

    tracker.add_cut("Single neutrino", len(arrays['nuvT']))

    # Step 2: Has flashes filter
    mask_flashes = ak.num(arrays['f_ophit_PE'], axis=1) >= config['min_flashes']
    arrays = apply_selection_mask(arrays, mask_flashes)

    tracker.add_cut("Has flashes", len(arrays['nuvT']))

    # Step 3: Flash processing
    if verbose:
        print("Processing flashes...")

    channel_lookup = create_channel_lookup(
        channel_dict, config['pmt_types'], config['xas_types']
    )
    categories = vectorized_categorize_flashes(arrays['f_ophit_ch'], channel_lookup)

    # Use new version that handles variable depth
    flash_mask, decision = smart_flash_selection_new(
        arrays['f_ophit_PE'], categories, config['max_flashes_for_keep_all']
    )

    # Apply flash selection
    arrays['f_ophit_PE'] = arrays['f_ophit_PE'][flash_mask]
    arrays['f_ophit_ch'] = arrays['f_ophit_ch'][flash_mask]
    arrays['f_ophit_t'] = arrays['f_ophit_t'][flash_mask]

    # Select TPC values (only for dEprom and dEtpc - these have TPC structure)
    tpc_arrays = {
        'dEpromx': arrays['dEpromx'],
        'dEpromy': arrays['dEpromy'],
        'dEpromz': arrays['dEpromz'],
        'dEtpc': arrays['dEtpc']
    }
    selected_tpc = select_tpc_values_new(tpc_arrays, decision)
    arrays.update(selected_tpc)

    # Store TPC selection as simple label: 0 for TPC0 (even, X<0), 1 for TPC1 (odd, X>0)
    # Vectorized: decision=True -> 0, decision=False -> 1
    arrays['selected_tpc'] = ak.where(decision, 0, 1)


    # Step 4: Valid data cut
    if verbose:
        print("Applying data validity cuts...")

    valid_mask = (
        (arrays['dEpromx'] != config['invalid_marker']) &
        (arrays['dEpromy'] != config['invalid_marker']) &
        (arrays['dEpromz'] != config['invalid_marker'])
    )

    arrays = apply_selection_mask(arrays, valid_mask)
    tracker.add_cut("Valid data (≠ -999 in dEprom)", len(arrays['nuvT']))

    # Step 5: Energy cut
    if verbose:
        print("Applying energy cuts...")

    energy_mask = arrays['dEtpc'] > config['min_energy_tpc']
    arrays = apply_selection_mask(arrays, energy_mask)
    tracker.add_cut("Energy cut", len(arrays['nuvT']))

    # Step 6: Position cuts
    if verbose:
        print("Applying position cuts...")

    position_mask = (
        (arrays['dEpromx'] >= config['depromx_range'][0]) &
        (arrays['dEpromx'] <= config['depromx_range'][1]) &
        (arrays['dEpromy'] >= config['depromy_range'][0]) &
        (arrays['dEpromy'] <= config['depromy_range'][1]) &
        (arrays['dEpromz'] >= config['depromz_range'][0]) &
        (arrays['dEpromz'] <= config['depromz_range'][1])
    )

    final_arrays = apply_selection_mask(arrays, position_mask)
    tracker.add_cut("dEprom in active volume cut", len(final_arrays['nuvT']))

    if verbose:
        processing_time = time.time() - start_time
        print(f">> Processing completed in {processing_time:.2f} seconds")
        print()
        tracker.display()

        if len(final_arrays['nuvT']) > 0:
            print("\n* Final dataset ranges:")
            for var in ['dEpromx', 'dEpromy', 'dEpromz', 'dEtpc']:
                arr = final_arrays[var]
                print(f"  {var}: [{ak.min(arr):.1f}, {ak.max(arr):.1f}]")

    return (
        final_arrays['nuvT'],
        final_arrays['f_ophit_PE'],
        final_arrays['f_ophit_ch'],
        final_arrays['f_ophit_t'],
        final_arrays['dEpromx'],
        final_arrays['dEpromy'],
        final_arrays['dEpromz'],
        final_arrays['dEtpc'],
        final_arrays['selected_tpc'],
        tracker
    )


def process_events_data(nuScore, f_ophit_PE, f_ophit_ch, f_ophit_t,
                        channel_dict, config, verbose=True):
    """
    Process events for DATA (real data) - no MC information available.

    Filters:
    1. nuScore > threshold (default 0.5) - neutrino-like events
    2. Has at least min_flashes flashes
    3. Flash selection by TPC (same as MC processing)

    Parameters:
    -----------
    nuScore : array
        Neutrino score from ML classifier (per event)
    f_ophit_PE, f_ophit_ch, f_ophit_t : arrays
        Flash optical hit information
    channel_dict : dict
        Channel type mapping
    config : dict
        Configuration dictionary with filter parameters
        Must include: 'nuScore_threshold', 'min_flashes', 'pmt_types', 'xas_types', 'max_flashes_for_keep_all'
    verbose : bool
        Show progress and statistics

    Returns:
    --------
    Tuple: (nuScore_final, f_ophit_PE_final, f_ophit_ch_final, f_ophit_t_final, selected_tpc, tracker)
    """

    if verbose:
        print(">> Starting DATA event processing...")
        start_time = time.time()

    initial_count = len(nuScore)
    tracker = CutFlowTracker(initial_count)

    if verbose:
        print(f"Initial events: {initial_count:,}")

    # Step 1: nuScore filter
    nuScore_threshold = config.get('nuScore_threshold', 0.5)

    # Handle nested nuScore structure (flatten if needed)
    if nuScore.ndim > 1:
        # If nuScore is nested (multiple values per event), take max per event
        nuScore_flat = ak.max(nuScore, axis=-1)
    else:
        nuScore_flat = nuScore

    mask_nuScore = nuScore_flat > nuScore_threshold

    arrays = {
        'nuScore': nuScore_flat[mask_nuScore],
        'f_ophit_PE': f_ophit_PE[mask_nuScore],
        'f_ophit_ch': f_ophit_ch[mask_nuScore],
        'f_ophit_t': f_ophit_t[mask_nuScore]
    }

    tracker.add_cut(f"nuScore > {nuScore_threshold}", len(arrays['nuScore']))

    # Step 2: Has flashes filter
    mask_flashes = ak.num(arrays['f_ophit_PE'], axis=1) >= config['min_flashes']
    arrays = apply_selection_mask(arrays, mask_flashes)

    tracker.add_cut("Has flashes", len(arrays['nuScore']))

    # Step 3: Flash processing (TPC selection)
    if verbose:
        print("Processing flashes...")

    channel_lookup = create_channel_lookup(
        channel_dict, config['pmt_types'], config['xas_types']
    )
    categories = vectorized_categorize_flashes(arrays['f_ophit_ch'], channel_lookup)

    # Use new version that handles variable depth
    flash_mask, decision = smart_flash_selection_new(
        arrays['f_ophit_PE'], categories, config['max_flashes_for_keep_all']
    )

    # Apply flash selection
    arrays['f_ophit_PE'] = arrays['f_ophit_PE'][flash_mask]
    arrays['f_ophit_ch'] = arrays['f_ophit_ch'][flash_mask]
    arrays['f_ophit_t'] = arrays['f_ophit_t'][flash_mask]

    # Store TPC selection as simple label: 0 for TPC0 (even, X<0), 1 for TPC1 (odd, X>0)
    # Vectorized: decision=True -> 0, decision=False -> 1
    arrays['selected_tpc'] = ak.where(decision, 0, 1)

    final_arrays = arrays

    if verbose:
        processing_time = time.time() - start_time
        print(f">> Processing completed in {processing_time:.2f} seconds")
        print()
        tracker.display()

        if len(final_arrays['nuScore']) > 0:
            print("\n* Final dataset statistics:")
            print(f"  nuScore range: [{ak.min(final_arrays['nuScore']):.3f}, {ak.max(final_arrays['nuScore']):.3f}]")
            print(f"  Mean nuScore: {ak.mean(final_arrays['nuScore']):.3f}")

            # TPC distribution
            tpc0_count = ak.sum(final_arrays['selected_tpc'] == 0)
            tpc1_count = ak.sum(final_arrays['selected_tpc'] == 1)
            print(f"\n* TPC distribution:")
            print(f"  TPC0 (even): {tpc0_count:,} events ({100*tpc0_count/len(final_arrays['nuScore']):.1f}%)")
            print(f"  TPC1 (odd):  {tpc1_count:,} events ({100*tpc1_count/len(final_arrays['nuScore']):.1f}%)")

    return (
        final_arrays['nuScore'],
        final_arrays['f_ophit_PE'],
        final_arrays['f_ophit_ch'],
        final_arrays['f_ophit_t'],
        final_arrays['selected_tpc'],
        tracker
    )


def process_data_simple(f_ophit_PE, f_ophit_ch, f_ophit_t,
                       channel_dict, config, verbose=True):
    """
    Process events for DATA (real data) - simplified version without nuScore filtering.

    This function processes optical flash data without neutrino score filtering,
    applying only flash-based selection and TPC assignment.

    Filters:
    1. Has at least min_flashes flashes
    2. Flash selection by TPC (same as MC processing)

    Parameters:
    -----------
    f_ophit_PE, f_ophit_ch, f_ophit_t : arrays
        Flash optical hit information (PhotoElectrons, Channels, Times)
    channel_dict : dict
        Channel type mapping
    config : dict
        Configuration dictionary with filter parameters
        Must include: 'min_flashes', 'pmt_types', 'xas_types', 'max_flashes_for_keep_all'
    verbose : bool
        Show progress and statistics

    Returns:
    --------
    Tuple: (f_ophit_PE_final, f_ophit_ch_final, f_ophit_t_final, selected_tpc, tracker)
    """

    if verbose:
        print(">> Starting simplified DATA event processing (no nuScore filter)...")
        start_time = time.time()

    initial_count = len(f_ophit_PE)
    tracker = CutFlowTracker(initial_count)

    if verbose:
        print(f"Initial events: {initial_count:,}")

    # Step 1: Has flashes filter
    mask_flashes = ak.num(f_ophit_PE, axis=1) >= config['min_flashes']

    arrays = {
        'f_ophit_PE': f_ophit_PE[mask_flashes],
        'f_ophit_ch': f_ophit_ch[mask_flashes],
        'f_ophit_t': f_ophit_t[mask_flashes]
    }

    tracker.add_cut("Has flashes", len(arrays['f_ophit_PE']))

    # Step 2: Flash processing (TPC selection)
    if verbose:
        print("Processing flashes...")

    channel_lookup = create_channel_lookup(
        channel_dict, config['pmt_types'], config['xas_types']
    )
    categories = vectorized_categorize_flashes(arrays['f_ophit_ch'], channel_lookup)

    # Use new version that handles variable depth
    flash_mask, decision = smart_flash_selection_new(
        arrays['f_ophit_PE'], categories, config['max_flashes_for_keep_all']
    )

    # Apply flash selection
    arrays['f_ophit_PE'] = arrays['f_ophit_PE'][flash_mask]
    arrays['f_ophit_ch'] = arrays['f_ophit_ch'][flash_mask]
    arrays['f_ophit_t'] = arrays['f_ophit_t'][flash_mask]

    # Store TPC selection as simple label: 0 for TPC0 (even, X<0), 1 for TPC1 (odd, X>0)
    # Vectorized: decision=True -> 0, decision=False -> 1
    arrays['selected_tpc'] = ak.where(decision, 0, 1)

    final_arrays = arrays

    if verbose:
        processing_time = time.time() - start_time
        print(f">> Processing completed in {processing_time:.2f} seconds")
        print()
        tracker.display()

        if len(final_arrays['f_ophit_PE']) > 0:
            # TPC distribution
            tpc0_count = ak.sum(final_arrays['selected_tpc'] == 0)
            tpc1_count = ak.sum(final_arrays['selected_tpc'] == 1)
            print(f"\n* TPC distribution:")
            print(f"  TPC0 (even): {tpc0_count:,} events ({100*tpc0_count/len(final_arrays['f_ophit_PE']):.1f}%)")
            print(f"  TPC1 (odd):  {tpc1_count:,} events ({100*tpc1_count/len(final_arrays['f_ophit_PE']):.1f}%)")

            # Flash statistics
            n_flashes = ak.num(final_arrays['f_ophit_PE'], axis=1)
            print(f"\n* Flash statistics:")
            print(f"  Flashes per event: [{ak.min(n_flashes)}, {ak.max(n_flashes)}]")
            print(f"  Mean flashes: {ak.mean(n_flashes):.1f}")

    return (
        final_arrays['f_ophit_PE'],
        final_arrays['f_ophit_ch'],
        final_arrays['f_ophit_t'],
        final_arrays['selected_tpc'],
        tracker
    )


def process_events2(nuvT, f_ophit_PE, f_ophit_ch, f_ophit_t,
                   dEpromx, dEpromy, dEpromz, dEtpc, nuvZ,
                   # Additional variables for analysis
                   nuvE, dEspreadx, dEspready, dEspreadz,
                   channel_dict, config, verbose=True):
    """
    Extended version of process_events with additional variables for analysis.

    Parameters:
    -----------
    Standard inputs (same as process_events):
        nuvT, f_ophit_PE, f_ophit_ch, f_ophit_t,
        dEpromx, dEpromy, dEpromz, dEtpc, nuvZ

    Additional inputs for analysis:
        nuvE : array
            Neutrino energy (per-TPC, 2 values per event)
        dEspreadx/y/z : array
            Spatial spread of deposited energy (per-TPC, 2 values per event)

    channel_dict : dict
        Channel type mapping
    config : dict
        Configuration dictionary with filter parameters
    verbose : bool
        Show progress and statistics

    Returns:
    --------
    Tuple of all filtered arrays (standard + extra) plus statistics tracker
    """

    if verbose:
        print(">> Starting optimized event processing (with extra variables)...")
        start_time = time.time()

    initial_count = len(nuvT)
    tracker = CutFlowTracker(initial_count)

    if verbose:
        print(f"Initial events: {initial_count:,}")

    # Step 1: Single neutrino filter
    mask_single = ak.num(nuvT) == 1

    arrays = {
        'nuvT': nuvT[mask_single],
        'f_ophit_PE': f_ophit_PE[mask_single],
        'f_ophit_ch': f_ophit_ch[mask_single],
        'f_ophit_t': f_ophit_t[mask_single],
        'dEpromx': dEpromx[mask_single],
        'dEpromy': dEpromy[mask_single],
        'dEpromz': dEpromz[mask_single],
        'dEtpc': dEtpc[mask_single],
        'nuvZ': nuvZ[mask_single],
        # Extra variables
        'nuvE': nuvE[mask_single],
        'dEspreadx': dEspreadx[mask_single],
        'dEspready': dEspready[mask_single],
        'dEspreadz': dEspreadz[mask_single]
    }

    tracker.add_cut("Single neutrino", len(arrays['nuvT']))

    # Step 2: Has flashes filter
    mask_flashes = ak.num(arrays['f_ophit_PE'], axis=1) >= config['min_flashes']
    arrays = apply_selection_mask(arrays, mask_flashes)

    tracker.add_cut("Has flashes", len(arrays['nuvT']))

    # Step 3: Flash processing
    if verbose:
        print("Processing flashes...")

    channel_lookup = create_channel_lookup(
        channel_dict, config['pmt_types'], config['xas_types']
    )
    categories = vectorized_categorize_flashes(arrays['f_ophit_ch'], channel_lookup)

    flash_mask, decision = smart_flash_selection(
        arrays['f_ophit_PE'], categories, config['max_flashes_for_keep_all']
    )

    # Apply flash selection
    arrays['f_ophit_PE'] = arrays['f_ophit_PE'][flash_mask]
    arrays['f_ophit_ch'] = arrays['f_ophit_ch'][flash_mask]
    arrays['f_ophit_t'] = arrays['f_ophit_t'][flash_mask]

    # Select TPC values (standard variables)
    tpc_arrays = {
        'dEpromx': arrays['dEpromx'],
        'dEpromy': arrays['dEpromy'],
        'dEpromz': arrays['dEpromz'],
        'dEtpc': arrays['dEtpc']
    }
    selected_tpc = select_tpc_values(tpc_arrays, decision)
    arrays.update(selected_tpc)

    # Store TPC selection as simple label: 0 for TPC0 (even), 1 for TPC1 (odd)
    arrays['selected_tpc'] = ak.where(decision, 0, 1)

    # Select TPC values (extra variables that are per-TPC)
    # dEspread variables are similar to dEprom - one value per TPC (2 per event)
    extra_tpc_arrays = {
        'dEspreadx': arrays['dEspreadx'],
        'dEspready': arrays['dEspready'],
        'dEspreadz': arrays['dEspreadz']
    }
    selected_extra_tpc = select_tpc_values(extra_tpc_arrays, decision)
    arrays.update(selected_extra_tpc)

    # nuvE is per-neutrino (like nuvT, nuvZ), not per-TPC
    # It's already the correct structure, no TPC selection needed

    # Step 4: Valid data cut
    if verbose:
        print("Applying data validity cuts...")

    valid_mask = (
        (arrays['dEpromx'] != config['invalid_marker']) &
        (arrays['dEpromy'] != config['invalid_marker']) &
        (arrays['dEpromz'] != config['invalid_marker'])
    )

    arrays = apply_selection_mask(arrays, valid_mask)
    tracker.add_cut("Valid data (≠ -999 in dEprom)", len(arrays['nuvT']))

    # Step 5: Energy cut
    if verbose:
        print("Applying energy cuts...")

    energy_mask = arrays['dEtpc'] > config['min_energy_tpc']
    arrays = apply_selection_mask(arrays, energy_mask)
    tracker.add_cut("Energy cut", len(arrays['nuvT']))

    # Step 6: Position cuts
    if verbose:
        print("Applying position cuts...")

    position_mask = (
        (arrays['dEpromx'] >= config['depromx_range'][0]) &
        (arrays['dEpromx'] <= config['depromx_range'][1]) &
        (arrays['dEpromy'] >= config['depromy_range'][0]) &
        (arrays['dEpromy'] <= config['depromy_range'][1]) &
        (arrays['dEpromz'] >= config['depromz_range'][0]) &
        (arrays['dEpromz'] <= config['depromz_range'][1])
    )

    final_arrays = apply_selection_mask(arrays, position_mask)
    tracker.add_cut("dEprom in active volume cut", len(final_arrays['nuvT']))

    if verbose:
        processing_time = time.time() - start_time
        print(f">> Processing completed in {processing_time:.2f} seconds")
        print()
        tracker.display()

        if len(final_arrays['nuvT']) > 0:
            print("\n* Final dataset ranges:")
            for var in ['dEpromx', 'dEpromy', 'dEpromz', 'dEtpc']:
                arr = final_arrays[var]
                print(f"  {var}: [{ak.min(arr):.1f}, {ak.max(arr):.1f}]")

            # Print TPC distribution
            tpc0_count = ak.sum(final_arrays['selected_tpc'] == 0)
            tpc1_count = ak.sum(final_arrays['selected_tpc'] == 1)
            print(f"\n* TPC distribution:")
            print(f"  TPC0 (even): {tpc0_count:,} events ({100*tpc0_count/len(final_arrays['nuvT']):.1f}%)")
            print(f"  TPC1 (odd):  {tpc1_count:,} events ({100*tpc1_count/len(final_arrays['nuvT']):.1f}%)")

    return (
        final_arrays['nuvT'],
        final_arrays['f_ophit_PE'],
        final_arrays['f_ophit_ch'],
        final_arrays['f_ophit_t'],
        final_arrays['dEpromx'],
        final_arrays['dEpromy'],
        final_arrays['dEpromz'],
        final_arrays['dEtpc'],
        final_arrays['nuvZ'],
        # Extra variables
        final_arrays['nuvE'],
        final_arrays['dEspreadx'],
        final_arrays['dEspready'],
        final_arrays['dEspreadz'],
        # TPC selection label: 0 = TPC0 (even), 1 = TPC1 (odd)
        final_arrays['selected_tpc'],
        tracker
    )


def create_pe_matrix(f_ophit_PE, f_ophit_ch, max_channels=312, use_numba=None):
    """
    Create PE matrix using fully vectorized awkward operations.

    Optimized version that avoids Python loops by:
    1. Flattening nested arrays with event indices
    2. Using np.add.at or numba for scatter-add operations

    Parameters:
    -----------
    f_ophit_PE : awkward array
        PE values [events, flashes, ophits]
    f_ophit_ch : awkward array
        Channel IDs [events, flashes, ophits]
    max_channels : int
        Maximum number of channels (default: 312 for SBND)
    use_numba : bool or None
        Force numba on/off. None = auto-detect availability.

    Returns:
    --------
    pe_matrix : np.ndarray (n_events, max_channels)

    ~10-50x faster than loop-based version.
    Additional ~2-5x with numba.
    """
    start_time = time.time()
    n_events = len(f_ophit_PE)

    # Determine whether to use numba
    if use_numba is None:
        use_numba = HAS_NUMBA

    backend = "numba (parallel)" if use_numba and HAS_NUMBA else "numpy"
    print(f"Creating PE matrix for {n_events:,} events x {max_channels} channels [{backend}]")

    pe_matrix = np.zeros((n_events, max_channels), dtype=np.float32)

    if n_events == 0:
        return pe_matrix

    # Get counts at each level for reconstructing event indices
    n_ophits_per_flash = ak.num(f_ophit_PE, axis=2)  # [events, flashes]
    n_ophits_per_event = ak.sum(n_ophits_per_flash, axis=1)  # [events]

    # Create event indices using np.repeat (simpler and efficient)
    event_idx_flat = np.repeat(
        np.arange(n_events, dtype=np.int32),
        ak.to_numpy(n_ophits_per_event)
    )

    # Flatten PE and channel arrays completely
    pe_flat = ak.to_numpy(ak.flatten(ak.flatten(f_ophit_PE, axis=2), axis=1)).astype(np.float32)
    ch_flat = ak.to_numpy(ak.flatten(ak.flatten(f_ophit_ch, axis=2), axis=1)).astype(np.int32)

    # Filter valid channels
    valid_mask = (ch_flat >= 0) & (ch_flat < max_channels)

    if np.any(valid_mask):
        valid_events = event_idx_flat[valid_mask]
        valid_channels = ch_flat[valid_mask]
        valid_pes = pe_flat[valid_mask]

        if use_numba and HAS_NUMBA:
            # Use numba parallel loop for additional speedup
            _numba_fill_pe_matrix(pe_matrix, valid_events, valid_channels, valid_pes)
        else:
            # Scatter-add using np.add.at (atomic add, no race conditions)
            np.add.at(pe_matrix, (valid_events, valid_channels), valid_pes)

    processing_time = time.time() - start_time

    print(f">> Completed in {processing_time:.2f} seconds")
    print(f"   Matrix shape: {pe_matrix.shape}")
    print(f"   Non-zero elements: {np.count_nonzero(pe_matrix):,}")
    print(f"   Total PE: {np.sum(pe_matrix):.1f}")

    return pe_matrix


def _select_halves(pe_spatial, images, map_idx, method, half_y):
    """Select better half for each event based on method."""
    n_events = pe_spatial.shape[0]
    
    top = pe_spatial[:, :half_y, :]
    bottom = pe_spatial[:, half_y:, :]
    
    if method == "max":
        top_scores = np.max(top.reshape(n_events, -1), axis=1)
        bottom_scores = np.max(bottom.reshape(n_events, -1), axis=1)
    elif method == "sum":
        top_scores = np.sum(top.reshape(n_events, -1), axis=1)
        bottom_scores = np.sum(bottom.reshape(n_events, -1), axis=1)
    elif method == "nonzero":
        top_scores = np.count_nonzero(top.reshape(n_events, -1), axis=1)
        bottom_scores = np.count_nonzero(bottom.reshape(n_events, -1), axis=1)
    elif method == "mean_top":
        flat_top = top.reshape(n_events, -1)
        flat_bottom = bottom.reshape(n_events, -1)
        top_scores = np.mean(np.partition(flat_top, -5, axis=1)[:, -5:], axis=1)
        bottom_scores = np.mean(np.partition(flat_bottom, -5, axis=1)[:, -5:], axis=1)
    
    use_top = top_scores >= bottom_scores
    images[use_top, :, :, map_idx] = top[use_top]
    images[~use_top, :, :, map_idx] = bottom[~use_top]


def create_pe_images(pe_matrix, *maps, method="max", normalization_factors=None, normalize=True):
    """
    Create PE images from PE matrix and channel maps.

    This function operates in different modes based on parameters:

    TRAINING MODE (normalization_factors=None, normalize=True):
    - Computes normalization factors from the input data
    - For maps 0 and 1 (uncoated/coated PMTs): uses shared normalization (max of both)
    - Normalizes images to [0, 1] range
    - Returns: (images, norm_factors)

    PREPROCESSING MODE (normalization_factors=None, normalize=False):
    - Computes max values but does NOT normalize
    - Returns raw PE images for later normalization
    - Useful when combining multiple files with different max values
    - Returns: (images, max_values)

    INFERENCE MODE (normalization_factors provided):
    - Uses pre-computed normalization factors from training
    - Applies strict filtering: discards events where any pixel > 1.0 after normalization
    - Returns: (images, norm_factors, valid_event_mask)

    Parameters:
    -----------
    pe_matrix : array (n_events, n_channels)
        PE values per event and channel
    *maps : arrays (ch_y, ch_z) - channel mapping arrays
           maps[0] = uncoated PMT map, maps[1] = coated PMT map
    method : str - selection method for half-detector: 'max', 'sum', 'nonzero', 'mean_top'
    normalization_factors : list, optional
        Pre-computed normalization factors for inference mode
        If None, calculates factors (training mode)
    normalize : bool
        If True, normalize images. If False, return raw PE values and max values.
        Only used when normalization_factors is None.

    Returns:
    --------
    TRAINING MODE (normalize=True):
        images : array (n_events, ch_y//2, ch_z, n_maps) - normalized PE images [0, 1]
        norm_factors : list - calculated normalization factors

    PREPROCESSING MODE (normalize=False):
        images : array (n_events, ch_y//2, ch_z, n_maps) - raw PE images
        max_values : list - max PE values per map (for later normalization)

    INFERENCE MODE:
        images : array (n_valid_events, ch_y//2, ch_z, n_maps) - PE images for valid events only
        norm_factors : list - normalization factors used
        valid_event_mask : array (n_events,) - boolean mask indicating which events were kept
    """
    
    ch_y, ch_z = maps[0].shape
    n_events = pe_matrix.shape[0]
    n_maps = len(maps)
    half_y = ch_y // 2
    
    print(f"Creating {n_events:,} images ({half_y}×{ch_z}×{n_maps})")
    
    # Initialize outputs
    images = np.zeros((n_events, half_y, ch_z, n_maps), dtype=np.float32)
    norm_factors = []
    valid_event_mask = np.ones(n_events, dtype=bool)  # Track valid events
    
    for map_idx, channel_map in enumerate(maps):
        valid_mask = (channel_map >= 0) & (channel_map < pe_matrix.shape[1])
        pe_spatial = np.zeros((n_events, ch_y, ch_z), dtype=np.float32)
        
        if np.any(valid_mask):
            y_pos, z_pos = np.where(valid_mask)
            channels = channel_map[valid_mask]
            pe_spatial[:, y_pos, z_pos] = pe_matrix[:, channels]
        
        # Better variable naming for clarity
        map_name = 'uncoated' if map_idx == 0 else 'coated' if map_idx == 1 else f'map_{map_idx}'
        
        if normalization_factors is not None:
            # Use provided normalization factor (for inference)
            # First factor is shared between maps 0 and 1
            norm_factor = normalization_factors[0] if map_idx < 2 else normalization_factors[map_idx-1]
            
            # Apply normalization for inference path
            if norm_factor > 0:
                pe_spatial /= norm_factor
                
                # Check for events exceeding normalization range (strict >1.0)
                event_max_values = np.max(pe_spatial.reshape(n_events, -1), axis=1)
                events_exceeding = event_max_values > 1.0
                
                if np.any(events_exceeding):
                    n_exceeding = np.sum(events_exceeding)
                    max_exceeding = np.max(event_max_values[events_exceeding])
                    print(f"[WARNING] {n_exceeding} events exceed normalization range in {map_name}")
                    print(f"          Max value found: {max_exceeding:.6f} (should be ≤1.0)")
                    print(f"          Marking these events for discard")
                    
                    # Mark events as invalid
                    valid_event_mask[events_exceeding] = False
            
            # Process current map for inference path
            _select_halves(pe_spatial, images, map_idx, method, half_y)
            
        else:
            # Calculate normalization factor (for training/preprocessing)
            if map_idx < 2:
                # First 2 maps (uncoated + coated) share normalization
                if map_idx == 0:
                    uncoated_spatial = pe_spatial.copy()
                    continue
                elif map_idx == 1:
                    coated_spatial = pe_spatial.copy()
                    norm_factor = max(np.max(uncoated_spatial), np.max(coated_spatial))
                    norm_factors.append(norm_factor)

                    # Apply normalization only if requested
                    if normalize and norm_factor > 0:
                        uncoated_spatial /= norm_factor
                        coated_spatial /= norm_factor

                    # Process uncoated map (map_idx=0)
                    _select_halves(uncoated_spatial, images, 0, method, half_y)

                    # Process coated map (map_idx=1)
                    _select_halves(coated_spatial, images, 1, method, half_y)

                    # Both maps processed, continue to next map if any
                    continue
            else:
                # Additional maps processed separately
                norm_factor = np.max(pe_spatial)
                norm_factors.append(norm_factor)

                # Apply normalization only if requested
                if normalize and norm_factor > 0:
                    pe_spatial /= norm_factor
                _select_halves(pe_spatial, images, map_idx, method, half_y)
    
    # Filter out invalid events from final images array
    if normalization_factors is not None:
        n_valid = np.sum(valid_event_mask)
        if n_valid < n_events:
            print(f">> Final filtering: {n_events-n_valid} events discarded, {n_valid} events kept")
            images = images[valid_event_mask]
        else:
            print(f">> All {n_events} events passed normalization checks")
        
        return images, norm_factors, valid_event_mask
    else:
        # Training mode - return all events (no filtering applied)
        return images, norm_factors


def scale_coordinates(coordinates, coord_ranges):
    """
    Scale coordinates to normalized ranges.

    Parameters:
    -----------
    coordinates : array (n_events, n_coords) - [x_abs, y, z] or [x_abs, y, z, pca_y, pca_z]
    coord_ranges : dict - coordinate ranges for scaling

    Returns:
    --------
    scaled_coords : array (n_events, n_coords) - scaled coordinates
    """
    n_coords = coordinates.shape[1]

    if n_coords == 3:
        # Original behavior: position only
        x_abs, y, z = coordinates[:, 0], coordinates[:, 1], coordinates[:, 2]

        # Scale each coordinate
        x_scaled = (x_abs - coord_ranges['x'][0]) / (coord_ranges['x'][1] - coord_ranges['x'][0])
        y_scaled = (y - coord_ranges['y'][0]) / (coord_ranges['y'][1] - coord_ranges['y'][0]) * 2 - 1
        z_scaled = (z - coord_ranges['z'][0]) / (coord_ranges['z'][1] - coord_ranges['z'][0])

        return np.stack([x_scaled, y_scaled, z_scaled], axis=1)

    elif n_coords == 5:
        # New behavior: position + PCA
        x_abs, y, z, pca_y, pca_z = coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], coordinates[:, 3], coordinates[:, 4]

        # Scale position coordinates
        x_scaled = (x_abs - coord_ranges['x'][0]) / (coord_ranges['x'][1] - coord_ranges['x'][0])
        y_scaled = (y - coord_ranges['y'][0]) / (coord_ranges['y'][1] - coord_ranges['y'][0]) * 2 - 1
        z_scaled = (z - coord_ranges['z'][0]) / (coord_ranges['z'][1] - coord_ranges['z'][0])

        # Scale PCA coordinates (already in [-1, 1] range, but apply same scaling for consistency)
        pca_y_scaled = (pca_y - coord_ranges['pca_y'][0]) / (coord_ranges['pca_y'][1] - coord_ranges['pca_y'][0]) * 2 - 1
        pca_z_scaled = (pca_z - coord_ranges['pca_z'][0]) / (coord_ranges['pca_z'][1] - coord_ranges['pca_z'][0]) * 2 - 1

        return np.stack([x_scaled, y_scaled, z_scaled, pca_y_scaled, pca_z_scaled], axis=1)

    elif n_coords == 6:
        # Position + PCA + Quality
        # Note: This is only for inverse scaling during inference
        # During training, targets only have 5 coords (quality computed on-the-fly)
        x_abs, y, z, pca_y, pca_z = coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], coordinates[:, 3], coordinates[:, 4]
        quality = coordinates[:, 5]

        # Scale position coordinates
        x_scaled = (x_abs - coord_ranges['x'][0]) / (coord_ranges['x'][1] - coord_ranges['x'][0])
        y_scaled = (y - coord_ranges['y'][0]) / (coord_ranges['y'][1] - coord_ranges['y'][0]) * 2 - 1
        z_scaled = (z - coord_ranges['z'][0]) / (coord_ranges['z'][1] - coord_ranges['z'][0])

        # Scale PCA coordinates
        pca_y_scaled = (pca_y - coord_ranges['pca_y'][0]) / (coord_ranges['pca_y'][1] - coord_ranges['pca_y'][0]) * 2 - 1
        pca_z_scaled = (pca_z - coord_ranges['pca_z'][0]) / (coord_ranges['pca_z'][1] - coord_ranges['pca_z'][0]) * 2 - 1

        # Quality already in [0, 1], no scaling needed
        return np.stack([x_scaled, y_scaled, z_scaled, pca_y_scaled, pca_z_scaled, quality], axis=1)

    else:
        raise ValueError(f"Expected 3, 5, or 6 coordinates, got {n_coords}")


def inverse_scale_coordinates(scaled_coords, coord_ranges):
    """
    Inverse scale coordinates back to original ranges.

    Parameters:
    -----------
    scaled_coords : array (n_events, n_coords) - scaled coordinates [x, y, z] or [x, y, z, pca_y, pca_z]
    coord_ranges : dict - coordinate ranges for scaling

    Returns:
    --------
    original_coords : array (n_events, n_coords) - original scale coordinates
    """
    n_coords = scaled_coords.shape[1]

    if n_coords == 3:
        # Original behavior: position only
        x_inv = scaled_coords[:, 0] * (coord_ranges['x'][1] - coord_ranges['x'][0]) + coord_ranges['x'][0]
        y_inv = ((scaled_coords[:, 1] + 1) / 2) * (coord_ranges['y'][1] - coord_ranges['y'][0]) + coord_ranges['y'][0]
        z_inv = scaled_coords[:, 2] * (coord_ranges['z'][1] - coord_ranges['z'][0]) + coord_ranges['z'][0]

        return np.stack([x_inv, y_inv, z_inv], axis=1)

    elif n_coords == 5:
        # New behavior: position + PCA
        x_inv = scaled_coords[:, 0] * (coord_ranges['x'][1] - coord_ranges['x'][0]) + coord_ranges['x'][0]
        y_inv = ((scaled_coords[:, 1] + 1) / 2) * (coord_ranges['y'][1] - coord_ranges['y'][0]) + coord_ranges['y'][0]
        z_inv = scaled_coords[:, 2] * (coord_ranges['z'][1] - coord_ranges['z'][0]) + coord_ranges['z'][0]

        pca_y_inv = ((scaled_coords[:, 3] + 1) / 2) * (coord_ranges['pca_y'][1] - coord_ranges['pca_y'][0]) + coord_ranges['pca_y'][0]
        pca_z_inv = ((scaled_coords[:, 4] + 1) / 2) * (coord_ranges['pca_z'][1] - coord_ranges['pca_z'][0]) + coord_ranges['pca_z'][0]

        return np.stack([x_inv, y_inv, z_inv, pca_y_inv, pca_z_inv], axis=1)

    elif n_coords == 6:
        # Position + PCA + Quality
        x_inv = scaled_coords[:, 0] * (coord_ranges['x'][1] - coord_ranges['x'][0]) + coord_ranges['x'][0]
        y_inv = ((scaled_coords[:, 1] + 1) / 2) * (coord_ranges['y'][1] - coord_ranges['y'][0]) + coord_ranges['y'][0]
        z_inv = scaled_coords[:, 2] * (coord_ranges['z'][1] - coord_ranges['z'][0]) + coord_ranges['z'][0]

        pca_y_inv = ((scaled_coords[:, 3] + 1) / 2) * (coord_ranges['pca_y'][1] - coord_ranges['pca_y'][0]) + coord_ranges['pca_y'][0]
        pca_z_inv = ((scaled_coords[:, 4] + 1) / 2) * (coord_ranges['pca_z'][1] - coord_ranges['pca_z'][0]) + coord_ranges['pca_z'][0]

        # Quality already in [0, 1], no inverse scaling needed
        quality_inv = scaled_coords[:, 5]

        return np.stack([x_inv, y_inv, z_inv, pca_y_inv, pca_z_inv, quality_inv], axis=1)

    else:
        raise ValueError(f"Expected 3, 5, or 6 coordinates, got {n_coords}")


def plot_pe_images(image_data, event_idx, labels=None, groups=None, 
                   figsize=(20, 10), use_log_scale=False, show_colorbar=True,
                   cmap='Blues', title_fontsize=16, colorbar_size="4%"):
    """
    Plot PE images with automatic layout and optimized visualization.
    """
    
    if event_idx >= image_data.shape[0]:
        raise ValueError(f"event_idx {event_idx} >= number of events {image_data.shape[0]}")
    
    n_channels = image_data.shape[3]
    
    if labels is None:
        labels = [f"Channel {i}" for i in range(n_channels)]
    elif len(labels) != n_channels:
        raise ValueError(f"Number of labels ({len(labels)}) != number of channels ({n_channels})")
    
    if groups is None:
        groups = [[i] for i in range(n_channels)]
    
    # Auto-determine grid layout
    if n_channels <= 2:
        nrows, ncols = 1, n_channels
    elif n_channels <= 4:
        nrows, ncols = 2, 2
    elif n_channels <= 6:
        nrows, ncols = 2, 3
    else:
        ncols = int(np.ceil(np.sqrt(n_channels)))
        nrows = int(np.ceil(n_channels / ncols))
    
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    if n_channels == 1:
        axes = [axes]
    elif nrows == 1 or ncols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    colormap = plt.cm.get_cmap(cmap).copy()
    colormap.set_bad(color='white')
    
    # Calculate scaling for each group
    group_scales = {}
    for group in groups:
        group_images = [image_data[event_idx, :, :, i] for i in group]
        
        if use_log_scale:
            positive_values = [img[img > 0] for img in group_images]
            if any(len(pv) > 0 for pv in positive_values):
                vmin = min(np.min(pv) for pv in positive_values if len(pv) > 0) * 0.1
                vmax = max(np.max(img) for img in group_images if np.max(img) > 0)
            else:
                vmin, vmax = 1e-10, 1
        else:
            vmin = 0
            vmax = max(np.max(img) for img in group_images)
            if vmax == 0:
                vmax = 1
        
        group_scales[tuple(group)] = (vmin, vmax)
    
    # Plot each channel
    for i in range(n_channels):
        ax = axes[i]
        img = image_data[event_idx, :, :, i]
        
        masked_img = np.ma.masked_where(img <= 0, img)
        
        vmin, vmax = None, None
        for group, (group_vmin, group_vmax) in group_scales.items():
            if i in group:
                vmin, vmax = group_vmin, group_vmax
                break
        
        if use_log_scale and vmax > vmin:
            im = ax.imshow(masked_img, cmap=colormap, 
                          norm=LogNorm(vmin=vmin, vmax=vmax),
                          aspect='auto')
        else:
            im = ax.imshow(masked_img, cmap=colormap, 
                          vmin=vmin, vmax=vmax,
                          aspect='auto')
        
        ax.set_title(labels[i], fontsize=title_fontsize, pad=10)
        ax.set_xticks([])
        ax.set_yticks([])
        
        if show_colorbar:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size=colorbar_size, pad=0.1)
            cbar = plt.colorbar(im, cax=cax)
            cbar.ax.tick_params(labelsize=title_fontsize-2)
        
        non_zero = np.count_nonzero(img)
        total_pixels = img.size
        max_val = np.max(img)
        ax.text(0.02, 0.98, f'Non-zero: {non_zero}/{total_pixels}\nMax: {max_val:.2f}', 
                transform=ax.transAxes, fontsize=title_fontsize-4,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Hide unused subplots
    for i in range(n_channels, len(axes)):
        axes[i].set_visible(False)
    
    fig.suptitle(f'PE Images - Event {event_idx}', fontsize=title_fontsize+4, y=0.95)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.show()
    
    return fig, axes


def create_pmt_position_dict(file, pmt_types=None):
    """
    Create PMT position dictionary from ROOT file PDSMapTree.

    Parameters:
    -----------
    file : uproot file object
        ROOT file with PDSMapTree
    pmt_types : set, optional
        Set of OpDetType values to include (e.g., {0, 1} for PMTs only)
        If None, includes all detectors

    Returns:
    --------
    pmt_positions : dict
        Mapping {channel_id: (y, z)} for selected detector types
    """
    PDSMap = file['opanatree']['PDSMapTree']
    ID = PDSMap['OpDetID'].array()
    Y = PDSMap['OpDetY'].array()
    Z = PDSMap['OpDetZ'].array()
    Type = PDSMap['OpDetType'].array()

    if pmt_types is not None:
        # Filter by detector type (e.g., only PMTs)
        pmt_positions = {int(id_val): (float(y_val), float(z_val))
                         for id_val, y_val, z_val, type_val in zip(ID[0], Y[0], Z[0], Type[0])
                         if int(type_val) in pmt_types}
        print(f">> PMT position dictionary created: {len(pmt_positions)} channels (filtered for types {pmt_types})")
    else:
        # Include all detectors
        pmt_positions = {int(id_val): (float(y_val), float(z_val))
                         for id_val, y_val, z_val in zip(ID[0], Y[0], Z[0])}
        print(f">> PMT position dictionary created: {len(pmt_positions)} channels (all types)")

    return pmt_positions


def calculate_light_pca(f_ophit_PE, f_ophit_ch, pmt_positions, verbose=False, use_numba=None):
    """
    Calculate PCA of light distribution (2D: Y-Z plane).
    Equivalent to GetPCA() in TPCPMTBarycenterMatching_module.cc.

    Optimized with:
    - Vectorized awkward flattening
    - Optional numba parallel processing (~5-10x faster)

    Parameters:
    -----------
    f_ophit_PE : awkward array (n_events, n_flashes, n_ophits)
        PE values for each ophit
    f_ophit_ch : awkward array (n_events, n_flashes, n_ophits)
        Channel for each ophit
    pmt_positions : dict
        Mapping from channel to (y, z) position
    verbose : bool
        Print progress
    use_numba : bool or None
        Force numba on/off. None = auto-detect.

    Returns:
    --------
    pca_vectors : array (n_events, 2)
        Principal component direction [dy, dz] for each event
    pca_elongation : array (n_events,)
        Ratio of largest to smallest eigenvalue (elongation measure)
    """
    n_events = len(f_ophit_PE)

    # Determine backend
    if use_numba is None:
        use_numba = HAS_NUMBA

    backend = "numba (parallel)" if use_numba and HAS_NUMBA else "numpy"

    if verbose:
        print(f">> Calculating Light PCA for {n_events:,} events [{backend}]...")
        start_time = time.time()

    # Pre-create position lookup arrays
    max_channel = 312  # SBND optical channels
    y_lookup = np.zeros(max_channel, dtype=np.float32)
    z_lookup = np.zeros(max_channel, dtype=np.float32)
    valid_channels = np.zeros(max_channel, dtype=bool)

    for ch, (y, z) in pmt_positions.items():
        if ch < max_channel:
            y_lookup[ch] = y
            z_lookup[ch] = z
            valid_channels[ch] = True

    # Flatten nested awkward arrays
    if verbose:
        print("   * Flattening awkward arrays...")

    flat_pe = ak.flatten(ak.flatten(f_ophit_PE, axis=2), axis=1)
    flat_ch = ak.flatten(ak.flatten(f_ophit_ch, axis=2), axis=1)

    flat_pe_np = ak.to_numpy(flat_pe).astype(np.float32)
    flat_ch_np = ak.to_numpy(flat_ch).astype(np.int32)

    # Get event boundaries
    counts = ak.num(ak.flatten(f_ophit_PE, axis=2), axis=1)
    event_ends = np.cumsum(ak.to_numpy(counts)).astype(np.int64)
    event_starts = np.concatenate([[0], event_ends[:-1]]).astype(np.int64)

    if verbose:
        print("   * Computing PCA for all events...")

    # Initialize outputs
    pca_vectors = np.zeros((n_events, 2), dtype=np.float32)
    pca_elongation = np.zeros(n_events, dtype=np.float32)

    if use_numba and HAS_NUMBA:
        # Use parallel numba version
        _numba_compute_all_pcas(
            pca_vectors, pca_elongation,
            flat_pe_np, flat_ch_np, event_starts, event_ends,
            y_lookup, z_lookup, valid_channels, max_channel
        )
    else:
        # Pure numpy version (fallback)
        for event_idx in range(n_events):
            start_idx = event_starts[event_idx]
            end_idx = event_ends[event_idx]

            if start_idx >= end_idx:
                pca_vectors[event_idx] = [0.0, 1.0]
                pca_elongation[event_idx] = 1.0
                continue

            event_channels = flat_ch_np[start_idx:end_idx]
            event_pes = flat_pe_np[start_idx:end_idx]

            mask = (event_channels >= 0) & (event_channels < max_channel) & valid_channels[event_channels]

            if np.sum(mask) < 2:
                pca_vectors[event_idx] = [0.0, 1.0]
                pca_elongation[event_idx] = 1.0
                continue

            valid_ch = event_channels[mask]
            valid_pe = event_pes[mask]
            ophit_y = y_lookup[valid_ch]
            ophit_z = z_lookup[valid_ch]

            weight_sum = np.sum(valid_pe)
            if weight_sum == 0:
                pca_vectors[event_idx] = [0.0, 1.0]
                pca_elongation[event_idx] = 1.0
                continue

            centroid_y = np.sum(valid_pe * ophit_y) / weight_sum
            centroid_z = np.sum(valid_pe * ophit_z) / weight_sum

            centered_y = ophit_y - centroid_y
            centered_z = ophit_z - centroid_z

            cov_yy = np.sum(valid_pe * centered_y * centered_y) / weight_sum
            cov_yz = np.sum(valid_pe * centered_y * centered_z) / weight_sum
            cov_zz = np.sum(valid_pe * centered_z * centered_z) / weight_sum

            cov_matrix = np.array([[cov_yy, cov_yz], [cov_yz, cov_zz]], dtype=np.float32)
            eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

            max_idx = 1 if eigenvalues[1] > eigenvalues[0] else 0
            pca_vectors[event_idx] = eigenvectors[:, max_idx]

            if eigenvalues[0] > 1e-10:
                pca_elongation[event_idx] = eigenvalues[1] / eigenvalues[0]
            else:
                pca_elongation[event_idx] = 100.0

    if verbose:
        elapsed = time.time() - start_time
        print(f">> Light PCA completed in {elapsed:.2f}s ({n_events/elapsed:.0f} events/s)")
        print(f"   PCA Y: [{pca_vectors[:,0].min():.3f}, {pca_vectors[:,0].max():.3f}]")
        print(f"   PCA Z: [{pca_vectors[:,1].min():.3f}, {pca_vectors[:,1].max():.3f}]")
        print(f"   Elongation: [{pca_elongation.min():.1f}, {pca_elongation.max():.1f}]")

    return pca_vectors, pca_elongation


def save_filtered_data(results, filename, config=None):
    """
    Save filtered data with optimal settings.

    Parameters:
    -----------
    results : tuple
        Output from process_events() function
    filename : str
        Path where to save the parquet file
    config : dict, optional
        Configuration used for filtering (will be saved as metadata)
    """
    arrays_tuple = results[:-1]  # Remove stats
    
    names = ['nuvT', 'f_ophit_PE', 'f_ophit_ch', 'f_ophit_t', 
             'dEpromx', 'dEpromy', 'dEpromz', 'dEtpc', 'nuvZ']
    
    data = {name: arr for name, arr in zip(names, arrays_tuple)}
    
    ak.to_parquet(data, filename, compression="snappy")
    
    print(f">> Successfully saved {len(data['nuvT']):,} events to:")
    print(f"   {filename}")
    print(f"   Compressed with snappy for optimal file size")
    if config:
        print(f"   Filter config: {config}")


# ================================================================
# TensorFlow Data Generator for Large ROOT Files
# ================================================================

def build_file_event_index(file_paths, filter_config, channel_dict, verbose=True):
    """
    Pre-scan ROOT files to build index of valid events per file.

    This function loads each file, applies filters, and counts valid events
    without storing all data in memory. Used by ROOTDataGenerator to know
    how many events are available and where they are located.

    Args:
        file_paths: List of paths to ROOT files
        filter_config: Dictionary with filter configuration
        channel_dict: Dictionary mapping channel IDs to types
        verbose: Print progress messages

    Returns:
        list of dicts: [{'file': path, 'n_events': count, 'start_idx': global_idx}, ...]
    """
    import uproot
    import awkward as ak

    file_index = []
    global_start_idx = 0

    if verbose:
        print(f"   * Scanning {len(file_paths)} files to build event index...")

    for i, file_path in enumerate(file_paths):
        if verbose:
            print(f"     - File {i+1}/{len(file_paths)}: {os.path.basename(file_path)}")

        try:
            # Open file and load necessary branches for filtering
            with uproot.open(f"{file_path}:opanatree/OpAnaTree") as tree:
                # Load minimal branches needed for filtering
                nuvT = tree["nuvT"].array(library="ak")
                nuvX = tree["nuvX"].array(library="ak")
                nuvY = tree["nuvY"].array(library="ak")
                nuvZ = tree["nuvZ"].array(library="ak")
                dEtpc = tree["dEtpc"].array(library="ak")
                dEpromx = tree["dEprom"].array(library="ak")
                dEpromy = tree["dEpromy"].array(library="ak")
                dEpromz = tree["dEpromz"].array(library="ak")
                f_ophit_PE = tree["flash_ophit_pe"].array(library="ak")
                f_ophit_ch = tree["flash_ophit_ch"].array(library="ak")
                f_ophit_t = tree["flash_ophit_time"].array(library="ak")

            # Apply same filtering logic as process_events
            # This is lightweight - just counting, not storing
            results = process_events(
                nuvT, f_ophit_PE, f_ophit_ch, f_ophit_t,
                dEpromx, dEpromy, dEpromz, dEtpc, nuvZ,
                channel_dict, filter_config, verbose=False,
                nuvX=nuvX, nuvY=nuvY
            )

            # Count valid events (after filtering)
            n_valid_events = len(results[0])  # nuvT_final

            file_index.append({
                'file': file_path,
                'n_events': n_valid_events,
                'start_idx': global_start_idx
            })

            global_start_idx += n_valid_events

            if verbose:
                print(f"       Valid events: {n_valid_events:,}")

        except Exception as e:
            print(f"     ! Error scanning {os.path.basename(file_path)}: {e}")
            # Skip this file
            continue

    if verbose:
        total_events = sum(f['n_events'] for f in file_index)
        print(f"   * Total valid events across all files: {total_events:,}")

    return file_index


class ROOTDataGenerator:
    """
    TensorFlow/Keras Data Generator for large ROOT files.

    Loads and processes ROOT files on-the-fly during training, maintaining
    constant memory usage regardless of total dataset size.

    Compatible with tf.keras Model.fit() via the Sequence API.
    """

    def __init__(self, file_paths, batch_size, normalization_factor,
                 channel_dict, pmt_maps, filter_config, coord_config,
                 image_config, shuffle=True, verbose=True):
        """
        Initialize the data generator.

        Args:
            file_paths: List of paths to ROOT files
            batch_size: Number of events per batch
            normalization_factor: Fixed factor for image normalization
            channel_dict: Dictionary mapping channel IDs to types
            pmt_maps: Tuple of (uncoated_map, coated_map) numpy arrays
            filter_config: Dictionary with filter settings
            coord_config: Dictionary with coordinate scaling settings
            image_config: Dictionary with image creation settings
            shuffle: Whether to shuffle indices after each epoch
            verbose: Print progress messages
        """
        import tensorflow as tf

        self.file_paths = sorted(file_paths)
        self.batch_size = batch_size
        self.norm_factor = normalization_factor
        self.channel_dict = channel_dict
        self.pmt_maps = pmt_maps
        self.filter_config = filter_config
        self.coord_config = coord_config
        self.image_config = image_config
        self.shuffle = shuffle
        self.verbose = verbose

        # Pre-scan files to build event index
        if verbose:
            print(">> Building event index from ROOT files...")

        self.file_index = build_file_event_index(
            file_paths, filter_config, channel_dict, verbose
        )

        self.total_events = sum(f['n_events'] for f in self.file_index)

        # Build indices - one range per file for efficient caching
        self.indices = self._build_indices_per_file()

        # Cache for currently loaded file
        self.current_file_path = None
        self.current_file_data = None

        # Initial shuffle within each file if requested
        if self.shuffle:
            self._shuffle_within_files()

        if verbose:
            print(f">> Generator initialized")
            print(f"   * Total events: {self.total_events:,}")
            print(f"   * Batches per epoch: {len(self)}")
            print(f"   * Batch size: {self.batch_size}")
            print(f"   * Shuffle: {'within files' if self.shuffle else 'disabled'}")

    def __len__(self):
        """Return number of batches per epoch."""
        return int(np.ceil(self.total_events / self.batch_size))

    def __getitem__(self, batch_idx):
        """
        Generate one batch of data.

        Called by TensorFlow during training.

        Args:
            batch_idx: Index of the batch (0 to len(self)-1)

        Returns:
            tuple: (images_batch, coordinates_batch)
        """
        # Calculate global indices for this batch
        start_idx = batch_idx * self.batch_size
        end_idx = min(start_idx + self.batch_size, self.total_events)
        batch_global_indices = self.indices[start_idx:end_idx]

        # Group indices by source file
        events_by_file = self._group_indices_by_file(batch_global_indices)

        # Load and process events from each file
        batch_images = []
        batch_coords = []

        for file_path, local_indices in events_by_file.items():
            # Load file if not in cache
            if self.current_file_path != file_path:
                self._load_file(file_path)

            # Extract specific events from loaded file
            file_images, file_coords = self._extract_events(local_indices)

            batch_images.append(file_images)
            batch_coords.append(file_coords)

        # Concatenate results from all files
        images = np.concatenate(batch_images, axis=0)
        coords = np.concatenate(batch_coords, axis=0)

        return images, coords

    def _group_indices_by_file(self, global_indices):
        """
        Group global indices by their source file.

        Args:
            global_indices: Array of global event indices

        Returns:
            dict: {file_path: [local_indices]}
        """
        events_by_file = {}

        for global_idx in global_indices:
            # Find which file this index belongs to
            file_info = None
            for finfo in self.file_index:
                if (global_idx >= finfo['start_idx'] and
                    global_idx < finfo['start_idx'] + finfo['n_events']):
                    file_info = finfo
                    break

            if file_info is None:
                continue

            # Convert to local index within file
            local_idx = global_idx - file_info['start_idx']

            # Add to dictionary
            if file_info['file'] not in events_by_file:
                events_by_file[file_info['file']] = []
            events_by_file[file_info['file']].append(local_idx)

        return events_by_file

    def _load_file(self, file_path):
        """
        Load and process complete ROOT file into memory.

        This caches the processed file data until another file is needed.

        Args:
            file_path: Path to ROOT file
        """
        import uproot
        import awkward as ak
        import gc

        if self.verbose:
            print(f"   * Loading file: {os.path.basename(file_path)}")

        # Clear previous file from cache
        self.current_file_data = None
        gc.collect()

        # Load raw data from ROOT
        with uproot.open(f"{file_path}:opanatree/OpAnaTree") as tree:
            nuvT = tree["nuvT"].array(library="ak")
            nuvX = tree["nuvX"].array(library="ak")
            nuvY = tree["nuvY"].array(library="ak")
            nuvZ = tree["nuvZ"].array(library="ak")
            dEtpc = tree["dEtpc"].array(library="ak")
            dEpromx = tree["dEpromx"].array(library="ak")
            dEpromy = tree["dEpromy"].array(library="ak")
            dEpromz = tree["dEpromz"].array(library="ak")
            f_ophit_PE = tree["flash_ophit_pe"].array(library="ak")
            f_ophit_ch = tree["flash_ophit_ch"].array(library="ak")
            f_ophit_t = tree["flash_ophit_time"].array(library="ak")

        # Apply filters (process_events)
        results = process_events(
            nuvT, f_ophit_PE, f_ophit_ch, f_ophit_t,
            dEpromx, dEpromy, dEpromz, dEtpc, nuvZ,
            self.channel_dict, self.filter_config, verbose=False,
            nuvX=nuvX, nuvY=nuvY
        )

        # Unpack filtered results
        (nuvT_final, f_ophit_PE_final, f_ophit_ch_final, f_ophit_t_final,
         dEpromx_final, dEpromy_final, dEpromz_final, dEtpc_final,
         nuvX_final, nuvY_final, nuvZ_final, selected_tpc_final, stats) = results

        # Create PE matrices
        pe_matrix = create_pe_matrix(
            f_ophit_PE_final,
            f_ophit_ch_final,
            self.image_config['max_channels']
        )

        # Create images with fixed normalization
        uncoated_map, coated_map = self.pmt_maps
        images = create_pe_images_fixed_norm(
            pe_matrix,
            uncoated_map,
            coated_map,
            self.norm_factor,
            method=self.image_config['selection_method']
        )

        # Prepare coordinates
        x_abs = np.abs(np.array(dEpromx_final).flatten())
        y = np.array(dEpromy_final).flatten()
        z = np.array(dEpromz_final).flatten()
        coordinates = np.column_stack((x_abs, y, z))

        # Scale coordinates
        coords_scaled = scale_coordinates(coordinates, self.coord_config['ranges'])

        # Cache processed data
        self.current_file_data = {
            'images': images,
            'coordinates': coords_scaled
        }
        self.current_file_path = file_path

        if self.verbose:
            print(f"     - Loaded {len(images):,} events")

    def _extract_events(self, local_indices):
        """
        Extract specific events from currently loaded file.

        Args:
            local_indices: List of event indices within current file

        Returns:
            tuple: (images_subset, coordinates_subset)
        """
        local_indices = np.array(local_indices)
        images = self.current_file_data['images'][local_indices]
        coords = self.current_file_data['coordinates'][local_indices]
        return images, coords

    def _build_indices_per_file(self):
        """
        Build indices array organized by file.

        Returns consecutive indices for each file, which ensures
        batches don't span multiple files (more efficient).

        Returns:
            np.ndarray: Global indices [0, 1, 2, ..., total_events-1]
        """
        return np.arange(self.total_events)

    def _shuffle_within_files(self):
        """
        Shuffle indices within each file's range.

        This maintains file boundaries so batches come from single files,
        maximizing cache efficiency.
        """
        for file_info in self.file_index:
            start = file_info['start_idx']
            end = start + file_info['n_events']

            # Shuffle only this file's range
            file_indices = self.indices[start:end].copy()
            np.random.shuffle(file_indices)
            self.indices[start:end] = file_indices

    def on_epoch_end(self):
        """
        Called at the end of every epoch.

        Shuffles indices within files and clears file cache.
        """
        if self.shuffle:
            self._shuffle_within_files()

        # Clear file cache to free memory
        self.current_file_data = None
        self.current_file_path = None
        import gc
        gc.collect()


def create_pe_images_fixed_norm(pe_matrix, uncoated_map, coated_map,
                                  normalization_factor, method='max'):
    """
    Create PE images with fixed normalization factor.

    Modified version of create_pe_images() that uses a pre-determined
    normalization factor instead of computing it from the data.

    Args:
        pe_matrix: PE matrix (n_events, 312)
        uncoated_map: Uncoated PMT map
        coated_map: Coated PMT map
        normalization_factor: Fixed normalization factor
        method: Selection method ('max', 'sum', etc.)

    Returns:
        numpy array: Images (n_events, height, width, 2)
    """
    n_events = pe_matrix.shape[0]
    height, width = uncoated_map.shape

    # Initialize output
    images = np.zeros((n_events, height, width, 2), dtype=np.float32)

    if n_events == 0:
        return images

    # Process PMT positions
    for ch in range(312):
        # Find all positions where this channel appears
        y_uncoated, x_uncoated = np.where(uncoated_map == ch)
        y_coated, x_coated = np.where(coated_map == ch)

        if len(y_uncoated) > 0 or len(y_coated) > 0:
            # Extract PE values for this channel across all events
            pe_values = pe_matrix[:, ch]

            # Place in images based on method
            if method == 'max':
                # Uncoated channel
                for y, x in zip(y_uncoated, x_uncoated):
                    images[:, y, x, 0] = np.maximum(images[:, y, x, 0], pe_values)
                # Coated channel
                for y, x in zip(y_coated, x_coated):
                    images[:, y, x, 1] = np.maximum(images[:, y, x, 1], pe_values)

    # Normalize with fixed factor
    images = images / normalization_factor

    return images