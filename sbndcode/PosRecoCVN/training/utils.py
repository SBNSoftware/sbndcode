"""
Utility functions for PosRecoCVN training pipeline.
Contains data processing, filtering, and visualization functions.
"""

import awkward as ak
import numpy as np
import time
import pandas as pd
from IPython.display import display
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
    """Fast flash categorization for irregular arrays."""
    first_channels = f_ophit_ch[:, :, 0]
    
    if not channel_lookup:
        return ak.full_like(first_channels, -1)
    
    def categorize_channel(ch):
        if ch in channel_lookup:
            return channel_lookup[ch]
        return -1
    
    categories = ak.Array([
        [categorize_channel(ch) for ch in event_channels]
        for event_channels in ak.to_list(first_channels)
    ])
    
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
    Select TPC values based on decision - simple version using loops.

    Arrays have structure [n_events, 2] where index 0=TPC0, index 1=TPC1
    Decision: True=TPC0 (even), False=TPC1 (odd)
    """
    # Convert decision to list
    decision_list = ak.to_list(decision)

    result = {}
    for name, array in tpc_dict.items():
        array_list = ak.to_list(array)

        selected_values = []
        for event_idx, dec in enumerate(decision_list):
            event_values = array_list[event_idx]

            if event_values is None or len(event_values) < 2:
                selected_values.append(0.0)  # Default
            else:
                # dec=True -> select index 0 (TPC0), dec=False -> select index 1 (TPC1)
                idx = 0 if dec else 1
                selected_values.append(float(event_values[idx]))

        result[name] = ak.Array(selected_values)

    return result


def process_events(nuvT, f_ophit_PE, f_ophit_ch, f_ophit_t,
                  dEpromx, dEpromy, dEpromz, dEtpc, nuvZ,
                  channel_dict, config, verbose=True):
    """
    Process events with optimized pipeline.

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
        print(">> Starting optimized event processing...")
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
        'nuvZ': nuvZ[mask_single]
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

    # Select TPC values (only for dEprom and dEtpc - these have TPC structure)
    tpc_arrays = {
        'dEpromx': arrays['dEpromx'],
        'dEpromy': arrays['dEpromy'],
        'dEpromz': arrays['dEpromz'],
        'dEtpc': arrays['dEtpc']
    }
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
    tracker.add_cut("Position cut", len(final_arrays['nuvT']))
    
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
        final_arrays['nuvZ'],
        final_arrays['selected_tpc'],
        tracker
    )


def smart_flash_selection_new(f_ophit_PE, categories, max_flashes_for_keep_all=2):
    """
    Simple flash selection - uses basic Python loops to avoid UnionArray issues.

    Logic:
    1. Sum PEs per flash
    2. Sum flashes by TPC (even channels = TPC0, odd channels = TPC1)
    3. Pick TPC with more total PE
    """
    # Sum PE per flash (axis=2 sums ophits within each flash)
    sum_pe = ak.sum(f_ophit_PE, axis=2)

    # Convert everything to lists for simple processing
    sum_pe_list = ak.to_list(sum_pe)
    categories_list = ak.to_list(categories)

    n_events = len(sum_pe_list)
    decisions = []
    selection_masks = []

    for event_idx in range(n_events):
        event_sum_pe = sum_pe_list[event_idx]
        event_categories = categories_list[event_idx]

        if event_sum_pe is None or len(event_sum_pe) == 0:
            decisions.append(True)  # Default to TPC0 (even)
            selection_masks.append([])
            continue

        n_flashes = len(event_sum_pe)

        # Sum PEs by TPC
        sum_even = 0.0
        sum_odd = 0.0

        for flash_idx in range(n_flashes):
            pe = float(event_sum_pe[flash_idx]) if event_sum_pe[flash_idx] is not None else 0.0
            cat = event_categories[flash_idx]

            # Categories: 0,2 = even (TPC0), 1,3 = odd (TPC1)
            if cat in [0, 2]:
                sum_even += pe
            elif cat in [1, 3]:
                sum_odd += pe

        # Decide which TPC to use
        decision = sum_even >= sum_odd  # True = TPC0 (even), False = TPC1 (odd)
        decisions.append(decision)

        # Create selection mask
        if n_flashes <= max_flashes_for_keep_all:
            # Keep all flashes
            mask = [True] * n_flashes
        else:
            # Keep only flashes from selected TPC
            if decision:  # Keep even (TPC0)
                mask = [cat in [0, 2] for cat in event_categories]
            else:  # Keep odd (TPC1)
                mask = [cat in [1, 3] for cat in event_categories]

        selection_masks.append(mask)

    # Convert back to awkward arrays
    selection_mask = ak.Array(selection_masks)
    decision = ak.Array(decisions)

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
    decision_list = ak.to_list(decision)
    arrays['selected_tpc'] = ak.Array([0 if d else 1 for d in decision_list])


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
    tracker.add_cut("Position cut", len(final_arrays['nuvT']))

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
    decision_list = ak.to_list(decision)
    arrays['selected_tpc'] = ak.Array([0 if d else 1 for d in decision_list])

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
    tracker.add_cut("Position cut", len(final_arrays['nuvT']))

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


def create_pe_matrix(f_ophit_PE, f_ophit_ch, max_channels=312):
    """Create PE matrix using optimized approach for irregular arrays."""

    start_time = time.time()
    n_events = len(f_ophit_PE)

    print(f"Creating PE matrix for {n_events:,} events x {max_channels} channels")

    pe_matrix = np.zeros((n_events, max_channels), dtype=np.float32)
    
    for i in range(n_events):
        if i % 10000 == 0:
            print(f"  Processing event {i:,}/{n_events:,} ({100*i/n_events:.1f}%)")
        
        for j in range(len(f_ophit_PE[i])):
            pe_flash = ak.to_numpy(f_ophit_PE[i][j])
            ch_flash = ak.to_numpy(f_ophit_ch[i][j])
            
            valid_mask = (ch_flash >= 0) & (ch_flash < max_channels)
            if np.any(valid_mask):
                valid_channels = ch_flash[valid_mask]
                valid_pes = pe_flash[valid_mask]
                
                binned = np.bincount(valid_channels, weights=valid_pes, minlength=max_channels)
                pe_matrix[i] += binned[:max_channels]
    
    processing_time = time.time() - start_time
    
    print(f">> Completed in {processing_time:.2f} seconds")
    print(f"Matrix shape: {pe_matrix.shape}")
    print(f"Non-zero elements: {np.count_nonzero(pe_matrix):,}")
    print(f"Total PE: {np.sum(pe_matrix):.1f}")
    
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


def create_pe_images(pe_matrix, *maps, method="max", normalization_factors=None):
    """
    Create PE images from PE matrix and channel maps.
    
    This function operates in two distinct modes based on the normalization_factors parameter:
    
    TRAINING MODE (normalization_factors=None):
    - Computes normalization factors from the input data
    - For maps 0 and 1 (uncoated/coated PMTs): uses shared normalization (max of both)
    - For additional maps: uses individual normalization per map
    - Processes all events without filtering
    - Returns: (images, norm_factors)
    
    INFERENCE MODE (normalization_factors provided):
    - Uses pre-computed normalization factors from training
    - Applies strict filtering: discards events where any pixel > 1.0 after normalization
    - This prevents the model from seeing out-of-distribution data (values > training range)
    - Events exceeding the normalization range are marked as invalid
    - Returns: (images, norm_factors, valid_event_mask)
    
    Image Processing Details:
    - Converts PE matrix to spatial images using channel maps
    - Applies half-detector selection based on specified method
    - Maps are processed in order: uncoated PMTs, coated PMTs, additional detectors
    - Final images have shape (n_valid_events, ch_y//2, ch_z, n_maps)
    
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
    
    Returns:
    --------
    TRAINING MODE:
        images : array (n_events, ch_y//2, ch_z, n_maps) - PE images
        norm_factors : list - calculated normalization factors
    
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
            # Calculate normalization factor (for training)
            if map_idx < 2:
                # First 2 maps (uncoated + coated) share normalization
                if map_idx == 0:
                    uncoated_spatial = pe_spatial.copy()
                    continue
                elif map_idx == 1:
                    coated_spatial = pe_spatial.copy()
                    norm_factor = max(np.max(uncoated_spatial), np.max(coated_spatial))
                    norm_factors.append(norm_factor)
                    
                    # Apply normalization and process both maps
                    if norm_factor > 0:
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
                
                # Apply normalization and process current map
                if norm_factor > 0:
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


def calculate_light_pca(f_ophit_PE, f_ophit_ch, pmt_positions, verbose=False):
    """
    Calculate PCA of light distribution (2D: Y-Z plane) - VECTORIZED VERSION.
    Equivalent to GetPCA() in TPCPMTBarycenterMatching_module.cc but ~100x faster.

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

    Returns:
    --------
    pca_vectors : array (n_events, 2)
        Principal component direction [dy, dz] for each event
    pca_elongation : array (n_events,)
        Ratio of largest to smallest eigenvalue (elongation measure)
    """
    n_events = len(f_ophit_PE)

    if verbose:
        print(f">> Calculating Light PCA for {n_events:,} events (vectorized)...")
        import time
        start_time = time.time()

    # Pre-create position lookup arrays for faster access
    # Use max channel ID that could appear in data (not just in pmt_positions)
    max_channel = 312  # SBND has 312 optical channels total
    y_lookup = np.zeros(max_channel, dtype=np.float32)
    z_lookup = np.zeros(max_channel, dtype=np.float32)
    valid_channels = np.zeros(max_channel, dtype=bool)

    for ch, (y, z) in pmt_positions.items():
        if ch < max_channel:
            y_lookup[ch] = y
            z_lookup[ch] = z
            valid_channels[ch] = True

    # Flatten the nested awkward arrays for batch processing
    if verbose:
        print("   * Flattening awkward arrays...")

    flat_pe = ak.flatten(ak.flatten(f_ophit_PE, axis=2), axis=1)
    flat_ch = ak.flatten(ak.flatten(f_ophit_ch, axis=2), axis=1)

    # Convert to numpy for fast indexing
    flat_pe_np = ak.to_numpy(flat_pe)
    flat_ch_np = ak.to_numpy(flat_ch).astype(int)

    # Get event boundaries
    counts = ak.num(ak.flatten(f_ophit_PE, axis=2), axis=1)
    event_ends = np.cumsum(ak.to_numpy(counts))
    event_starts = np.concatenate([[0], event_ends[:-1]])

    if verbose:
        print("   * Computing PCA for all events...")

    # Initialize outputs
    pca_vectors = np.zeros((n_events, 2), dtype=np.float32)
    pca_elongation = np.zeros(n_events, dtype=np.float32)

    # Process in batches for efficiency
    batch_size = 1000
    n_batches = (n_events + batch_size - 1) // batch_size

    for batch_idx in range(n_batches):
        start_ev = batch_idx * batch_size
        end_ev = min((batch_idx + 1) * batch_size, n_events)

        if verbose and batch_idx % 10 == 0:
            elapsed = time.time() - start_time
            progress = start_ev / n_events
            if progress > 0:
                eta = elapsed / progress - elapsed
                print(f"   Batch {batch_idx+1}/{n_batches} | "
                      f"Progress: {100*progress:.1f}% | "
                      f"ETA: {eta:.1f}s")

        # Process each event in batch
        for event_idx in range(start_ev, end_ev):
            start_idx = event_starts[event_idx]
            end_idx = event_ends[event_idx]

            if start_idx >= end_idx:  # No ophits
                pca_vectors[event_idx] = [0.0, 1.0]
                pca_elongation[event_idx] = 1.0
                continue

            # Get channels and PEs for this event
            event_channels = flat_ch_np[start_idx:end_idx]
            event_pes = flat_pe_np[start_idx:end_idx]

            # Filter valid channels using lookup
            mask = (event_channels < max_channel) & valid_channels[event_channels]

            if np.sum(mask) < 2:  # Need at least 2 points
                pca_vectors[event_idx] = [0.0, 1.0]
                pca_elongation[event_idx] = 1.0
                continue

            # Get positions using vectorized lookup
            valid_ch = event_channels[mask]
            valid_pe = event_pes[mask]
            ophit_y = y_lookup[valid_ch]
            ophit_z = z_lookup[valid_ch]

            # Weighted centroid
            weight_sum = np.sum(valid_pe)
            if weight_sum == 0:
                pca_vectors[event_idx] = [0.0, 1.0]
                pca_elongation[event_idx] = 1.0
                continue

            centroid_y = np.sum(valid_pe * ophit_y) / weight_sum
            centroid_z = np.sum(valid_pe * ophit_z) / weight_sum

            # Center points
            centered_y = ophit_y - centroid_y
            centered_z = ophit_z - centroid_z

            # Weighted covariance matrix (vectorized)
            cov_yy = np.sum(valid_pe * centered_y * centered_y) / weight_sum
            cov_yz = np.sum(valid_pe * centered_y * centered_z) / weight_sum
            cov_zz = np.sum(valid_pe * centered_z * centered_z) / weight_sum

            # Eigenvalues and eigenvectors
            cov_matrix = np.array([[cov_yy, cov_yz], [cov_yz, cov_zz]], dtype=np.float32)
            eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)  # eigh is faster for symmetric

            # Get principal component (largest eigenvalue)
            max_idx = 1 if eigenvalues[1] > eigenvalues[0] else 0
            pca_vectors[event_idx] = eigenvectors[:, max_idx]

            # Elongation
            if eigenvalues[0] > 1e-10:
                pca_elongation[event_idx] = eigenvalues[1] / eigenvalues[0]
            else:
                pca_elongation[event_idx] = 100.0

    if verbose:
        elapsed = time.time() - start_time
        print(f">> Light PCA calculation completed in {elapsed:.2f}s")
        print(f"   Average rate: {n_events/elapsed:.0f} events/s")
        print(f"   PCA vector range: Y=[{pca_vectors[:,0].min():.3f}, {pca_vectors[:,0].max():.3f}], "
              f"Z=[{pca_vectors[:,1].min():.3f}, {pca_vectors[:,1].max():.3f}]")
        print(f"   Elongation range: [{pca_elongation.min():.1f}, {pca_elongation.max():.1f}]")

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