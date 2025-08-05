"""
Utility functions for PosRecoCNN training pipeline.
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
        print("🚀 Starting optimized event processing...")
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
    
    # Select TPC values
    tpc_arrays = {
        'dEpromx': arrays['dEpromx'],
        'dEpromy': arrays['dEpromy'],
        'dEpromz': arrays['dEpromz'], 
        'dEtpc': arrays['dEtpc']
    }
    selected_tpc = select_tpc_values(tpc_arrays, decision)
    arrays.update(selected_tpc)
    
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
        print(f"✅ Processing completed in {processing_time:.2f} seconds")
        print()
        tracker.display()
        
        if len(final_arrays['nuvT']) > 0:
            print("\n📊 Final dataset ranges:")
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
    
    print(f"✅ Completed in {processing_time:.2f} seconds")
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
    
    Parameters:
    -----------
    pe_matrix : array (n_events, n_channels)
    *maps : arrays (ch_y, ch_z) - channel mapping arrays  
    method : str - selection method: 'max', 'sum', 'nonzero', 'mean_top'
    normalization_factors : list, optional - pre-computed normalization factors
    
    Returns:
    --------
    images : array (n_events, ch_y//2, ch_z, n_maps) - PE images
    norm_factors : list - normalization factors used (for inference)
    """
    
    ch_y, ch_z = maps[0].shape
    n_events = pe_matrix.shape[0]
    n_maps = len(maps)
    half_y = ch_y // 2
    
    print(f"Creating {n_events:,} images ({half_y}×{ch_z}×{n_maps})")
    
    images = np.zeros((n_events, half_y, ch_z, n_maps), dtype=np.float32)
    norm_factors = []
    
    for map_idx, channel_map in enumerate(maps):
        valid_mask = (channel_map >= 0) & (channel_map < pe_matrix.shape[1])
        pe_spatial = np.zeros((n_events, ch_y, ch_z), dtype=np.float32)
        
        if np.any(valid_mask):
            y_pos, z_pos = np.where(valid_mask)
            channels = channel_map[valid_mask]
            pe_spatial[:, y_pos, z_pos] = pe_matrix[:, channels]
        
        # Normalization
        if normalization_factors is not None:
            # Use provided normalization factor (for inference)
            # First factor is shared between maps 0 and 1
            norm_factor = normalization_factors[0] if map_idx < 2 else normalization_factors[map_idx-1]
            print(f"[Norm] Using provided factor for map {map_idx}: {norm_factor:.3f}")
        else:
            # Calculate normalization factor (for training)
            if map_idx < 2:
                # First 2 maps together
                if map_idx == 0:
                    first_map = pe_spatial.copy()
                    continue
                elif map_idx == 1:
                    norm_factor = max(np.max(first_map), np.max(pe_spatial))
                    print(f"[Norm] Shared max value for maps 0 & 1: {norm_factor:.3f}")
                    norm_factors.append(norm_factor)
                    if norm_factor > 0:
                        first_map /= norm_factor
                        pe_spatial /= norm_factor
                    _select_halves(first_map, images, 0, method, half_y)
            else:
                # Other maps separately
                norm_factor = np.max(pe_spatial)
                print(f"[Norm] Max value for map {map_idx}: {norm_factor:.3f}")
                norm_factors.append(norm_factor)
        
        # Apply normalization and process current map
        if map_idx != 0:
            if norm_factor > 0:
                pe_spatial /= norm_factor
            _select_halves(pe_spatial, images, map_idx, method, half_y)
    
    return images, norm_factors


def scale_coordinates(coordinates, coord_ranges):
    """
    Scale coordinates to normalized ranges.
    
    Parameters:
    -----------
    coordinates : array (n_events, 3) - [x_abs, y, z]
    coord_ranges : dict - coordinate ranges for scaling
    
    Returns:
    --------
    scaled_coords : array (n_events, 3) - scaled coordinates
    """
    x_abs, y, z = coordinates[:, 0], coordinates[:, 1], coordinates[:, 2]
    
    # Scale each coordinate
    x_scaled = (x_abs - coord_ranges['x'][0]) / (coord_ranges['x'][1] - coord_ranges['x'][0])
    y_scaled = (y - coord_ranges['y'][0]) / (coord_ranges['y'][1] - coord_ranges['y'][0]) * 2 - 1
    z_scaled = (z - coord_ranges['z'][0]) / (coord_ranges['z'][1] - coord_ranges['z'][0])
    
    return np.stack([x_scaled, y_scaled, z_scaled], axis=1)


def inverse_scale_coordinates(scaled_coords, coord_ranges):
    """
    Inverse scale coordinates back to original ranges.
    
    Parameters:
    -----------
    scaled_coords : array (n_events, 3) - scaled coordinates
    coord_ranges : dict - coordinate ranges for scaling
    
    Returns:
    --------
    original_coords : array (n_events, 3) - original scale coordinates
    """
    x_inv = scaled_coords[:, 0] * (coord_ranges['x'][1] - coord_ranges['x'][0]) + coord_ranges['x'][0]
    y_inv = ((scaled_coords[:, 1] + 1) / 2) * (coord_ranges['y'][1] - coord_ranges['y'][0]) + coord_ranges['y'][0]
    z_inv = scaled_coords[:, 2] * (coord_ranges['z'][1] - coord_ranges['z'][0]) + coord_ranges['z'][0]
    
    return np.stack([x_inv, y_inv, z_inv], axis=1)


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
    
    print(f"✅ Successfully saved {len(data['nuvT']):,} events to:")
    print(f"   {filename}")
    print(f"   Compressed with snappy for optimal file size")
    if config:
        print(f"   Filter config: {config}")