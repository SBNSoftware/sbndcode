"""
Configuration settings for Dual TPC CNN Flash Matcher.

Output: 14 values per event
  TPC0: [x0, y0, z0, dx0, dy0, dz0, conf0]
  TPC1: [x1, y1, z1, dx1, dy1, dz1, conf1]
"""

import os
from pathlib import Path

BASE_DIR = Path(__file__).parent
PROJECT_ROOT = BASE_DIR  # Changed: use current directory (copia_de_seguridad) instead of parent

# Data paths
DATA_CONFIG = {
    'training_file': '/exp/sbnd/data/users/svidales/AI_nuvT_project_support/mcdata/v10_06_00_02/mc_MCP2025B_02_prodgenie_corsika_proton_rockbox_sbnd_CV_reco2_sbnd_30k_training.root',
    'test_file': '/exp/sbnd/app/users/svidales/larsoft_v10_06_00_02/barycenter_v0212/opana_tree.root',

    'keys_to_load': [
        # Optical hits for images
        'flash_ophit_pe', 'flash_ophit_ch', 'flash_ophit_time',
        # Filters and selection
        'nuvT',                             # Neutrino time
        'nuvX', 'nuvY', 'nuvZ',             # Neutrino vertex (for active volume filter)
        # Pre-computed from C++ (position, direction, TPC indicator)
        'dEpromx', 'dEpromy', 'dEpromz',    # Centroid (position)
        'dEdirx', 'dEdiry', 'dEdirz',       # PCA direction
        'dEspreadx', 'dEspready', 'dEspreadz',  # PCA spread/dispersion
        'dEtpc',                            # Energy per TPC (>0 = valid)
    ],

    'pmt_maps': {
        'coated': PROJECT_ROOT / 'pmt_maps' / 'coatedPMT_map.csv',
        'uncoated': PROJECT_ROOT / 'pmt_maps' / 'uncoatedPMT_map.csv'
    }
}

# Loading configuration
LOADING_CONFIG = {
    'load_all': True,           # True = all events, False = fraction
    'fraction': 1.0,            # If load_all=False, fraction to load
    'start_event': 0            # Starting event
}

# Event filtering configuration (same as backup PosRecoCVN_Training)
FILTER_CONFIG = {
    # Energy cuts
    'min_energy_tpc': 50.0,     # MeV

    # Active volume cuts (applied BEFORE single neutrino cut)
    'active_volume': {
        'x_range': (-200.0, 200.0),  # cm
        'y_range': (-200.0, 200.0),  # cm
        'z_range': (0.0, 500.0)      # cm
    },

    # Position cuts (for dEprom - applied AFTER energy cuts)
    'depromx_range': (-200.0, 200.0),
    'depromy_range': (-200.0, 200.0),
    'depromz_range': (0.0, 500.0),

    # Flash parameters
    'min_flashes': 1,
    'max_flashes_for_keep_all': 2,
    'invalid_marker': -999.0,

    # Channel types
    'pmt_types': {0, 1},
    'xas_types': {2, 3}
}

# Image: 4 channels
IMAGE_CONFIG = {
    'max_channels': 312,
    'n_channels': 4,
    'channel_names': ['TPC0_uncoated', 'TPC0_coated', 'TPC1_uncoated', 'TPC1_coated'],
    'selection_method': 'max',  # Method for selecting channel values: 'max', 'mean', or 'sum'
}

# Output: 14 values (útil para inferencia)
MODEL_OUTPUT_CONFIG = {
    'n_outputs': 14,
    'tpc0': {'pos': (0, 3), 'dir': (3, 6), 'conf': 6},
    'tpc1': {'pos': (7, 10), 'dir': (10, 13), 'conf': 13}
}

# Coordinate scaling
COORD_CONFIG = {
    'scale_ranges': {
        'x': (0, 200),
        'y': (-200, 200),
        'z': (0, 500),
        'dir': (-1, 1),
    }
}

# Training
TRAINING_CONFIG = {
    'train_fraction': 0.85,
    'val_fraction': 0.15,
    'epochs': 50,
    'batch_size': 32,
    'patience': 5,
    'reduce_lr_patience': 3,
    'reduce_lr_factor': 0.5,
    'min_lr': 1e-6,
    'learning_rate': 1e-3,
    'loss_weights': {'position': 1.0, 'direction': 1.0, 'confidence': 1.0},
    'huber_delta': 0.1,
}

# Model
MODEL_CONFIG = {
    'head_units': 128,
    'dropout_rate': 0.3,
    'weights_file': '/tmp/weights_dual_tpc.keras',
    'export_path': '/exp/sbnd/data/users/svidales',
    'model_name_template': 'v{date}_dual_tpc_{n_events}'
}


def get_model_export_path(n_events=None):
    """Generate model export path with current date."""
    from datetime import datetime
    date_str = datetime.now().strftime('%m%d')
    if n_events is None:
        model_name = f"v{date_str}_dual_tpc"
    else:
        n_str = f"{n_events//1000}k" if n_events >= 1000 else str(n_events)
        model_name = MODEL_CONFIG['model_name_template'].format(date=date_str, n_events=n_str)
    return os.path.join(MODEL_CONFIG['export_path'], model_name)


def validate_paths():
    """Validate required paths exist."""
    for map_name, map_path in DATA_CONFIG['pmt_maps'].items():
        if not map_path.exists():
            print(f"ERROR: Missing {map_name} map: {map_path}")
            return False
    print(">> Paths validated")
    return True


def print_config_summary():
    """Print configuration summary."""
    print("=" * 50)
    print("Dual TPC CNN Flash Matcher")
    print("=" * 50)
    print(f"Input: 4 channels ({', '.join(IMAGE_CONFIG['channel_names'])})")
    print(f"Output: 14 values (7 per TPC: pos + dir + conf)")
    print(f"Training: {TRAINING_CONFIG['epochs']} epochs, batch={TRAINING_CONFIG['batch_size']}")
    print("=" * 50)
