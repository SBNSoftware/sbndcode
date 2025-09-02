"""
Configuration settings for PosRecoCNN training pipeline.
Centralized configuration management.
"""

import os
from pathlib import Path

# Get the base directory
BASE_DIR = Path(__file__).parent
PROJECT_ROOT = BASE_DIR.parent

# Data paths configuration
DATA_CONFIG = {
    # Root file paths - update these as needed
    'training_file': '/exp/sbnd/data/users/svidales/AI_nuvT_project_support/mcdata/v10_06_00_02/mc_MCP2025B_02_prodgenie_corsika_proton_rockbox_sbnd_CV_reco2_sbnd_30k_training.root',
    'test_file': '/exp/sbnd/data/users/svidales/AI_nuvT_project_support/mcdata/v10_06_00_02/mc_MCP2025B_02_prodgenie_corsika_proton_rockbox_sbnd_CV_reco2_sbnd_8k_test.root',
    
    # Keys to load from ROOT files
    'keys_to_load': [
        'flash_ophit_pe', 'flash_ophit_ch', 'flash_ophit_time', 
        'nuvT', 'dEpromx', 'dEpromy', 'dEpromz', 'dEtpc', 'nuvZ'
    ],
    
    # PMT map paths (relative to project root)
    'pmt_maps': {
        'coated': PROJECT_ROOT / 'pmt_maps' / 'coated_pmt_map_realistic_flipped.csv',
        'uncoated': PROJECT_ROOT / 'pmt_maps' / 'uncoated_pmt_map_realistic_flipped.csv'
    }
}

# Loading configuration
LOADING_CONFIG = {
    'load_all': True,           # True = all events, False = fraction
    'fraction': 0.5,            # If load_all=False, fraction to load
    'start_event': 0            # Starting event
}

# Event filtering configuration
FILTER_CONFIG = {
    # Energy cuts
    'min_energy_tpc': 50.0,     # MeV
    
    # Position cuts  
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

# Image creation configuration
IMAGE_CONFIG = {
    'max_channels': 312,
    'selection_method': 'max',  # 'max', 'sum', 'nonzero', 'mean_top'
}

# Coordinate scaling configuration
COORD_CONFIG = {
    'ranges': {
        'x': (0, 200),      # x_abs range
        'y': (-200, 200),   # y range  
        'z': (0, 500)       # z range
    }
}

# Training configuration
TRAINING_CONFIG = {
    # Data splits (training uses only train/validation - test is separate)
    'train_fraction': 0.85,  # Increased since no test split needed
    'val_fraction': 0.15,    # Validation for training monitoring
    
    # Model training parameters
    'epochs': 80,
    'batch_size': 32,       # Increased for better gradient estimates
    'patience': 8,          # Increased early stopping patience
    'reduce_lr_patience': 5, # More patience before reducing LR
    'reduce_lr_factor': 0.7, # Less aggressive LR reduction
    'min_lr': 1e-6,
    'dropout_rate': 0.2,    # Reduced dropout to prevent overfitting
    
    # Optimizer
    'optimizer': 'adam',
    'loss': 'mean_squared_error',
    'metrics': ['mse']
}

# Model saving configuration
MODEL_CONFIG = {
    'weights_file': '/tmp/weights_nuvT.hdf5.keras',
    'export_path': '/exp/sbnd/data/users/svidales/AI_nuvT_project_support/cnn_models/v0508_trained_w_30k',
    'model_name_template': 'v{date}_trained_w_{n_events}'  # Will be formatted with date and event count
}

# Visualization configuration
PLOT_CONFIG = {
    'figsize': (20, 10),
    'title_fontsize': 16,
    'colorbar_size': '4%',
    'cmap': 'Blues',
    'use_log_scale': False,
    'show_colorbar': True
}

# Analysis configuration
ANALYSIS_CONFIG = {
    'bias_analysis': {
        'hist_bins': 200,
        'fit_ranges': {
            'X': (65, 130),
            'Y': (65, 130), 
            'Z': (65, 135)
        },
        'hist_ranges': {
            'X': (-250, 250),
            'Y': (-250, 250),
            'Z': (-250, 250)
        }
    },
    
    'reco_truth_comparison': {
        'bins': 80,
        'cmap': 'viridis'
    }
}

# Utility functions for configuration
def get_model_export_path(n_events=None):
    """Generate model export path with current date and event count."""
    from datetime import datetime
    date_str = datetime.now().strftime('%m%d')
    
    if n_events is None:
        model_name = f"v{date_str}_trained"
    else:
        if n_events >= 1000:
            n_str = f"{n_events//1000}k"
        else:
            n_str = str(n_events)
        model_name = MODEL_CONFIG['model_name_template'].format(date=date_str, n_events=n_str)
    
    return os.path.join(MODEL_CONFIG['export_path'], model_name)


def validate_paths():
    """Validate that all required paths exist."""
    missing_paths = []
    
    # Check PMT map files
    for map_name, map_path in DATA_CONFIG['pmt_maps'].items():
        if not map_path.exists():
            missing_paths.append(f"PMT map '{map_name}': {map_path}")
    
    # Check model export directory
    export_dir = Path(MODEL_CONFIG['export_path'])
    if not export_dir.exists():
        print(f">> WARNING: Model export directory doesn't exist: {export_dir}")
        print("   Will be created when saving model.")
    
    if missing_paths:
        print("* ERROR: Missing required files:")
        for path in missing_paths:
            print(f"   {path}")
        return False
    else:
        print(">> SUCCESS: All configuration paths validated successfully")
        return True


def print_config_summary():
    """Print a summary of current configuration."""
    print(">> Configuration Summary")
    print("=" * 50)
    print(f"PMT Maps Directory: {DATA_CONFIG['pmt_maps']['coated'].parent}")
    print(f"Model Export Path: {MODEL_CONFIG['export_path']}")
    print(f"Training Data: {Path(DATA_CONFIG['training_file']).name}")
    print(f"Test Data: {Path(DATA_CONFIG['test_file']).name}")
    print()
    print("Filter Settings:")
    print(f"  Energy cut: >{FILTER_CONFIG['min_energy_tpc']} MeV")
    print(f"  Position ranges: X{FILTER_CONFIG['depromx_range']}, Y{FILTER_CONFIG['depromy_range']}, Z{FILTER_CONFIG['depromz_range']}")
    print()
    print("Training Settings:")
    print(f"  Epochs: {TRAINING_CONFIG['epochs']}")
    print(f"  Batch size: {TRAINING_CONFIG['batch_size']}")
    print(f"  Train/Val split: {TRAINING_CONFIG['train_fraction']:.0%}/{TRAINING_CONFIG['val_fraction']:.0%} (Test data handled separately in inference)")


if __name__ == "__main__":
    print_config_summary()
    validate_paths()