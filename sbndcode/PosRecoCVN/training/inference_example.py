"""
Example script for running inference with trained PosRecoCNN model.
Demonstrates how to use the training pipeline for new data.
"""

import os
import json
import numpy as np
import awkward as ak
import uproot
import tensorflow as tf
from pathlib import Path

# Import local modules
from config import *
from utils import *

def load_inference_config(model_dir):
    """Load normalization factors and config from trained model directory."""
    norm_file = os.path.join(model_dir, 'normalization_factors.json')
    
    if not os.path.exists(norm_file):
        raise FileNotFoundError(f"Normalization file not found: {norm_file}")
    
    with open(norm_file, 'r') as f:
        config = json.load(f)
    
    return config

def run_inference_pipeline(data_file, model_path, verbose=True):
    """
    Complete inference pipeline for new data.
    
    Parameters:
    -----------
    data_file : str
        Path to ROOT file with new data
    model_path : str  
        Path to exported TensorFlow model directory
    verbose : bool
        Print progress information
        
    Returns:
    --------
    predictions : array
        Predicted coordinates in original scale
    """
    
    if verbose:
        print("🔮 Starting inference pipeline...")
    
    # Load inference configuration
    inference_config = load_inference_config(os.path.dirname(model_path))
    normalization_factors = inference_config['normalization_factors']
    coord_ranges = inference_config['coord_ranges']
    
    if verbose:
        print(f"📋 Loaded normalization factors: {normalization_factors}")
    
    # Load new data
    if verbose:
        print(f"📂 Loading data from: {Path(data_file).name}")
    
    file = uproot.open(data_file)
    optree = file['opanatree']['OpAnaTree']
    
    # Load required arrays
    arrays = [optree[key].array() for key in DATA_CONFIG['keys_to_load']]
    f_ophit_PE, f_ophit_ch, f_ophit_t, nuvT, dEpromx, dEpromy, dEpromz, dEtpc, nuvZ = arrays
    
    # Create channel dictionary
    PDSMap = file['opanatree']['PDSMapTree']
    ID = PDSMap['OpDetID'].array()
    Type = PDSMap['OpDetType'].array()
    channel_dict = {id_val: int(type_val) for id_val, type_val in zip(ID[0], Type[0])}
    
    # Process events (apply same filters as training)
    if verbose:
        print("🔄 Processing events...")
    
    results = process_events(
        nuvT, f_ophit_PE, f_ophit_ch, f_ophit_t,
        dEpromx, dEpromy, dEpromz, dEtpc, nuvZ,
        channel_dict, FILTER_CONFIG, verbose=verbose
    )
    
    # Unpack results
    (nuvT_final, f_ophit_PE_final, f_ophit_ch_final, f_ophit_t_final,
     dEpromx_final, dEpromy_final, dEpromz_final, dEtpc_final, nuvZ_final, stats) = results
    
    if len(nuvT_final) == 0:
        print("⚠️  No events passed filters!")
        return np.array([])
    
    # Create PE matrix
    if verbose:
        print("🔢 Creating PE matrix...")
    pe_matrix = create_pe_matrix(f_ophit_PE_final, f_ophit_ch_final, IMAGE_CONFIG['max_channels'])
    
    # Load PMT maps
    uncoated_map = np.loadtxt(DATA_CONFIG['pmt_maps']['uncoated'], delimiter=",", dtype=int)
    coated_map = np.loadtxt(DATA_CONFIG['pmt_maps']['coated'], delimiter=",", dtype=int)
    
    # Create images using SAVED normalization factors (critical for inference!)
    if verbose:
        print("🖼️ Creating PE images with saved normalization...")
    
    images, _ = create_pe_images(
        pe_matrix, uncoated_map, coated_map,
        method=IMAGE_CONFIG['selection_method'],
        normalization_factors=normalization_factors  # Use training normalization!
    )
    
    # Load model and predict
    if verbose:
        print(f"🤖 Loading model from: {Path(model_path).name}")
    
    model = tf.saved_model.load(model_path)
    
    if verbose:
        print(f"🔮 Making predictions for {len(images)} events...")
    
    # Make predictions (model expects float32)
    predictions_scaled = model(images.astype(np.float32)).numpy()
    
    # Convert back to original coordinates
    predictions_original = inverse_scale_coordinates(predictions_scaled, coord_ranges)
    
    if verbose:
        print("✅ Inference completed successfully")
        print(f"📊 Predicted coordinate ranges:")
        print(f"  X: [{np.min(predictions_original[:, 0]):.1f}, {np.max(predictions_original[:, 0]):.1f}] cm")
        print(f"  Y: [{np.min(predictions_original[:, 1]):.1f}, {np.max(predictions_original[:, 1]):.1f}] cm") 
        print(f"  Z: [{np.min(predictions_original[:, 2]):.1f}, {np.max(predictions_original[:, 2]):.1f}] cm")
    
    return predictions_original, (dEpromx_final, dEpromy_final, dEpromz_final)

def compare_predictions_with_truth(predictions, truth_coords, coord_names=['X', 'Y', 'Z']):
    """Compare predictions with truth values and calculate metrics."""
    
    # Convert truth to numpy arrays
    truth_array = np.column_stack([
        np.abs(np.array(truth_coords[0]).flatten()),  # X absolute
        np.array(truth_coords[1]).flatten(),          # Y
        np.array(truth_coords[2]).flatten()           # Z
    ])
    
    # Calculate metrics
    diff = predictions - truth_array
    mse = np.mean(diff**2, axis=0)
    rmse = np.sqrt(mse)
    mae = np.mean(np.abs(diff), axis=0)
    
    print("\n📊 Inference Performance Metrics:")
    print("=" * 40)
    for i, coord in enumerate(coord_names):
        print(f"{coord} coordinate:")
        print(f"  RMSE: {rmse[i]:.2f} cm")
        print(f"  MAE:  {mae[i]:.2f} cm")
        print(f"  Bias: {np.mean(diff[:, i]):.2f} cm")
        print()
    
    overall_rmse = np.sqrt(np.mean(mse))
    print(f"Overall RMSE: {overall_rmse:.2f} cm")
    
    return {
        'rmse': rmse,
        'mae': mae,
        'bias': np.mean(diff, axis=0),
        'overall_rmse': overall_rmse
    }

if __name__ == "__main__":
    # Example usage
    print("🎯 PosRecoCNN Inference Example")
    print("=" * 40)
    
    # Update these paths for your setup
    TEST_FILE = DATA_CONFIG['test_file']
    MODEL_PATH = get_model_export_path()  # This will generate the expected path
    
    print(f"Data file: {Path(TEST_FILE).name}")
    print(f"Model path: {Path(MODEL_PATH).name}")
    
    # Check if model exists
    if not os.path.exists(MODEL_PATH):
        print(f"❌ Model not found at: {MODEL_PATH}")
        print("💡 Train a model first using the training notebook")
        exit(1)
    
    try:
        # Run inference
        predictions, truth_coords = run_inference_pipeline(TEST_FILE, MODEL_PATH)
        
        if len(predictions) > 0:
            # Compare with truth
            metrics = compare_predictions_with_truth(predictions, truth_coords)
            print("\n✅ Inference example completed successfully!")
        else:
            print("⚠️  No valid predictions to evaluate")
            
    except Exception as e:
        print(f"❌ Error during inference: {e}")
        raise