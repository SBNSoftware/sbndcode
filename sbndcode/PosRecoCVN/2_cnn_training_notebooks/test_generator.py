#!/usr/bin/env python3
"""
Test script for ROOTDataGenerator.

Verifies that the generator can load files and produce batches correctly.
Run this before training to ensure everything is set up properly.
"""

import os
import sys
import glob
import numpy as np
import uproot

# Add current directory to path
sys.path.append(os.path.dirname(__file__))

from config import *
from utils import ROOTDataGenerator

def main():
    print("=" * 70)
    print("  ROOTDataGenerator Test Script")
    print("=" * 70)

    # Get file paths
    file_path_config = DATA_CONFIG['training_file']

    if '*' in file_path_config:
        file_paths = sorted(glob.glob(file_path_config))
    else:
        file_paths = [file_path_config]

    print(f"\n>> Found {len(file_paths)} ROOT files:")
    for fp in file_paths:
        print(f"   * {os.path.basename(fp)}")

    if len(file_paths) == 0:
        print("\n!! ERROR: No files found. Check your config.py path.")
        return 1

    # Load channel dictionary
    print("\n>> Loading channel dictionary...")
    with uproot.open(f"{file_paths[0]}:opanatree/PDSMapTree") as pdsmap_tree:
        ids = pdsmap_tree["OpDetID"].array(library="np")
        types = pdsmap_tree["OpDetType"].array(library="np")
        channel_dict = {int(id_val): int(type_val) for id_val, type_val in zip(ids, types)}

    print(f"   * Loaded {len(channel_dict)} channels")

    # Load PMT maps
    print("\n>> Loading PMT maps...")
    uncoated_map = np.loadtxt(DATA_CONFIG['pmt_maps']['uncoated'], delimiter=",", dtype=int)
    coated_map = np.loadtxt(DATA_CONFIG['pmt_maps']['coated'], delimiter=",", dtype=int)
    pmt_maps = (uncoated_map, coated_map)
    print(f"   * Uncoated map shape: {uncoated_map.shape}")
    print(f"   * Coated map shape: {coated_map.shape}")

    # Use only first file for quick test
    test_files = file_paths[:1]
    print(f"\n>> Testing with {len(test_files)} file(s) for quick validation")

    # Create generator
    print("\n>> Creating generator...")
    try:
        generator = ROOTDataGenerator(
            file_paths=test_files,
            batch_size=TRAINING_CONFIG['batch_size'],
            normalization_factor=IMAGE_CONFIG['normalization_factor'],
            channel_dict=channel_dict,
            pmt_maps=pmt_maps,
            filter_config=FILTER_CONFIG,
            coord_config=COORD_CONFIG,
            image_config=IMAGE_CONFIG,
            shuffle=False,  # No shuffle for test
            verbose=True
        )
    except Exception as e:
        print(f"\n!! ERROR creating generator: {e}")
        import traceback
        traceback.print_exc()
        return 1

    print(f"\n>> Generator created successfully!")
    print(f"   * Total events: {generator.total_events:,}")
    print(f"   * Batches: {len(generator)}")
    print(f"   * Batch size: {generator.batch_size}")

    # Test loading first batch
    print("\n>> Testing batch loading...")
    try:
        images, coords = generator[0]
        print(f"   * Batch 0 loaded successfully")
        print(f"     - Images shape: {images.shape}")
        print(f"     - Coords shape: {coords.shape}")
        print(f"     - Images dtype: {images.dtype}")
        print(f"     - Coords dtype: {coords.dtype}")
        print(f"     - Images range: [{np.min(images):.6f}, {np.max(images):.6f}]")
        print(f"     - Coords range: [{np.min(coords):.6f}, {np.max(coords):.6f}]")
    except Exception as e:
        print(f"\n!! ERROR loading batch: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Test loading second batch
    if len(generator) > 1:
        print("\n>> Testing second batch...")
        try:
            images2, coords2 = generator[1]
            print(f"   * Batch 1 loaded successfully")
            print(f"     - Images shape: {images2.shape}")
            print(f"     - Coords shape: {coords2.shape}")
        except Exception as e:
            print(f"\n!! ERROR loading second batch: {e}")
            import traceback
            traceback.print_exc()
            return 1

    # Summary
    print("\n" + "=" * 70)
    print("  ✅ All tests passed!")
    print("=" * 70)
    print("\nThe generator is working correctly. You can now use it in training.")
    print("\nNext steps:")
    print("  1. Update your notebook to use the generator (see example cells)")
    print("  2. Run training with: model.fit(train_generator, validation_data=val_generator)")
    print("  3. Monitor memory usage: watch -n 5 'free -h'")
    print("=" * 70)

    return 0

if __name__ == "__main__":
    sys.exit(main())
