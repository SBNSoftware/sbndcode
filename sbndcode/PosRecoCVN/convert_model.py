#!/usr/bin/env python3
"""
Script to convert Keras model to TensorFlow SavedModel format
Run this in your SBND environment where TensorFlow is available
"""

import tensorflow as tf
import os
import sys

def convert_keras_to_savedmodel():
    keras_path = "/exp/sbnd/data/users/svidales/AI_nuvT_project_support/cnn_models/model_resnet18_50energy_pmts_v2406_part1_reentrenarconpart2.keras"
    savedmodel_path = "/exp/sbnd/data/users/svidales/AI_nuvT_project_support/cnn_models/resnet18_modeltest_v2406"
    
    print("🔄 Converting Keras model to SavedModel format")
    print(f"Input:  {keras_path}")
    print(f"Output: {savedmodel_path}")
    print("=" * 70)
    
    try:
        # Load Keras model
        print("Loading Keras model...")
        model = tf.keras.models.load_model(keras_path)
        print("✅ Model loaded successfully")
        
        # Print model info
        print("\n📋 Model Summary:")
        model.summary()
        
        print(f"\n🔍 Model Details:")
        print(f"Input shape: {model.input_shape}")
        print(f"Output shape: {model.output_shape}")
        print(f"Model type: {type(model)}")
        
        # Create output directory
        os.makedirs(os.path.dirname(savedmodel_path), exist_ok=True)
        
        # Save as SavedModel
        print(f"\n💾 Saving as SavedModel...")
        model.save(savedmodel_path, save_format='tf')
        print("✅ SavedModel created successfully!")
        
        # Verify saved model
        print(f"\n🔍 Verifying saved model...")
        loaded_model = tf.saved_model.load(savedmodel_path)
        
        # Check signatures
        signatures = list(loaded_model.signatures.keys())
        print(f"Available signatures: {signatures}")
        
        if 'serving_default' in signatures:
            serving_default = loaded_model.signatures['serving_default']
            input_names = list(serving_default.structured_input_signature[1].keys())
            output_names = list(serving_default.structured_outputs.keys())
            print(f"Input tensor names: {input_names}")
            print(f"Output tensor names: {output_names}")
        
        print("\n" + "=" * 70)
        print("🎉 Conversion completed successfully!")
        print(f"\n📁 SavedModel created at: {savedmodel_path}")
        print(f"\n⚙️  Update your FCL with:")
        print(f'   ModelPath: "{savedmodel_path}"')
        print(f'   RunInference: true')
        
        return True
        
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        return False

if __name__ == "__main__":
    convert_keras_to_savedmodel()