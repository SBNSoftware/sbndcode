"""
Custom loss functions for PosRecoCNN training.
Includes chi2-based loss that combines position and PCA direction.
"""

import tensorflow as tf
import numpy as np


def cosine_similarity_loss(y_true_pca, y_pred_pca):
    """
    Calculate loss based on cosine similarity between true and predicted PCA vectors.

    Parameters:
    -----------
    y_true_pca : tensor, shape (batch, 2)
        True PCA components [pca_y, pca_z]
    y_pred_pca : tensor, shape (batch, 2)
        Predicted PCA components [pca_y, pca_z]

    Returns:
    --------
    loss : tensor, shape (batch,)
        Angular loss in radians
    """
    # Normalize vectors
    y_true_norm = tf.nn.l2_normalize(y_true_pca, axis=-1)
    y_pred_norm = tf.nn.l2_normalize(y_pred_pca, axis=-1)

    # Cosine similarity (dot product of normalized vectors)
    cosine = tf.reduce_sum(y_true_norm * y_pred_norm, axis=-1)

    # Clamp to avoid numerical issues with acos
    cosine = tf.clip_by_value(cosine, -1.0, 1.0)

    # Take absolute value (direction doesn't matter, only alignment)
    cosine = tf.abs(cosine)

    # Convert to angle in radians
    angle_rad = tf.acos(cosine)

    return angle_rad


def chi2_loss_with_pca(x_error=10.0, y_error=10.0, z_error=15.0, angle_error=30.0):
    """
    Create a chi2-based loss function that combines position MSE and PCA angle.

    Based on TPCPMTBarycenterMatching chi2 calculation:
    chi2 = (x-x_true)²/σ_x² + (y-y_true)²/σ_y² + (z-z_true)²/σ_z² + angle²/σ_angle²

    Parameters:
    -----------
    x_error : float
        Expected error in X coordinate (cm)
    y_error : float
        Expected error in Y coordinate (cm)
    z_error : float
        Expected error in Z coordinate (cm)
    angle_error : float
        Expected error in angle (degrees)

    Returns:
    --------
    loss_fn : function
        Loss function compatible with Keras/TensorFlow
    """
    # Convert angle error from degrees to radians
    angle_error_rad = angle_error * (np.pi / 180.0)

    def loss(y_true, y_pred):
        """
        Custom loss function.

        Parameters:
        -----------
        y_true : tensor, shape (batch, 5)
            Ground truth [x, y, z, pca_y, pca_z] (scaled 0-1)
        y_pred : tensor, shape (batch, 5)
            Predictions [x, y, z, pca_y, pca_z] (scaled 0-1)

        Returns:
        --------
        loss : tensor, shape ()
            Mean chi2 loss across batch
        """
        # Split coordinates and PCA
        coords_true = y_true[:, :3]  # [x, y, z]
        coords_pred = y_pred[:, :3]

        pca_true = y_true[:, 3:5]    # [pca_y, pca_z]
        pca_pred = y_pred[:, 3:5]

        # Position chi2 terms (scaled coordinates, so errors are in scaled space)
        # We need to scale errors to match the 0-1 coordinate space
        # For now, use relative weights based on expected errors
        x_weight = 1.0 / (x_error ** 2)
        y_weight = 1.0 / (y_error ** 2)
        z_weight = 1.0 / (z_error ** 2)

        position_loss = (
            x_weight * tf.square(coords_true[:, 0] - coords_pred[:, 0]) +
            y_weight * tf.square(coords_true[:, 1] - coords_pred[:, 1]) +
            z_weight * tf.square(coords_true[:, 2] - coords_pred[:, 2])
        )

        # Angular loss (PCA direction)
        angle_rad = cosine_similarity_loss(pca_true, pca_pred)
        angle_weight = 1.0 / (angle_error_rad ** 2)
        angle_loss = angle_weight * tf.square(angle_rad)

        # Combined chi2
        chi2 = position_loss + angle_loss

        return tf.reduce_mean(chi2)

    return loss


def simple_mse_with_pca():
    """
    Simple MSE over all 5 outputs (position + PCA) without arbitrary weights.

    This treats all outputs equally in scaled coordinate space:
    - Position coordinates: scaled to [0,1] or [-1,1]
    - PCA components: already normalized to [-1,1]

    Returns:
    --------
    loss_fn : function
        Loss function compatible with Keras/TensorFlow
    """
    def loss(y_true, y_pred):
        # Simple MSE over all 5 coordinates
        return tf.reduce_mean(tf.square(y_true - y_pred))

    return loss


def simple_mse_with_pca_and_quality(quality_weight=0.5, sigma_x=5.0, sigma_y=6.5, sigma_z=8.0, sigma_pca=0.1):
    """
    MSE loss with quality score output (6 outputs total).

    Loss = MSE(coords) + quality_weight * MSE(quality_pred - quality_target)

    Quality target is computed on-the-fly as: exp(-error_normalized²/2)
    Error is normalized by expected resolutions (sigma_x, sigma_y, sigma_z, sigma_pca).

    During training:
    - y_true has 5 values: [x, y, z, pca_y, pca_z] (scaled)
    - y_pred has 6 values: [x, y, z, pca_y, pca_z, quality] (scaled + quality ∈ [0,1])

    During inference:
    - Network outputs quality without needing ground truth
    - Quality indicates prediction confidence based on learned patterns

    Parameters:
    -----------
    quality_weight : float
        Weight for quality loss term (default 0.5)
    sigma_x, sigma_y, sigma_z : float
        Expected resolutions in cm for position coordinates
    sigma_pca : float
        Expected resolution for PCA components (in normalized units)

    Returns:
    --------
    loss_fn : function
        Loss function compatible with Keras/TensorFlow
    """
    # Precompute normalization factors
    # We'll work in scaled space, so these are approximate relative weights
    # The actual normalization happens through the exponential

    def loss(y_true, y_pred):
        # y_true: [x, y, z, pca_y, pca_z] (5 values, scaled)
        # y_pred: [x, y, z, pca_y, pca_z, quality] (6 values)

        coords_true = y_true[:, :5]
        coords_pred = y_pred[:, :5]
        quality_pred = y_pred[:, 5]

        # Coordinate MSE (simple MSE over all 5 coords)
        coord_loss = tf.reduce_mean(tf.square(coords_true - coords_pred))

        # Calculate quality target based on prediction error
        # Separate position and PCA for better normalization
        pos_error_sq = tf.square(coords_true[:, :3] - coords_pred[:, :3])  # [x, y, z]
        pca_error_sq = tf.square(coords_true[:, 3:5] - coords_pred[:, 3:5])  # [pca_y, pca_z]

        # Normalized error (approximate, since we're in scaled space)
        # Average over position and PCA separately
        pos_error_normalized = tf.reduce_mean(pos_error_sq, axis=1)
        pca_error_normalized = tf.reduce_mean(pca_error_sq, axis=1)

        # Total normalized error
        total_error = tf.sqrt(pos_error_normalized + pca_error_normalized)

        # Quality target: exp(-error²/2)
        # High quality (→1) when error is small
        # Low quality (→0) when error is large
        quality_target = tf.exp(-total_error / 2.0)

        # Quality loss
        quality_loss = tf.reduce_mean(tf.square(quality_pred - quality_target))

        # Combined loss
        total_loss = coord_loss + quality_weight * quality_loss

        return total_loss

    return loss


def weighted_mse_with_pca(position_weight=1.0, angle_weight=1.0):
    """
    Weighted MSE loss combining position and angle with configurable weights.

    Parameters:
    -----------
    position_weight : float
        Weight for position MSE
    angle_weight : float
        Weight for angle loss

    Returns:
    --------
    loss_fn : function
        Loss function compatible with Keras/TensorFlow
    """
    def loss(y_true, y_pred):
        # Split coordinates and PCA
        coords_true = y_true[:, :3]
        coords_pred = y_pred[:, :3]
        pca_true = y_true[:, 3:5]
        pca_pred = y_pred[:, 3:5]

        # Position MSE
        position_mse = tf.reduce_mean(tf.square(coords_true - coords_pred))

        # Angular loss
        angle_rad = cosine_similarity_loss(pca_true, pca_pred)
        angle_mse = tf.reduce_mean(tf.square(angle_rad))

        # Weighted combination
        total_loss = position_weight * position_mse + angle_weight * angle_mse

        return total_loss

    return loss


# Custom metrics for monitoring
def position_mse(y_true, y_pred):
    """MSE for position only (first 3 coordinates)."""
    coords_true = y_true[:, :3]
    coords_pred = y_pred[:, :3]
    return tf.reduce_mean(tf.square(coords_true - coords_pred))


def angle_mae(y_true, y_pred):
    """Mean absolute angle error in degrees."""
    pca_true = y_true[:, 3:5]
    pca_pred = y_pred[:, 3:5]
    angle_rad = cosine_similarity_loss(pca_true, pca_pred)
    angle_deg = angle_rad * (180.0 / np.pi)
    return tf.reduce_mean(angle_deg)


def angle_mse(y_true, y_pred):
    """MSE for angle in radians."""
    pca_true = y_true[:, 3:5]
    pca_pred = y_pred[:, 3:5]
    angle_rad = cosine_similarity_loss(pca_true, pca_pred)
    return tf.reduce_mean(tf.square(angle_rad))


def quality_mse(y_true, y_pred):
    """
    MSE for quality score (6th output).

    Computes quality target on-the-fly based on prediction error,
    then returns MSE between predicted quality and target quality.
    """
    # Check if quality output exists
    if y_pred.shape[1] < 6:
        return tf.constant(0.0)

    # Extract coordinates and quality prediction
    coords_true = y_true[:, :5]
    coords_pred = y_pred[:, :5]
    quality_pred = y_pred[:, 5]

    # Calculate quality target (same as in loss function)
    pos_error_sq = tf.square(coords_true[:, :3] - coords_pred[:, :3])
    pca_error_sq = tf.square(coords_true[:, 3:5] - coords_pred[:, 3:5])
    pos_error_normalized = tf.reduce_mean(pos_error_sq, axis=1)
    pca_error_normalized = tf.reduce_mean(pca_error_sq, axis=1)
    total_error = tf.sqrt(pos_error_normalized + pca_error_normalized)
    quality_target = tf.exp(-total_error / 2.0)

    # Quality MSE
    return tf.reduce_mean(tf.square(quality_pred - quality_target))


def quality_mean(y_true, y_pred):
    """
    Mean predicted quality (useful for monitoring).
    """
    if y_pred.shape[1] < 6:
        return tf.constant(0.0)

    quality_pred = y_pred[:, 5]
    return tf.reduce_mean(quality_pred)


def huber_loss_with_pca(delta=0.1):
    """
    Huber loss for position + PCA prediction (5 outputs).

    Huber loss is more robust to outliers than MSE:
    - For small errors (|error| <= delta): behaves like MSE (quadratic)
    - For large errors (|error| > delta): behaves like MAE (linear)

    This prevents outliers from dominating the gradient and helps avoid
    the model collapsing to mean predictions for difficult events.

    Parameters:
    -----------
    delta : float
        Threshold for switching from quadratic to linear loss.
        In scaled coordinate space (0-1), typical values:
        - 0.05: Very robust (more like MAE for outliers)
        - 0.1: Balanced (recommended starting point)
        - 0.2: Less robust (more like MSE)

    Returns:
    --------
    loss_fn : function
        Loss function compatible with Keras/TensorFlow
    """

    def loss(y_true, y_pred):
        """
        Huber loss over all 5 coordinates: [x, y, z, pca_y, pca_z]

        Parameters:
        -----------
        y_true : tensor, shape (batch, 5)
            Ground truth [x, y, z, pca_y, pca_z] (scaled 0-1)
        y_pred : tensor, shape (batch, 5)
            Predictions [x, y, z, pca_y, pca_z] (scaled 0-1)

        Returns:
        --------
        loss : tensor, shape ()
            Mean Huber loss across batch and all coordinates
        """
        # Calculate error
        error = y_true - y_pred
        abs_error = tf.abs(error)

        # Huber loss formula:
        # if |error| <= delta: loss = 0.5 * error^2
        # if |error| > delta:  loss = delta * (|error| - 0.5 * delta)
        quadratic = tf.minimum(abs_error, delta)
        linear = abs_error - quadratic
        huber = 0.5 * tf.square(quadratic) + delta * linear

        # Mean over all coordinates and batch
        return tf.reduce_mean(huber)

    return loss


def huber_position_pca_separate(delta_pos=0.1, delta_pca=0.05, weight_pos=1.0, weight_pca=1.0):
    """
    Huber loss with separate deltas and weights for position and PCA.

    This allows different robustness levels for position vs PCA direction,
    which can be useful if they have different error characteristics.

    Parameters:
    -----------
    delta_pos : float
        Huber threshold for position coordinates (x, y, z)
    delta_pca : float
        Huber threshold for PCA components (pca_y, pca_z)
    weight_pos : float
        Weight for position loss
    weight_pca : float
        Weight for PCA loss

    Returns:
    --------
    loss_fn : function
        Loss function compatible with Keras/TensorFlow
    """

    def loss(y_true, y_pred):
        # Split position and PCA
        pos_true = y_true[:, :3]  # [x, y, z]
        pos_pred = y_pred[:, :3]
        pca_true = y_true[:, 3:5]  # [pca_y, pca_z]
        pca_pred = y_pred[:, 3:5]

        # Position Huber loss
        pos_error = tf.abs(pos_true - pos_pred)
        pos_quadratic = tf.minimum(pos_error, delta_pos)
        pos_linear = pos_error - pos_quadratic
        pos_loss = tf.reduce_mean(0.5 * tf.square(pos_quadratic) + delta_pos * pos_linear)

        # PCA Huber loss
        pca_error = tf.abs(pca_true - pca_pred)
        pca_quadratic = tf.minimum(pca_error, delta_pca)
        pca_linear = pca_error - pca_quadratic
        pca_loss = tf.reduce_mean(0.5 * tf.square(pca_quadratic) + delta_pca * pca_linear)

        # Weighted combination
        return weight_pos * pos_loss + weight_pca * pca_loss

    return loss
