"""
Chang's method for SNR estimation in MRI images.

Reference:
Chang, L.-C., Rohde, G. K., & Pierpaoli, C. (2005).
An automatic method for estimating noise-induced signal variance in 
magnitude-reconstructed magnetic resonance images.
Medical Imaging, 1136-1142.
"""
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Optional
from scipy import stats


np.seterr(divide='ignore', invalid='ignore')


def gaussian_kernel(x: np.ndarray) -> np.ndarray:
    """
    Gaussian kernel function for kernel density estimation.
    
    Args:
        x: Input array
        
    Returns:
        Gaussian kernel values
    """
    return np.exp(-0.5 * x**2) / np.sqrt(2 * np.pi)


def estimate_noise_std_histogram(img_normalized: np.ndarray, 
                                  n_bins: int) -> float:
    """
    Estimate noise standard deviation from histogram peak.
    
    Args:
        img_normalized: Normalized image data (0-1 range)
        n_bins: Number of histogram bins
        
    Returns:
        Estimated standard deviation (normalized)
    """
    bin_counts, bin_edges = np.histogram(img_normalized, bins=n_bins)
    
    if bin_counts.max() == 0:
        return 0.0
    
    # Find the mode (most frequent value)
    max_bin_idx = np.argmax(bin_counts)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    return bin_centers[max_bin_idx]


def kernel_density_estimation(data: np.ndarray, 
                              eval_points: np.ndarray,
                              bandwidth: float,
                              sample_size: Optional[int] = 10000) -> np.ndarray:
    """
    Perform kernel density estimation using Gaussian kernel.
    
    Optimized version that samples data if too large for performance.
    
    Args:
        data: Input data points
        eval_points: Points at which to evaluate the density
        bandwidth: Bandwidth parameter for KDE
        sample_size: Maximum number of samples to use (None for all data)
        
    Returns:
        Density estimates at eval_points
    """
    n = len(data)
    
    # Sample data if too large for performance
    if sample_size and n > sample_size:
        indices = np.random.choice(n, sample_size, replace=False)
        data_sampled = data[indices]
        n_effective = sample_size
    else:
        data_sampled = data
        n_effective = n
    
    # Vectorized KDE computation
    # Shape: (n_eval_points, n_data_points)
    diff = (eval_points[:, np.newaxis] - data_sampled[np.newaxis, :]) / bandwidth
    
    # Apply Gaussian kernel and sum
    kernel_values = gaussian_kernel(diff)
    density = kernel_values.sum(axis=1) / (n_effective * bandwidth)
    
    return density


def calculate_snr_chang(img: np.ndarray, 
                       histogram_factor: float = 1.0,
                       show_plot: bool = False,
                       sample_size: Optional[int] = 10000) -> Tuple[np.ndarray, float, float]:
    """
    Calculate SNR map using Chang's method.
    
    Args:
        img: Input MRI image (2D or 3D array)
        histogram_factor: Factor to multiply number of histogram bins
        show_plot: Whether to display the SNR map
        sample_size: Number of samples for KDE (None to use all data)
        
    Returns:
        Tuple of (snr_map, estimated_std, estimated_std_normalized)
    """
    # Normalize image to [0, 1]
    img_float = img.astype(np.float64)
    img_max = img_float.max()
    
    if img_max == 0:
        return np.zeros_like(img), 0.0, 0.0
    
    img_flat = img_float.ravel()
    img_normalized = img_flat / img_max
    
    # Calculate number of bins using Sturges' rule (improved)
    n_pixels = img_normalized.size
    n_bins = int(np.ceil(np.sqrt(n_pixels)) * histogram_factor)
    
    # Initial noise estimate from histogram
    initial_std = estimate_noise_std_histogram(img_normalized, n_bins)
    
    # Silverman's rule of thumb for bandwidth
    bandwidth = 1.06 * n_pixels**(-0.2) * initial_std
    
    # Avoid zero bandwidth
    if bandwidth < 1e-6:
        bandwidth = 0.01
    
    # Kernel density estimation
    eval_points = np.linspace(0, 1, n_bins)
    density = kernel_density_estimation(
        img_normalized, 
        eval_points, 
        bandwidth,
        sample_size
    )
    
    # Find mode of the density estimate
    max_density_idx = np.argmax(density)
    estimated_std_normalized = eval_points[max_density_idx]
    
    # Scale back to original image units
    # Note: Division by 10 is from original implementation
    estimated_std = estimated_std_normalized * img_max / 10
    
    # Calculate SNR map
    # SNR = sqrt(|S^2 - σ^2|) / σ
    signal_squared = img_float**2
    noise_squared = estimated_std**2
    
    snr_map = np.sqrt(np.abs(signal_squared - noise_squared)) / estimated_std
    
    # Handle potential infinities
    snr_map = np.nan_to_num(snr_map, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Optional plotting
    if show_plot:
        _plot_snr_map(snr_map, img.shape)
    
    return snr_map, estimated_std, estimated_std_normalized


def _plot_snr_map(snr_map: np.ndarray, original_shape: tuple) -> None:
    """
    Plot SNR map for 2D or 3D images.
    
    Args:
        snr_map: SNR map array
        original_shape: Original image shape
    """
    plt.figure(figsize=(10, 8))
    
    if len(original_shape) == 2:
        plt.imshow(snr_map, cmap='viridis')
        plt.colorbar(label='SNR')
        plt.title("SNR Map (Chang's Method)")
    elif len(original_shape) == 3:
        # Show middle slice
        snr_reshaped = snr_map.reshape(original_shape)
        middle_slice = original_shape[2] // 2
        plt.imshow(snr_reshaped[:, :, middle_slice], cmap='viridis')
        plt.colorbar(label='SNR')
        plt.title(f"SNR Map - Slice {middle_slice} (Chang's Method)")
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.tight_layout()
    plt.show()


def calculate_snr_scipy(img: np.ndarray, 
                       show_plot: bool = False) -> Tuple[np.ndarray, float, float]:
    """
    Alternative implementation using scipy's gaussian_kde for comparison.
    This is often faster for large datasets.
    
    Args:
        img: Input MRI image
        show_plot: Whether to display the SNR map
        
    Returns:
        Tuple of (snr_map, estimated_std, estimated_std_normalized)
    """
    img_float = img.astype(np.float64)
    img_max = img_float.max()
    
    if img_max == 0:
        return np.zeros_like(img), 0.0, 0.0
    
    img_normalized = (img_float / img_max).ravel()
    
    # Use scipy's gaussian_kde (Scott's rule for bandwidth)
    kde = stats.gaussian_kde(img_normalized, bw_method='scott')
    
    # Evaluate density
    eval_points = np.linspace(0, 1, 1000)
    density = kde(eval_points)
    
    # Find mode
    max_density_idx = np.argmax(density)
    estimated_std_normalized = eval_points[max_density_idx]
    estimated_std = estimated_std_normalized * img_max / 10
    
    # Calculate SNR map
    snr_map = np.sqrt(np.abs(img_float**2 - estimated_std**2)) / estimated_std
    snr_map = np.nan_to_num(snr_map, nan=0.0, posinf=0.0, neginf=0.0)
    
    if show_plot:
        _plot_snr_map(snr_map, img.shape)
    
    return snr_map, estimated_std, estimated_std_normalized


# Legacy function name for backward compatibility
def calcSNR(img: np.ndarray, show: int = 0, fac: float = 1.0) -> Tuple[np.ndarray, float, float]:
    """
    Legacy function signature for backward compatibility.
    
    Args:
        img: Input image
        show: Whether to show plot (0 or 1)
        fac: Histogram factor
        
    Returns:
        Tuple of (snr_map, estimated_std, estimated_std_normalized)
    """
    return calculate_snr_chang(img, histogram_factor=fac, show_plot=bool(show))


if __name__ == "__main__":
    # Example usage and performance comparison
    print("Testing SNR calculation methods...")
    
    # Create synthetic test image
    np.random.seed(42)
    test_img = np.random.randn(128, 128) * 10 + 100
    test_img = np.abs(test_img)  # Magnitude image
    
    # Test custom implementation
    import time
    start = time.time()
    snr_map1, std1, std_norm1 = calculate_snr_chang(test_img, show_plot=False)
    time1 = time.time() - start
    
    # Test scipy implementation
    start = time.time()
    snr_map2, std2, std_norm2 = calculate_snr_scipy(test_img, show_plot=False)
    time2 = time.time() - start
    
    print(f"\nCustom implementation:")
    print(f"  Time: {time1:.4f}s")
    print(f"  Estimated STD: {std1:.4f}")
    print(f"  Mean SNR: {snr_map1.mean():.4f}")
    
    print(f"\nScipy implementation:")
    print(f"  Time: {time2:.4f}s")
    print(f"  Estimated STD: {std2:.4f}")
    print(f"  Mean SNR: {snr_map2.mean():.4f}")
