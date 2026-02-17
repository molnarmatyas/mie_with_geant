import numpy as np
import pandas as pd
import cv2
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from mpl_toolkits.mplot3d import Axes3D

def fit_ellipse(data, title="Fitted Ellipse"):
    """Fits an ellipse to a given 2D intensity profile and plots it."""
    points = np.column_stack(np.where(data > data.max() * 0.2))  # Threshold for ellipse fitting
    if len(points) < 5:
        raise ValueError("Not enough points to fit an ellipse")
    ellipse = cv2.fitEllipse(points)
    
    # Plot the intensity data and fitted ellipse
    plt.figure(figsize=(6, 6))
    plt.imshow(data, cmap='gray', origin='lower')
    y0, z0 = ellipse[0]
    a, b = ellipse[1]
    theta = np.deg2rad(ellipse[2])
    
    t = np.linspace(0, 2 * np.pi, 100)
    ellipse_y = y0 + a/2 * np.cos(t) * np.cos(theta) - b/2 * np.sin(t) * np.sin(theta)
    ellipse_z = z0 + a/2 * np.cos(t) * np.sin(theta) + b/2 * np.sin(t) * np.cos(theta)
    
    plt.plot(ellipse_z, ellipse_y, 'r', linewidth=2)
    plt.title(title)
    plt.colorbar()
    plt.show()
    
    # Define zoomed-in region
    zoom_size = 30
    y_min, y_max = int(max(0, y0 - zoom_size)), int(min(data.shape[0], y0 + zoom_size))
    z_min, z_max = int(max(0, z0 - zoom_size)), int(min(data.shape[1], z0 + zoom_size))
    zoomed_data = data[y_min:y_max, z_min:z_max]
    
    # Adjust ellipse points for zoomed-in view
    mask = (ellipse_y >= y_min) & (ellipse_y < y_max) & (ellipse_z >= z_min) & (ellipse_z < z_max)
    zoomed_ellipse_y = ellipse_y[mask]
    zoomed_ellipse_z = ellipse_z[mask]
    
    # Magnified view for ellipse fitting assessment
    plt.figure(figsize=(6, 6))
    plt.imshow(zoomed_data, cmap='gray', origin='lower', extent=[z_min, z_max, y_min, y_max])
    plt.plot(zoomed_ellipse_z, zoomed_ellipse_y, 'r', linewidth=2)
    plt.title(f"Magnified {title}")
    plt.colorbar()
    plt.show()
    
    return ellipse  # (center, (major_axis, minor_axis), angle)

def compute_transformation(ellipse1, ellipse2, dx, pixel_size_mm, align_centers=True, isotropic_scale=True):
    """
    Computes the affine transformation from ellipse1 to ellipse2.
    
    Parameters:
    - ellipse1: First ellipse parameters (center, axes, angle)
    - ellipse2: Second ellipse parameters (center, axes, angle)
    - dx: Distance between measurement planes (in mm)
    - pixel_size_mm: Size of one pixel in mm (e.g., 0.0044 for 4.4 µm pixels)
    - align_centers: If True, aligns ellipse centers before computing transformation.
                     This ensures accurate photon directions even when ellipses are shifted.
    - isotropic_scale: If True, uses the same scale factor for both y and z directions.
                       This is recommended for optical systems where the beam expands uniformly.
                       If False, uses separate scale factors for each axis (may cause distortion).
    
    Returns: (scale_y, scale_z, rotation, translation, dx, center1, center2, align_centers, pixel_size_mm)
    """
    (y1, z1), (a1, b1), theta1 = ellipse1
    (y2, z2), (a2, b2), theta2 = ellipse2
    
    if isotropic_scale:
        # Use geometric mean of the axis ratios for uniform scaling
        scale_avg = np.sqrt((a2 / a1) * (b2 / b1))
        scale_y = scale_avg
        scale_z = scale_avg
        print(f"Using isotropic (uniform) scaling: {scale_avg:.6f}")
    else:
        scale_y, scale_z = a2 / a1, b2 / b1
        print(f"Using anisotropic (non-uniform) scaling")
    
    rotation = R.from_euler('xyz', [0, 0, np.deg2rad(theta2 - theta1)]).as_matrix()
    
    if align_centers:
        # When aligning centers, translation only accounts for the change in ellipse size/rotation
        # The actual center positions are handled separately in the mapping functions
        translation = np.array([0.0, 0.0])
        print("Using center alignment mode:")
        print(f"  Ellipse 1 center: ({y1:.2f}, {z1:.2f}) pixels")
        print(f"  Ellipse 2 center: ({y2:.2f}, {z2:.2f}) pixels")
        print(f"  Center offset: ({y2-y1:.2f}, {z2-z1:.2f}) pixels = ({(y2-y1)*pixel_size_mm:.4f}, {(z2-z1)*pixel_size_mm:.4f}) mm")
    else:
        # Original behavior: translation is the difference in center positions
        translation = np.array([y2 - y1, z2 - z1])
        print("Using original translation mode:")
        print(f"  Translation: {translation} pixels")
    
    print(f"Scale factors: y: {scale_y:.6f}, z: {scale_z:.6f}")
    print(f"Rotation angle (deg): {theta2 - theta1:.6f}")
    print(f"Ellipse 1 axes: a={a1:.2f}, b={b1:.2f} pixels, angle={theta1:.2f}°")
    print(f"Ellipse 2 axes: a={a2:.2f}, b={b2:.2f} pixels, angle={theta2:.2f}°")
    print(f"Pixel size: {pixel_size_mm} mm = {pixel_size_mm*1000} µm")
    print(f"Plane separation: {dx} mm")
    
    return scale_y, scale_z, rotation, translation, dx, (y1, z1), (y2, z2), align_centers, pixel_size_mm

def compute_photon_directions(data1, transform, backwards=False):
    """
    Computes photon directions based on mapping from the first measurement plane to the second.
    Only returns directions for pixels above the intensity threshold.
    
    When align_centers=True, directions are computed relative to the beam center,
    effectively removing artificial shifts from measurement imperfections.
    """
    scale_y, scale_z, rotation, translation, dx, center1, center2, align_centers, pixel_size_mm = transform
    y1_center, z1_center = center1
    y2_center, z2_center = center2
    
    points = np.column_stack(np.where(data1 > data1.max() * 0.2))
    
    directions = []
    mapped_points = []
    for y1, z1 in points:
        if align_centers:
            # Compute position relative to ellipse center (in pixels)
            y_rel = y1 - y1_center
            z_rel = z1 - z1_center
            
            # Position at plane 2 relative to ellipse 2 center (in pixels)
            y_rel2 = scale_y * y_rel
            z_rel2 = scale_z * z_rel
            
            # The displacement in y and z is due to beam expansion only (in pixels)
            dy_pixels = y_rel2 - y_rel  # = (scale_y - 1) * y_rel
            dz_pixels = z_rel2 - z_rel  # = (scale_z - 1) * z_rel
            
            # Convert to physical units (mm)
            dy_mm = dy_pixels * pixel_size_mm
            dz_mm = dz_pixels * pixel_size_mm
            
            # Direction vector: (dx, dy, dz) where dx is in mm
            direction = np.array([dx, dy_mm, dz_mm])
            direction /= np.linalg.norm(direction)
            
            # Store mapped point for visualization (absolute position)
            y2 = y_rel2 + y2_center
            z2 = z_rel2 + z2_center
        else:
            # Original behavior
            y2 = scale_y * (y1 - data1.shape[0] // 2) + translation[0] + data1.shape[0] // 2
            z2 = scale_z * (z1 - data1.shape[1] // 2) + translation[1] + data1.shape[1] // 2
            
            # Convert pixel displacement to mm
            dy_mm = (y2 - y1) * pixel_size_mm
            dz_mm = (z2 - z1) * pixel_size_mm
            
            direction = np.array([dx, dy_mm, dz_mm])
            direction /= np.linalg.norm(direction)
        
        if backwards:
            direction[1] *= -1 # Reverse y component
            direction[2] *= -1 # Reverse z component
        directions.append(direction)
        mapped_points.append((y2, z2))
    
    # Plot the photon mapping
    plt.figure(figsize=(8, 8))
    plt.scatter(points[:, 1], points[:, 0], c='blue', s=1, label='Original Points', alpha=0.5)
    mapped_points = np.array(mapped_points)
    plt.scatter(mapped_points[:, 1], mapped_points[:, 0], c='red', s=1, label='Mapped Points', alpha=0.5)
    if align_centers:
        plt.scatter([z1_center], [y1_center], c='blue', s=100, marker='x', linewidths=3, label='Ellipse 1 Center')
        plt.scatter([z2_center], [y2_center], c='red', s=100, marker='x', linewidths=3, label='Ellipse 2 Center')
    plt.legend()
    plt.title("Photon Mapping from First to Second Plane")
    plt.xlabel("Z coordinate (pixels)")
    plt.ylabel("Y coordinate (pixels)")
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.show()
    
    # Print some statistics about the directions
    directions_arr = np.array(directions)
    print(f"\nDirection statistics:")
    print(f"  Mean direction: ({directions_arr[:, 0].mean():.6f}, {directions_arr[:, 1].mean():.6f}, {directions_arr[:, 2].mean():.6f})")
    print(f"  Max angle from (1,0,0): {np.arccos(directions_arr[:, 0].min()) * 180 / np.pi:.4f}°")
    print(f"  Mean angle from (1,0,0): {np.arccos(directions_arr[:, 0].mean()) * 180 / np.pi:.4f}°")
    
    return directions_arr

def create_full_direction_matrix(data1, transform, backwards=False):
    """
    Creates a full matrix of photon directions for all pixels in data1.
    For pixels above threshold: uses computed photon directions.
    For pixels below threshold: uses straight propagation [1, 0, 0].
    
    Returns: 3D array of shape (height, width, 3) containing direction vectors at each pixel.
    """
    scale_y, scale_z, rotation, translation, dx, center1, center2, align_centers, pixel_size_mm = transform
    y1_center, z1_center = center1
    y2_center, z2_center = center2
    
    # Initialize matrix with straight propagation [1, 0, 0] for all pixels
    direction_matrix = np.zeros((data1.shape[0], data1.shape[1], 3))
    direction_matrix[:, :, 0] = 1.0  # Default: straight propagation in x-direction
    
    # Get points above threshold
    points = np.column_stack(np.where(data1 > data1.max() * 0.2))
    
    # Compute directions for pixels above threshold
    for y1, z1 in points:
        if align_centers:
            # Compute position relative to ellipse center (in pixels)
            y_rel = y1 - y1_center
            z_rel = z1 - z1_center
            
            # Position at plane 2 relative to ellipse 2 center (in pixels)
            y_rel2 = scale_y * y_rel
            z_rel2 = scale_z * z_rel
            
            # The displacement in y and z is due to beam expansion only (in pixels)
            dy_pixels = y_rel2 - y_rel  # = (scale_y - 1) * y_rel
            dz_pixels = z_rel2 - z_rel  # = (scale_z - 1) * z_rel
            
            # Convert to physical units (mm)
            dy_mm = dy_pixels * pixel_size_mm
            dz_mm = dz_pixels * pixel_size_mm
            
            # Direction vector: (dx, dy, dz) where dx is in mm
            direction = np.array([dx, dy_mm, dz_mm])
            direction /= np.linalg.norm(direction)
        else:
            # Original behavior
            y2 = scale_y * (y1 - data1.shape[0] // 2) + translation[0] + data1.shape[0] // 2
            z2 = scale_z * (z1 - data1.shape[1] // 2) + translation[1] + data1.shape[1] // 2
            
            # Convert pixel displacement to mm
            dy_mm = (y2 - y1) * pixel_size_mm
            dz_mm = (z2 - z1) * pixel_size_mm
            
            direction = np.array([dx, dy_mm, dz_mm])
            direction /= np.linalg.norm(direction)
        
        if backwards:
            direction[1] *= -1 # Reverse y component
            direction[2] *= -1 # Reverse z component
        
        direction_matrix[y1, z1] = direction
    
    return direction_matrix

def visualize_photon_directions_3d(data1, ellipse1, directions, transform, sample_rate=1):
    """
    Creates a 3D visualization of photon directions.
    
    Parameters:
    - data1: The initial intensity profile
    - ellipse1: The fitted ellipse parameters
    - directions: The computed photon directions (Nx3 array)
    - transform: The transformation parameters
    - sample_rate: Plot every nth point to avoid overcrowding (default: 1 for all points)
    """
    scale_y, scale_z, rotation, translation, dx, center1, center2, align_centers, pixel_size_mm = transform
    y1_center, z1_center = center1
    y2_center, z2_center = center2
    
    # Extract points where intensity is above threshold
    points = np.column_stack(np.where(data1 > data1.max() * 0.2))
    
    # Sample points if needed
    if sample_rate > 1:
        points = points[::sample_rate]
        directions_sampled = directions[::sample_rate]
    else:
        directions_sampled = directions
    
    # Create 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot original points at x=0 plane (blue) - convert to mm relative to center
    y_rel_mm = (points[:, 0] - y1_center) * pixel_size_mm
    z_rel_mm = (points[:, 1] - z1_center) * pixel_size_mm
    ax.scatter(np.zeros(len(points)), z_rel_mm, y_rel_mm, 
               c='blue', s=10, label='Original Points (x=0)', alpha=0.6)
    
    # Plot mapped points at x=dx plane (red)
    mapped_y_rel_mm = scale_y * y_rel_mm
    mapped_z_rel_mm = scale_z * z_rel_mm
    ax.scatter(np.ones(len(points)) * dx, mapped_z_rel_mm, mapped_y_rel_mm, 
               c='red', s=10, label=f'Mapped Points (x={dx}mm)', alpha=0.6)
    
    # Draw direction arrows (scaled for visibility)
    arrow_scale = dx * 0.5  # Scale relative to plane separation
    for i in range(len(points)):
        direction = directions_sampled[i]
        scaled_dir = direction * arrow_scale
        ax.quiver(0, z_rel_mm[i], y_rel_mm[i], 
                 scaled_dir[0], scaled_dir[2], scaled_dir[1],
                 color='green', arrow_length_ratio=0.1, linewidth=1, alpha=0.7)
    
    # Extract ellipse parameters for visualization
    (y0, z0), (a, b), theta = ellipse1
    
    # Draw ellipse outline at x=0 (in mm relative to center)
    t = np.linspace(0, 2 * np.pi, 100)
    theta_rad = np.deg2rad(theta)
    ellipse_y = (a/2 * np.cos(t) * np.cos(theta_rad) - b/2 * np.sin(t) * np.sin(theta_rad)) * pixel_size_mm
    ellipse_z = (a/2 * np.cos(t) * np.sin(theta_rad) + b/2 * np.sin(t) * np.cos(theta_rad)) * pixel_size_mm
    ax.plot(np.zeros(len(t)), ellipse_z, ellipse_y, 'b-', linewidth=2, label='Ellipse (x=0)')
    
    # Plot beam center
    ax.scatter([0], [0], [0], c='blue', s=100, marker='x', linewidths=3, label='Beam Center (x=0)')
    ax.scatter([dx], [0], [0], c='red', s=100, marker='x', linewidths=3, label=f'Beam Center (x={dx}mm)')
    
    # Set labels and title
    ax.set_xlabel('X (propagation direction) [mm]')
    ax.set_ylabel('Z axis [mm]')
    ax.set_zlabel('Y axis [mm]')
    ax.set_title('3D Photon Directions Visualization (relative to beam center)')
    ax.legend()
    
    # Set aspect ratio and adjust viewing angle
    ax.view_init(elev=20, azim=45)
    
    plt.tight_layout()
    plt.show()


def main(file1, file2, dx, output_file, pixel_size_mm, backwards=False, align_centers=True, isotropic_scale=True):
    """
    Main function to process the data and compute photon directions.
    
    Parameters:
    - file1: Path to first CSV file (e.g., measurement at 10mm)
    - file2: Path to second CSV file (e.g., measurement at 0mm)
    - dx: Distance between measurement planes (in mm)
    - output_file: Path for output CSV file
    - pixel_size_mm: Size of one pixel in mm (e.g., 0.0044 for 4.4 µm pixels)
    - backwards: If True, reverse photon directions
    - align_centers: If True, align ellipse centers before transformation (recommended for shifted ellipses)
    - isotropic_scale: If True, use uniform scaling (same factor for y and z). Recommended for optical systems.
    """
    data1 = pd.read_csv(file1, skiprows=1, header=None).to_numpy()
    data2 = pd.read_csv(file2, skiprows=1, header=None).to_numpy()
    
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.imshow(data1, cmap='gray', origin='lower')
    plt.title("Intensity Profile 1")
    plt.colorbar()
    
    plt.subplot(1, 2, 2)
    plt.imshow(data2, cmap='gray', origin='lower')
    plt.title("Intensity Profile 2")
    plt.colorbar()
    plt.show()
    
    ellipse1 = fit_ellipse(data1, "Fitted Ellipse 1")
    ellipse2 = fit_ellipse(data2, "Fitted Ellipse 2")
    
    transform = compute_transformation(ellipse1, ellipse2, dx, pixel_size_mm,
                                       align_centers=align_centers, 
                                       isotropic_scale=isotropic_scale)
    directions = compute_photon_directions(data1, transform, backwards=backwards)
    
    # Save sparse photon directions (only for pixels above threshold)
    np.savetxt(output_file, directions, fmt='%.6f', delimiter=',')
    print(f"Saved photon directions to {output_file}")
    
    # Create and save full direction matrix for all pixels
    direction_matrix = create_full_direction_matrix(data1, transform, backwards=backwards)
    
    # Save as numpy array (preserves 3D structure, efficient)
    matrix_file_npy = output_file.replace('.csv', '_full_matrix.npy')
    np.save(matrix_file_npy, direction_matrix)
    print(f"Saved full direction matrix to {matrix_file_npy} with shape {direction_matrix.shape}")
    
    # Save as CSV with coordinates (flattened: each row is y, z, dx, dy, dz)
    matrix_file_csv = output_file.replace('.csv', '_full_matrix.csv')
    height, width = data1.shape
    matrix_flat = []
    for y in range(height):
        for z in range(width):
            matrix_flat.append([y, z, direction_matrix[y, z, 0], 
                              direction_matrix[y, z, 1], direction_matrix[y, z, 2]])
    np.savetxt(matrix_file_csv, matrix_flat, fmt='%.6f', delimiter=',', 
               header='y,z,direction_x,direction_y,direction_z', comments='')
    print(f"Saved full direction matrix to {matrix_file_csv} (y, z, dx, dy, dz format)")
    
    # Save as CSV in sequential format (each row corresponds to a row in the original image)
    # Each pixel is represented by 3 consecutive values: direction_x, direction_y, direction_z
    matrix_file_sequential = output_file.replace('.csv', '_full_matrix_sequential.csv')
    sequential_data = []
    for y in range(height):
        row_data = []
        for z in range(width):
            row_data.extend([direction_matrix[y, z, 0], 
                           direction_matrix[y, z, 1], 
                           direction_matrix[y, z, 2]])
        sequential_data.append(row_data)
    np.savetxt(matrix_file_sequential, sequential_data, fmt='%.6f', delimiter=',')
    print(f"Saved sequential direction matrix to {matrix_file_sequential} ({height} rows × {width*3} columns)")

    # Visualize the photon directions in 3D (sample every 10th point to avoid overcrowding)
    visualize_photon_directions_3d(data1, ellipse1, directions, transform, sample_rate=10)
    
    return direction_matrix

# Example usage
# File1 is the first measurement plane, file2 is the second measurement plane
# Maps photons from file1 plane to file2 plane (direction: file1 -> file2)
# 
# IMPORTANT: You need to set pixel_size_mm to match your camera!
# Common values:
#   - 3.45 µm pixels: pixel_size_mm = 0.00345
#   - 4.4 µm pixels:  pixel_size_mm = 0.0044
#   - 5.5 µm pixels:  pixel_size_mm = 0.0055
#
# Set backwards=True to reverse photon directions (from file2 -> file1)
# Set align_centers=True to align ellipse centers before transformation (recommended!)
# Set isotropic_scale=True for uniform scaling (recommended for optical beams)

pixel_size_mm = 0.0045  # <-- CHANGE THIS to match your camera's pixel size in mm

main("N958260004_10.csv", "N958260004_0.csv", dx=10, output_file="photon_directions.csv", 
     pixel_size_mm=pixel_size_mm,
     backwards=False, align_centers=True, isotropic_scale=True) 
# _10 is 10 mm from cell, _0 is at cell
# Here: mapping from 10mm plane to 0mm plane, so photons travel towards the cell
