#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test file to visualize LinkCSC interpolation with random paths

@author: Claude Sonnet 4.5
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from spatialmath import SE3, SO3
from LinkCSC import LinkCSC

def random_unit_vector():
    """Generate a random unit vector in 3D"""
    vec = np.random.randn(3)
    return vec / np.linalg.norm(vec)

def random_SE3_pose():
    """Generate a random SE3 pose with random position and orientation"""
    # Random position in a reasonable range
    position = np.random.uniform(-5, 5, 3)
    
    # Random orientation using random rotation
    # Generate random rotation by creating a random axis and angle
    axis = random_unit_vector()
    angle = np.random.uniform(0, 2*np.pi)
    
    rotation = SO3.AngleAxis(angle, axis)
    
    return SE3.Rt(rotation, position)

def create_random_link(r=1.0, min_separation=3.0, max_attempts=10):
    """Create a random LinkCSC with reasonable parameters"""
    for attempt in range(max_attempts):
        try:
            start_pose = random_SE3_pose()
            
            # Create end pose with some minimum separation to ensure valid path
            separation = np.random.uniform(min_separation, min_separation + 3)
            end_position = start_pose.t + random_unit_vector() * separation
            end_direction = random_unit_vector()
            
            # Create rotation matrix with end_direction as first column
            # Use Gram-Schmidt to create orthonormal basis
            v1 = end_direction
            v2_temp = random_unit_vector()
            v2 = v2_temp - np.dot(v2_temp, v1) * v1
            v2 = v2 / np.linalg.norm(v2)
            v3 = np.cross(v1, v2)
            
            R_end = np.column_stack([v1, v2, v3])
            end_pose = SE3.Rt(SO3(R_end), end_position)
            
            # Use smaller maxAnglePerElbow to avoid theta >= pi issues
            link = LinkCSC(r, start_pose, end_pose, maxAnglePerElbow=np.pi/3)
            
            # Test that interpolation works
            _ = link.interpolate(count=5)
            
            return link
        except Exception as e:
            if attempt == max_attempts - 1:
                raise RuntimeError(f"Failed to create valid link after {max_attempts} attempts: {e}")
            continue
    
    raise RuntimeError("Failed to create valid link")

def plot_link_interpolation(link, ax, num_points=50, label="", color=None, use_density=False):
    """Plot the interpolated points of a LinkCSC
    
    Parameters:
    -----------
    link : LinkCSC
        The link to interpolate
    ax : matplotlib axis
        The axis to plot on
    num_points : int or float
        If use_density=False: number of points to generate
        If use_density=True: density (points per unit length)
    label : str
        Label for the plot
    color : array-like, optional
        Color for the plot
    use_density : bool
        If True, interpret num_points as density instead of count
    """
    if use_density:
        points = link.interpolate(density=num_points)
    else:
        points = link.interpolate(count=num_points)
    
    if color is None:
        color = np.random.rand(3,)
    
    # Plot the interpolated points
    ax.plot(points[:, 0], points[:, 1], points[:, 2], 
            'o-', label=label, markersize=3, alpha=0.7, color=color)
    
    # Mark start and end points
    ax.plot([link.StartDubinsPose.t[0]], [link.StartDubinsPose.t[1]], [link.StartDubinsPose.t[2]], 
            'go', markersize=10, label=f'{label} Start' if label else 'Start')
    ax.plot([link.EndDubinsPose.t[0]], [link.EndDubinsPose.t[1]], [link.EndDubinsPose.t[2]], 
            'ro', markersize=10, label=f'{label} End' if label else 'End')
    
    return points

def main():
    """Create and visualize several random LinkCSC interpolations"""
    np.random.seed(42)  # For reproducibility
    
    # Create figure with 3D axis
    fig = plt.figure(figsize=(18, 14))
    
    # Create multiple links to compare count vs density approaches
    num_links = 3
    links = []
    
    print("Creating test links...")
    for i in range(num_links):
        try:
            link = create_random_link(r=1.0, min_separation=3.0)
            links.append(link)
            print(f"Successfully created link {i+1}")
        except Exception as e:
            print(f"Error creating link {i+1}: {e}")
    
    # Define varying point counts and densities
    point_counts = [10, 30, 50]  # Low, medium, high point counts
    densities = [5.0, 10.0, 20.0]  # Low, medium, high densities (points per unit)
    
    # Plot each link with varying count and density approaches
    for link_idx, link in enumerate(links):
        total_length = link.lengthC1 + link.lengthS + link.lengthC2
        
        # Row 1: Varying point counts
        for count_idx, count in enumerate(point_counts):
            ax = fig.add_subplot(num_links, 6, link_idx*6 + count_idx + 1, projection='3d')
            try:
                points_count = plot_link_interpolation(link, ax, num_points=count, 
                                                       label=f"Count={count}",
                                                       color=plt.cm.viridis(count_idx/len(point_counts)),
                                                       use_density=False)
                
                info_text = (f"Link {link_idx+1} - Count\n"
                            f"Points requested: {count}\n"
                            f"Points generated: {len(points_count)}\n"
                            f"Length: {total_length:.2f}")
                
                ax.text2D(0.05, 0.95, info_text, transform=ax.transAxes, 
                          fontsize=7, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
                
                ax.set_xlabel('X', fontsize=8)
                ax.set_ylabel('Y', fontsize=8)
                ax.set_zlabel('Z', fontsize=8)
                ax.set_title(f'Link {link_idx+1} - Count={count}', fontsize=9)
                ax.grid(True)
                
                # Set equal aspect ratio
                if len(points_count) > 0:
                    max_range = np.array([points_count[:, 0].max()-points_count[:, 0].min(),
                                         points_count[:, 1].max()-points_count[:, 1].min(),
                                         points_count[:, 2].max()-points_count[:, 2].min()]).max() / 2.0
                    if max_range > 0:
                        mid_x = (points_count[:, 0].max()+points_count[:, 0].min()) * 0.5
                        mid_y = (points_count[:, 1].max()+points_count[:, 1].min()) * 0.5
                        mid_z = (points_count[:, 2].max()+points_count[:, 2].min()) * 0.5
                        ax.set_xlim(mid_x - max_range, mid_x + max_range)
                        ax.set_ylim(mid_y - max_range, mid_y + max_range)
                        ax.set_zlim(mid_z - max_range, mid_z + max_range)
                
            except Exception as e:
                print(f"Error plotting count={count} for link {link_idx+1}: {e}")
                ax.text(0.5, 0.5, 0.5, f"Error: {str(e)[:40]}", 
                       ha='center', va='center', transform=ax.transAxes, fontsize=8)
        
        # Row 2: Varying densities
        for density_idx, density in enumerate(densities):
            ax = fig.add_subplot(num_links, 6, link_idx*6 + 3 + density_idx + 1, projection='3d')
            try:
                points_density = plot_link_interpolation(link, ax, num_points=density, 
                                                         label=f"Density={density}",
                                                         color=plt.cm.plasma(density_idx/len(densities)),
                                                         use_density=True)
                
                expected_points = int(total_length * density)
                info_text = (f"Link {link_idx+1} - Density\n"
                            f"Density: {density} pts/unit\n"
                            f"Points generated: {len(points_density)}\n"
                            f"Expected ≈ {expected_points}")
                
                ax.text2D(0.05, 0.95, info_text, transform=ax.transAxes, 
                          fontsize=7, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
                
                ax.set_xlabel('X', fontsize=8)
                ax.set_ylabel('Y', fontsize=8)
                ax.set_zlabel('Z', fontsize=8)
                ax.set_title(f'Link {link_idx+1} - Density={density}', fontsize=9)
                ax.grid(True)
                
                # Set equal aspect ratio
                if len(points_density) > 0:
                    max_range = np.array([points_density[:, 0].max()-points_density[:, 0].min(),
                                         points_density[:, 1].max()-points_density[:, 1].min(),
                                         points_density[:, 2].max()-points_density[:, 2].min()]).max() / 2.0
                    if max_range > 0:
                        mid_x = (points_density[:, 0].max()+points_density[:, 0].min()) * 0.5
                        mid_y = (points_density[:, 1].max()+points_density[:, 1].min()) * 0.5
                        mid_z = (points_density[:, 2].max()+points_density[:, 2].min()) * 0.5
                        ax.set_xlim(mid_x - max_range, mid_x + max_range)
                        ax.set_ylim(mid_y - max_range, mid_y + max_range)
                        ax.set_zlim(mid_z - max_range, mid_z + max_range)
                
            except Exception as e:
                print(f"Error plotting density={density} for link {link_idx+1}: {e}")
                ax.text(0.5, 0.5, 0.5, f"Error: {str(e)[:40]}", 
                       ha='center', va='center', transform=ax.transAxes, fontsize=8)
    
    plt.tight_layout()
    plt.show()
    
    print("Visualization complete!")

if __name__ == "__main__":
    main()
