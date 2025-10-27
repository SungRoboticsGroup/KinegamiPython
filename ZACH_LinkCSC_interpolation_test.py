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

def plot_link_interpolation(link, ax, num_points=50, label="", color=None):
    """Plot the interpolated points of a LinkCSC"""
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
    fig = plt.figure(figsize=(14, 10))
    
    # Create multiple subplots
    num_links = 4
    
    for i in range(num_links):
        ax = fig.add_subplot(2, 2, i+1, projection='3d')
        
        try:
            # Create a random link
            print(f"Creating link {i+1}...")
            link = create_random_link(r=1.0, min_separation=3.0)
            
            # Plot interpolation
            points = plot_link_interpolation(link, ax, num_points=50, 
                                            label=f"Link {i+1}", 
                                            color=plt.cm.viridis(i/num_links))
            
            # Add some info text
            info_text = (f"Link {i+1}\n"
                        f"r = {link.r:.2f}\n"
                        f"Length: {link.lengthC1 + link.lengthS + link.lengthC2:.2f}\n"
                        f"θ1 = {np.rad2deg(link.path.theta1):.1f}°\n"
                        f"θ2 = {np.rad2deg(link.path.theta2):.1f}°")
            
            ax.text2D(0.05, 0.95, info_text, transform=ax.transAxes, 
                     fontsize=8, verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_title(f'LinkCSC Interpolation {i+1}')
            ax.legend(fontsize=8)
            ax.grid(True)
            
            # Set equal aspect ratio
            max_range = np.array([points[:, 0].max()-points[:, 0].min(),
                                 points[:, 1].max()-points[:, 1].min(),
                                 points[:, 2].max()-points[:, 2].min()]).max() / 2.0
            mid_x = (points[:, 0].max()+points[:, 0].min()) * 0.5
            mid_y = (points[:, 1].max()+points[:, 1].min()) * 0.5
            mid_z = (points[:, 2].max()+points[:, 2].min()) * 0.5
            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)
            
            print(f"Successfully created link {i+1}")
            
        except Exception as e:
            print(f"Error creating link {i+1}: {e}")
            ax.text(0.5, 0.5, 0.5, f"Error: {str(e)[:50]}", 
                   ha='center', va='center', transform=ax.transAxes)
    
    plt.tight_layout()
    plt.show()
    
    print("Visualization complete!")

if __name__ == "__main__":
    main()
