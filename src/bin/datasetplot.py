# %%
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
import sys 

# python3 datasetplot.py <dataset_num>
if len(sys.argv) < 2:
    print("Usage: python datasetplot.py <dataset_num>")
    sys.exit(1)

dataset_num = int(sys.argv[1])

if dataset_num == 1:
  GRID_R = 5
  GRID_C = 5
  T_MAX = 10
  inf_file = 'wumpus_location1.txt'
  traj_file = 'wumpus_trajectory1.txt' # Ground truth for D1
  d_dir = 'dataset1'
elif dataset_num == 2:
  GRID_R = 20
  GRID_C = 20
  T_MAX = 20
  inf_file = 'wumpus_location2.txt'
  traj_file = 'wumpus_trajectory2.txt' # Ground truth for D2
  d_dir = 'dataset2'
elif dataset_num == 3:
  GRID_R = 10
  GRID_C = 20
  T_MAX = 20
  inf_file = 'wumpus_location3.txt'
  traj_file = 'wumpus_trajectory3.txt'
  d_dir = 'dataset3'


print(f"Visualizing Dataset {dataset_num}")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up two levels to get to 'emdw'
BASE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))

# Define your paths using the calculated BASE_DIR
inferred_path_file = os.path.join(BASE_DIR, 'build', inf_file)
actual_path_file   = os.path.join(BASE_DIR, 'src', 'groundtruth', traj_file)
dataset_dir        = os.path.join(BASE_DIR, 'src', d_dir)

# tO see exactly what Python thinks the path is now
print(f"Calculated BASE_DIR: {BASE_DIR}")
print(f"Targeting build file: {inferred_path_file}")

# 4. Load the data
inferred_path = np.loadtxt(inferred_path_file)
actual_path   = np.loadtxt(actual_path_file)

def load_detections(t):
    # zfill(3) ensures t=1 becomes '001' and t=10 becomes '010'
    filename = f"data_file{str(t).zfill(3)}.txt"
    path = os.path.join(dataset_dir, filename)
    
    if os.path.exists(path):
        # Load and reshape based on the current dataset's shape
        return np.loadtxt(path).reshape(GRID_R, GRID_C)
    else:
        raise FileNotFoundError(f"Missing detection file: {path}")

def calculate_metrics(inferred, actual):
    # Accuracy: % of time steps with an exact (x, y) match
    matches = np.all(inferred == actual, axis=1) # np.all(..., axis=1) checks if both x and y match for each row
    accuracy = np.mean(matches) * 100
    
    # Euclidean Distance Error for each time step
    distances = np.linalg.norm(inferred - actual, axis=1)
    rmse = np.sqrt(np.mean(distances**2))
    
    return accuracy, rmse, matches

accuracy, rmse, matches = calculate_metrics(inferred_path, actual_path)

print(f"\n--- Quantitative Analysis: Dataset {dataset_num} ---")
print(f"Overall Accuracy: {accuracy:.2f}%")
print(f"Root Mean Square Error (RMSE): {rmse:.4f} units")

def calculate_spatial_metrics(inferred, actual):
    # Calculate Manhattan distance for every time step
    # L1 norm: |x1-x2| + |y1-y2|
    manhattan_distances = np.sum(np.abs(inferred - actual), axis=1)
    
    # Neighborhood Accuracy: % of time within 1 cell of the truth
    neighborhood_acc = np.mean(manhattan_distances <= 1) * 100
    
    # Mean Manhattan Distance (MMD)
    mmd = np.mean(manhattan_distances)
    
    return neighborhood_acc, mmd

# --- Integration ---
n_acc, mmd = calculate_spatial_metrics(inferred_path, actual_path)
print(f"Neighborhood Accuracy (within 1 cell): {n_acc:.2f}%")
print(f"Mean Manhattan Distance (MMD): {mmd:.2f} steps")

# Print which specific time steps failed
failed_steps = np.where(~matches)[0]
if len(failed_steps) > 0:
    print(f"Deviations detected at time steps: {failed_steps}")
else:
    print("Perfect tracking achieved!")

fig, ax = plt.subplots(figsize=(7, 7))

def update(t):
    ax.clear()

    # Load detection grid (ensure this function returns a shape of (GRID_R, GRID_C))
    det_grid = load_detections(t)

    # Draw detections (blue cells)
    cmap = plt.cm.colors.ListedColormap(['white', 'blue'])
    # extent=[left, right, bottom, top]
    ax.imshow(det_grid, cmap=cmap, origin='upper', extent=[0, GRID_C, GRID_R, 0])

    # Draw grid lines based on columns and rows
    ax.set_xticks(np.arange(0, GRID_C + 1, 1))
    ax.set_yticks(np.arange(0, GRID_R + 1, 1))
    ax.grid(which='both', color='gray', linestyle='-', linewidth=0.5)

    # Plot wumpus actual location (yellow diamond)
    # Scaled size using GRID_C as the reference for marker consistency
    act_x, act_y = actual_path[t]
    ax.scatter(act_x + 0.5, act_y + 0.5, marker='D', color='yellow', s=(3000/GRID_C), 
               label='Actual', edgecolors='black', zorder=5)

    # Plot inferred wumpus location (red circle)
    inf_x, inf_y = inferred_path[t]
    ax.scatter(inf_x + 0.5, inf_y + 0.5, marker='o', color='red', s=(2000/GRID_C), 
               label='Inferred', alpha=0.8, edgecolors='darkred', zorder=6)

    # Labels and Formatting
    ax.set_title(f"Wumpus Location Tracking: Time Step {t}", fontsize=14)
    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Y coordinate")
    
    # Place legend outside the plot
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0, 
              fontsize=10, markerscale=0.5, scatterpoints=1, frameon=True, edgecolor='gray')
    
    # Ensure limits match the grid dimensions
    ax.set_xlim(0, GRID_C)
    ax.set_ylim(GRID_R, 0) # Flipped for top-left origin

    fig.tight_layout()


# Save to a file
wumpus_trajectory = animation.FuncAnimation(fig, update, frames = T_MAX, repeat = True)
wumpus_trajectory.save('wumpus_tracking.gif', writer='pillow', fps=1)
print("Animation saved to wumpus_tracking.gif")