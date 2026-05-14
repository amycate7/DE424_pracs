# %%
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

GRID_SIZE = 5
T_MAX = 10
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# 2. Go up two levels to get to 'emdw'
# /home/amycate/devel/emdw/src/bin -> src -> emdw
BASE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))

# 3. Define your paths using the calculated BASE_DIR
inferred_path_file = os.path.join(BASE_DIR, 'build', 'wumpus_location1.txt')
actual_path_file   = os.path.join(BASE_DIR, 'src', 'groundtruth', 'wumpus_trajectory.txt')
dataset_dir        = os.path.join(BASE_DIR, 'src', 'dataset1')

# Debug Print: Let's see exactly what Python thinks the path is now
print(f"Calculated BASE_DIR: {BASE_DIR}")
print(f"Targeting build file: {inferred_path_file}")

# 4. Load the data
inferred_path = np.loadtxt(inferred_path_file)
actual_path   = np.loadtxt(actual_path_file)

def load_detections(t):
  # Load the 5x5 binary detection grid from dataset 1 
  filename = os.path.join(dataset_dir, f"data_file00{t}.txt")
  return np.loadtxt(filename).reshape(GRID_SIZE, GRID_SIZE)

fig, ax = plt.subplots(figsize=(7, 7))

def update(t):
  ax.clear()

  # Load detection grid
  det_grid = load_detections(t)

  # Draw detections (blue cells)
  cmap = plt.cm.colors.ListedColormap(['white', 'blue'])
  ax.imshow(det_grid, cmap = cmap, origin = 'upper', extent=[0, 5, 5, 0])

  # Draw grid lines
  ax.set_xticks(np.arange(0, GRID_SIZE + 1, 1))
  ax.set_yticks(np.arange(0, GRID_SIZE + 1, 1))
  ax.grid(which='both', color='gray', linestyle='-', linewidth=1)

  # Plot wumpus actual location (yellow diamond)
  act_x, act_y = actual_path[t]
  ax.scatter(act_x + 0.5, act_y + 0.5, marker='D', color='yellow', s=600, 
               label='Actual', edgecolors='black', zorder=5)

  # Plot inferred wumpus location (red circle)
  inf_x, inf_y = inferred_path[t]
  ax.scatter(inf_x + 0.5, inf_y + 0.5, marker='o', color='red', s=400, 
               label='Inferred', alpha=0.8, edgecolors='darkred', zorder=6)
  # Labels and Formatting
  ax.set_title(f"Wumpus Location Tracking: Time Step {t}", fontsize=14)
  ax.set_xlabel("X coordinate")
  ax.set_ylabel("Y coordinate")
  ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0, fontsize=10, markerscale=0.5, scatterpoints=1, frameon=True, edgecolor='gray')
    
  # Ensure (0,0) is top-left 
  ax.set_xlim(0, 5)
  ax.set_ylim(5, 0)

  fig.tight_layout() # Ensures that the legend isn't cropped out of gif

# Save to a file
wumpus_trajectory = animation.FuncAnimation(fig, update, frames = T_MAX, repeat = True)
wumpus_trajectory.save('wumpus_tracking.gif', writer='pillow', fps=2)
print("Animation saved to wumpus_tracking.gif")