import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import numpy as np

def read_monomer_size(directory_path):
    k_file_path = glob.glob(os.path.join(directory_path, "*.k"))
    if not k_file_path:
        return None
    k_file_path = k_file_path[0]
    with open(k_file_path, 'r') as k_file:
        lines = k_file.readlines()
        monomer_size_line = lines[3]
        monomer_size = float(monomer_size_line.strip().split()[3])
    return int(monomer_size * 1000)

def plot_results(directory_path):
    num_monomers, num_layers = map(int, directory_path.split('_')[-2:])
    monomer_size = read_monomer_size(directory_path)
    if monomer_size is None:
        print(f"Skipping {directory_path}: No .k file found.")
        return

    file_path = os.path.join(directory_path, f'gmm01s_{num_monomers:05d}.out')

    try:
        data = pd.read_csv(file_path, skiprows=18, delim_whitespace=True, header=0)
    except FileNotFoundError:
        print(f"Skipping {directory_path}: {file_path} not found.")
        return

    x = data['s.a.']
    y = data['pol.']
    return {'x': x, 'y': y, 'monomer_size': monomer_size, 'num_layers': num_layers, 'num_monomers': num_monomers}

# Initialize the directory path and lists
agg_models_path = os.path.expanduser('~/runs/Agg_models/NEW_RUNS/')
folders = glob.glob(os.path.join(agg_models_path, "AGG_*_*"))
monomer_sizes_list = [140, 160, 180, 200, 220, 240, 260]
num_layers_list = [2, 3, 4, 5, 6, 7]

# Filter folders
folders = [
    f for f in folders
    if f.split('_')[-2].isdigit() and int(f.split('_')[-2]) in monomer_sizes_list
    and f.split('_')[-1].isdigit() and int(f.split('_')[-1]) in num_layers_list
]

# Initialize a dictionary to accumulate data
accumulated_data = {(size, layers): {'x': [], 'y': []} for size in monomer_sizes_list for layers in num_layers_list}

# Loop through each folder and accumulate the data
for folder in folders:
    plot_data = plot_results(folder)
    if plot_data is None:
        continue
    key = (plot_data['monomer_size'], plot_data['num_layers'])
    if key not in accumulated_data:
        print(f"Unexpected monomer size or number of layers {key}, skipping.")
        continue
    key = (plot_data['monomer_size'], plot_data['num_layers'])
    accumulated_data[key]['x'].append(plot_data['x'])
    accumulated_data[key]['y'].append(plot_data['y'])

# Now process the accumulated data to calculate the mean
average_data = {}
for key, data in accumulated_data.items():
    # Assuming all x values are the same within each group
    if data['x']:
        mean_x = data['x'][0]
        # Stack the y values and calculate the mean across the first axis (vertical stack)
        mean_y = np.vstack(data['y']).mean(axis=0)
        average_data[key] = {'x': mean_x, 'y': mean_y}

# Create subplots
fig, axes = plt.subplots(len(num_layers_list), len(monomer_sizes_list), figsize=(20, 20), sharex='col', sharey='row')

# Plot the average data on the subplots
for (monomer_size, num_layers), plot_data in average_data.items():
    i = num_layers_list.index(num_layers)
    j = monomer_sizes_list.index(monomer_size)
    ax = axes[i, j]

    # Check if there is data to plot
    if 'x' in plot_data and 'y' in plot_data:
        ax.plot(plot_data['x'], plot_data['y'], 'o-')  # Adjust marker style as needed

# Labeling the axes
for j, monomer_size in enumerate(monomer_sizes_list):
    axes[-1, j].set_xlabel(f'Monomer size: {monomer_size}')

for i, num_layers in enumerate(num_layers_list):
    axes[i, 0].set_ylabel(f'Layers: {num_layers}')

# Set title and layout
plt.suptitle("Effect of Number of Layers and size of Monomers on Polarisation", fontsize=24)
plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2, w_pad=2)

# Save and show the figure
plt.savefig('grid_plot_layers_vs_size_mon_comp_final.pdf', format='pdf', bbox_inches='tight')
plt.show()
