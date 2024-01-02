import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def read_monomer_size(directory_path):
    k_file_path = glob.glob(os.path.join(directory_path, "*.k"))
    if not k_file_path:
        print(f"Skipping {directory_path}: No .k file found.")
        return None
    k_file_path = k_file_path[0]
    with open(k_file_path, 'r') as k_file:
        first_line = k_file.readline().strip()
        print(f"Reading monomer size from {k_file_path}: {first_line}")  # Debugging line
        monomer_size = float(k_file.readline().strip())
    return monomer_size

def plot_results(directory_path):
    num_monomers, num_layers = map(int, directory_path.split('_')[-2:])
    monomer_size = read_monomer_size(directory_path)
    if monomer_size is None:
        return

    file_path = glob.glob(os.path.join(directory_path, "gmm01s_00000.out"))
    if not file_path:
        print(f"Skipping {directory_path}: No .out file found.")
        return
    file_path = file_path[0]

    try:
        data = pd.read_csv(file_path, skiprows=18, delim_whitespace=True, header=0)
    except pd.errors.ParserError:
        print(f"Error parsing file: {file_path}")
        return
    except FileNotFoundError:
        print(f"Skipping {directory_path}: {file_path} not found.")
        return

    return {'x': data['s.a.'], 'y': data['pol.'], 'monomer_size': monomer_size, 'num_layers': num_layers, 'num_monomers': num_monomers}

agg_models_path = os.path.expanduser('~/runs/Agg_models/NEW_RUNS/single_run')
folders = glob.glob(os.path.join(agg_models_path, "MONO_*_*"))

monomer_sizes_list = sorted(set([int(folder.split('_')[-2]) for folder in folders]))
num_layers_list = sorted(set([int(folder.split('_')[-1]) for folder in folders]))

# Create a grid of subplots
fig, axes = plt.subplots(len(num_layers_list), len(monomer_sizes_list), figsize=(20, 20), sharex=True, sharey=True)

# Ensure axes are 2D for easier indexing
if len(num_layers_list) == 1:
    axes = axes[np.newaxis, :]
if len(monomer_sizes_list) == 1:
    axes = axes[:, np.newaxis]

avg_x_values = {}
avg_y_values = {}

# Compute average
for folder in folders:
    plot_data = plot_results(folder)
    if plot_data is None:
        continue

    i = num_layers_list.index(plot_data['num_layers'])
    j = monomer_sizes_list.index(plot_data['monomer_size'])

    # Store x and y values for averaging later
    avg_x_values[(i, j)] = avg_x_values.get((i, j), []) + list(plot_data['x'])
    avg_y_values[(i, j)] = avg_y_values.get((i, j), []) + list(plot_data['y'])

for key in avg_x_values:
    avg_x_values[key] = np.mean(avg_x_values[key])
    avg_y_values[key] = np.mean(avg_y_values[key])

# Plotting
for folder in folders:
    plot_data = plot_results(folder)
    if plot_data is None:
        continue

    i = num_layers_list.index(plot_data['num_layers'])
    j = monomer_sizes_list.index(plot_data['monomer_size'])

    # Plotting the actual data and the average values
    axes[i, j].plot(plot_data['x'], plot_data['y'], '.', label="Aggregate")
    axes[i, j].plot(avg_x_values[(i, j)], avg_y_values[(i, j)], 'bo', label="Average Aggregate")

    axes[i, j].grid(True)  # Add grid for better readability
    axes[i, j].set_title(f"Size: {plot_data['monomer_size']}, Layers: {plot_data['num_layers']}")
    axes[i, j].legend()  # To show the legend indicating which plot is Aggregate and which is the average

# Setting labels
for j, monomer_size in enumerate(monomer_sizes_list):
    axes[-1, j].set_xlabel(f's.a.')
for i, num_layers in enumerate(num_layers_list):
    axes[i, 0].set_ylabel(f'Polarization')

plt.suptitle("Effect of Number of Layers and Monomer Size on Polarisation", fontsize=24)
plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2, w_pad=2)
plt.savefig('grid_plot_layers_vs_size_single_mon.png', format='png', bbox_inches='tight')
#plt.show()

