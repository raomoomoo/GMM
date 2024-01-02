import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def read_monomer_size(directory_path):
    k_file_path = glob.glob(os.path.join(directory_path, "*.k"))
    if not k_file_path:
        return None
    k_file_path = k_file_path[0]
    with open(k_file_path, 'r') as k_file:
        lines = k_file.readlines()
        monomer_size_line = lines[3]
        monomer_size = float(monomer_size_line.strip().split()[3])
    return monomer_size

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

agg_models_path = os.path.expanduser('~/runs/Agg_models/NEW_RUNS/')
folders = glob.glob(os.path.join(agg_models_path, "AGG_*_*"))

# Define desired monomer sizes in nanometers and convert to microns
monomer_sizes_nm = [140, 160, 180, 200, 220, 240, 260]
monomer_sizes_list = [size / 1000 for size in monomer_sizes_nm]

# Filter the folders to only include those with the desired monomer sizes
folders = [f for f in folders if read_monomer_size(f) in monomer_sizes_list]

# Extract unique values for each parameter
num_monomers_list = sorted(set([int(f.split('_')[-2]) for f in folders if f.split('_')[-2].isdigit()]))
num_layers_list = sorted(set([int(f.split('_')[-1]) for f in folders]))

# Create the grid of subplots
fig, axes = plt.subplots(len(num_layers_list), len(monomer_sizes_list), figsize=(20, 20), sharex='col', sharey='row')

# Go through all the folders and create the subplots
for folder in folders:
    plot_data = plot_results(folder)
    if plot_data is None:
        continue

    i = num_layers_list.index(plot_data['num_layers'])
    j = monomer_sizes_list.index(plot_data['monomer_size'])

    axes[i, j].plot(plot_data['x'], plot_data['y'], '.')

# Set the x and y axis labels
for j, monomer_size in enumerate(monomer_sizes_list):
    axes[-1, j].set_xlabel(f'Monoomer size: {monomer_size}')

for i, num_layers in enumerate(num_layers_list):
    axes[i, 0].set_ylabel(f'Layers: {num_layers}')

# Set the overall title
plt.suptitle("Effect of Number of Layers and size of Monomers on Polarisation", fontsize=24)

# Save and show the plot
plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2, w_pad=2)
plt.savefig('grid_plot_layers_vs_size_3.pdf', format='pdf', bbox_inches='tight')
plt.show()

