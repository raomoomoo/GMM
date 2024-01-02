import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import re

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
    monomer_size = read_monomer_size(directory_path)
    if monomer_size is None or monomer_size not in monomer_sizes_list:
        print(f"Skipping {directory_path}: Invalid monomer size.")
        return

    num_monomers, num_layers = map(int, directory_path.split('_')[-2:])

    file_path = os.path.join(directory_path, f'gmm01s_{num_monomers:05d}.out')

    try:
        data = pd.read_csv(file_path, skiprows=18, delim_whitespace=True, header=0)
    except FileNotFoundError:
        print(f"Skipping {directory_path}: {file_path} not found.")
        return

    x = data['s.a.']
    y = data['pol.']

    return {'x': x, 'y': y, 'monomer_size': monomer_size, 'num_layers': num_layers, 'num_monomers': num_monomers}

#def get_single_monomer_data(equivalent_size, num_layers):
#    folder_path = os.path.join(single_monomer_path, f"MONO_{equivalent_size}_{num_layers}")
#    file_path = os.path.join(folder_path, "*.out")
#    
#    file_paths = glob.glob(file_path)
#    
#    if not file_paths:
#        return None, None
#    data = pd.read_csv(file_paths[0], skiprows=18, delim_whitespace=True, header=0)
#    return data['s.a.'], data['pol.']

def get_single_monomer_data(equivalent_size, num_layers):
    folder_path = os.path.join(single_monomer_path, f"MONO_{equivalent_size}_{num_layers}")
    file_path = os.path.join(folder_path, "gmm01s_00000.out")

    file_paths = glob.glob(file_path)

    if not file_paths:
        return None, None

    try:
        data = pd.read_csv(file_paths[0], skiprows=18, delim_whitespace=True, header=0)
        return data['s.a.'], data['pol.']
    except pd.errors.EmptyDataError:
        print(f"Warning: {file_paths[0]} appears to be empty or lacks columns. Skipping...")
        return None, None
    except pd.errors.ParserError:
        print(f"Warning: Error parsing {file_paths[0]}. Possible column mismatch. Skipping...")
        return None, None

agg_models_path = os.path.expanduser('~/runs/Agg_models/NEW_RUNS/')
single_monomer_path = os.path.expanduser('~/runs/Agg_models/NEW_RUNS/single_run')
folders = glob.glob(os.path.join(agg_models_path, "AGG_*_*"))

monomer_sizes_list = [140, 160, 180, 200, 220, 240, 260]
num_layers_list = [2, 3, 4, 5, 6, 7]

pattern = re.compile(r'AGG_\d+_(\d+)$')
valid_folders = [f for f in folders if pattern.match(os.path.basename(f)) and int(pattern.match(os.path.basename(f)).group(1)) in num_layers_list]

fig, axes = plt.subplots(len(num_layers_list), len(monomer_sizes_list), figsize=(20, 20), sharex='col', sharey='row')

for folder in valid_folders:
    plot_data = plot_results(folder)
    if plot_data is None:
        continue

    i = num_layers_list.index(plot_data['num_layers'])
    j = monomer_sizes_list.index(plot_data['monomer_size'])

    # Plot the aggregate data
    axes[i, j].plot(plot_data['x'], plot_data['y'], '.', label="Aggregate")

    # Plot the "equivalent" single monomer data
    single_monomer_x, single_monomer_y = get_single_monomer_data(plot_data['monomer_size'], plot_data['num_layers'])
    if single_monomer_x is not None and single_monomer_y is not None:
        axes[i, j].plot(single_monomer_x, single_monomer_y, '-', label="Single Monomer")
    axes[i, j].legend()

# Labeling and formatting
for j, monomer_size in enumerate(monomer_sizes_list):
    axes[-1, j].set_xlabel(f'Monomer size: {monomer_size}')
for i, num_layers in enumerate(num_layers_list):
    axes[i, 0].set_ylabel(f'Layers: {num_layers}')
plt.suptitle("Effect of Number of Layers and size of Monomers on Polarisation", fontsize=24)
plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2, w_pad=2)
plt.savefig('grid_plot_layers_vs_size_mon_comp.pdf', format='pdf', bbox_inches='tight')
plt.show()

