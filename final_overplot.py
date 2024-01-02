import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import re
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

def get_single_monomer_data(monomer_size, num_layers):
    directory_path = os.path.expanduser(f'~/runs/Agg_models/NEW_RUNS/single_run/MONO_{monomer_size}_{num_layers}')
    file_paths = glob.glob(os.path.join(directory_path, 'gmm01s_*.out'))

    if not file_paths:
        return None, None

    data = pd.read_csv(file_paths[0], skiprows=18, delim_whitespace=True, header=0)

    return data['s.a.'], data['pol.']

def get_average_aggregate_data(folders, monomer_size, num_layers):
    pol_sums = None
    count = 0
    x_data = None

    for folder in folders:
        plot_data = plot_results(folder)
        if plot_data and plot_data['monomer_size'] == monomer_size and plot_data['num_layers'] == num_layers:
            if pol_sums is None:
                pol_sums = plot_data['y'].to_numpy()
                x_data = plot_data['x']
            else:
                pol_sums += plot_data['y'].to_numpy()
            count += 1

    if count > 0 and x_data is not None:
        avg_pol = pol_sums / count
        return x_data, avg_pol
    else:
        return None, None


agg_models_path = os.path.expanduser('~/runs/Agg_models/NEW_RUNS/')
folders = glob.glob(os.path.join(agg_models_path, "AGG_*_*"))

monomer_sizes_list = [140, 160, 180, 200, 220, 240, 260]
num_layers_list = [2, 3, 4, 5, 6, 7]

pattern = re.compile(r'AGG_\d+_(\d+)$')
valid_folders = [f for f in folders if pattern.match(os.path.basename(f)) and int(pattern.match(os.path.basename(f)).group(1)) in num_layers_list]

fig, axes = plt.subplots(len(num_layers_list), len(monomer_sizes_list), figsize=(20, 20), sharex='col', sharey='row')

for monomer_size in monomer_sizes_list:
    for num_layers in num_layers_list:
        avg_agg_x, avg_agg_y = get_average_aggregate_data(valid_folders, monomer_size, num_layers)
        single_monomer_x, single_monomer_y = get_single_monomer_data(monomer_size, num_layers)

        i = num_layers_list.index(num_layers)
        j = monomer_sizes_list.index(monomer_size)

        if avg_agg_x is not None and avg_agg_y is not None:
            axes[i, j].plot(avg_agg_x, avg_agg_y, '-o', label='Average Aggregate')
        if single_monomer_x is not None and single_monomer_y is not None:
            axes[i, j].plot(single_monomer_x, single_monomer_y, '--', label='Single Monomer')
        axes[i, j].legend()


for j, monomer_size in enumerate(monomer_sizes_list):
    axes[-1, j].set_xlabel(f'Monomer size: {monomer_size}')
for i, num_layers in enumerate(num_layers_list):
    axes[i, 0].set_ylabel(f'Layers: {num_layers}')
    
plt.suptitle("Effect of Number of Layers and size of Monomers on Polarisation", fontsize=24)
#plt.tight_layout(rect=[0, 0.03, 1, 0.97])
plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2, w_pad=2)
plt.savefig('grid_plot_layers_vs_size_mon_overplot.pdf', format='pdf', bbox_inches='tight')
#plt.show()

