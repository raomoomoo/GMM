import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import numpy as np

# Your existing functions read_monomer_size and plot_results
# Modify the plot_results function to return the data

def plot_results(directory_path):
    num_monomers, num_layers = map(int, directory_path.split('_')[-2:])
    monomer_size = read_monomer_size(directory_path)
    file_path = os.path.join(directory_path, f'gmm01s_{num_monomers:05d}.out')
    try:
        data = pd.read_csv(file_path, skiprows=18, delim_whitespace=True, header=0)
    except FileNotFoundError:
        print(f"Skipping {directory_path}: {file_path} not found.")
        return

#    data = pd.read_csv(file_path, skiprows=18, delim_whitespace=True, header=0)

    x = data['s.a.']
    y = data['pol.']

    return {'x': x, 'y': y}

def read_monomer_size(directory_path):
    k_file_path = glob.glob(os.path.join(directory_path, "*.k"))
    if not k_file_path:
        return None
    k_file_path = k_file_path[0]
    with open(k_file_path, 'r') as k_file:
        for line in k_file:
            if len(line.strip().split()) == 6:
                monomer_size = float(line.strip().split()[3])
                break
    return monomer_size

# Rest of the code

agg_models_path = os.path.expanduser('~/runs/Agg_models/NEW_RUNS/')
folders = glob.glob(os.path.join(agg_models_path, "AGG_*_*"))

# Extract unique values for each parameter
num_monomers_list = sorted(set([int(f.split('_')[-2]) for f in folders]))
num_layers_list = sorted(set([int(f.split('_')[-1]) for f in folders]))
monomer_sizes = sorted(set([read_monomer_size(f) for f in folders if read_monomer_size(f) is not None]))

# Create a grid of plots for each combination of parameters
fig, axes = plt.subplots(len(monomer_sizes), len(num_layers_list), figsize=(20, 20))
mer_sizes = sorted(set([read_monomer_size(f) for f in folders if read_monomer_size(f) is not None]))

#for i, monomer_size in enumerate(monomer_sizes):
#    for j, num_layers in enumerate(num_layers_list):
#        for num_monomers in num_monomers_list:
#            directory = [f for f in folders if int(f.split('_')[-2]) == num_monomers and int(f.split('_')[-1]) == num_layers]
#            if directory:
#                data = plot_results(directory[0])
#                axes[i, j].plot(data['x'], data['y'], '.', label=f'GMM {num_monomers} monomers')
#                axes[i, j].set_title(f'Size: {monomer_size}, Layers: {num_layers}')
#                axes[i, j].set_xlabel('Scatter angle')
#                axes[i, j].set_ylabel('$-F_{12}/F_{11}$')
#                axes[i, j].legend()

for i, monomer_size in enumerate(monomer_sizes):
    for j, num_layers in enumerate(num_layers_list):
        for num_monomers in num_monomers_list:
            directory = [f for f in folders if int(f.split('_')[-2]) == num_monomers and int(f.split('_')[-1]) == num_layers]
            if directory:
                data = plot_results(directory[0])
                if data is not None:
                    if len(monomer_sizes) == 1 or len(num_layers_list) == 1:
                        ax = axes[max(i, j)]
                    else:
                        ax = axes[i, j]
                    ax.plot(data['x'], data['y'], '.', label=f'GMM {num_monomers} monomers')
                    ax.set_title(f'Size: {monomer_size}, Layers: {num_layers}')
                    ax.set_xlabel('Scatter angle')
                    ax.set_ylabel('$-F_{12}/F_{11}$')
                    ax.legend()


plt.tight_layout()
plt.savefig('parameter_combination_plots.eps', format='eps', bbox_inches='tight')
plt.show()

