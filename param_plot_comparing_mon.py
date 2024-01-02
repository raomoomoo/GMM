import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

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

agg_models_path = os.path.expanduser('~/runs/Agg_models/NEW_RUNS/')
folders = glob.glob(os.path.join(agg_models_path, "AGG_*_*"))

monomer_sizes_list = [140, 160, 180, 200, 220, 240, 260]
num_layers_list = [2, 3, 4, 5, 6, 7]

folders = [
    f for f in folders
    if f.split('_')[-2].isdigit() and int(f.split('_')[-2]) in monomer_sizes_list
    and f.split('_')[-1].isdigit() and int(f.split('_')[-1]) in num_layers_list
]

#j = monomer_sizes_list.index(plot_data['monomer_size'])
#print("Searching for:", plot_data['monomer_size'])
print("In list:", monomer_sizes_list)
#j = monomer_sizes_list.index(plot_data['monomer_size'])
fig, axes = plt.subplots(len(num_layers_list), len(monomer_sizes_list), figsize=(20, 20), sharex='col', sharey='row')

for folder in folders:
    print(f"Filtered folder: {folder}")
    plot_data = plot_results(folder)
    if plot_data is None:
        continue
    print(f"Processing folder: {folder}")
    print(f"Extracted monomer size: {plot_data['monomer_size']}")
    # The rest of the loop stays unchanged
    i = num_layers_list.index(plot_data['num_layers'])
    j = monomer_sizes_list.index(plot_data['monomer_size'])

    axes[i, j].plot(plot_data['x'], plot_data['y'], '.')

# The rest of the script stays unchanged
for j, monomer_size in enumerate(monomer_sizes_list):
    axes[-1, j].set_xlabel(f'Monoomer size: {monomer_size}')

for i, num_layers in enumerate(num_layers_list):
    axes[i, 0].set_ylabel(f'Layers: {num_layers}')

plt.suptitle("Effect of Number of Layers and size of Monomers on Polarisation", fontsize=24)
plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2, w_pad=2)
plt.savefig('grid_plot_layers_vs_size_mon_comp.pdf', format='pdf', bbox_inches='tight')
plt.show()
a
