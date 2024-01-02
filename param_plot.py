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

# Extract unique values for each parameter
#num_monomers_list = sorted(set([int(f.split('_')[-2]) for f in folders]))
#num_monomers_list = sorted(set([int(f.split('_')[-2]) for f in folders if f.split('_')[-2].isdigit()]))
#num_layers_list = sorted(set([int(f.split('_')[-1]) for f in folders]))
#monomer_sizes_list = sorted(set([read_monomer_size(f) for f in folders if read_monomer_size(f) is not None]))
monomer_sizes_list = [140, 160, 180, 200, 220, 240, 260]
num_layers_list = [2, 3, 4, 5, 6, 7]


#folders = [f for f in folders if int(f.split('_')[-2]) in monomer_sizes_list and int(f.split('_')[-1]) in num_layers_list]
folders = [
    f for f in folders 
    if f.split('_')[-2].isdigit() and int(f.split('_')[-2]) in monomer_sizes_list 
    and f.split('_')[-1].isdigit() and int(f.split('_')[-1]) in num_layers_list
]

# Create the grid of subplots
fig, axes = plt.subplots(len(num_layers_list), len(monomer_sizes_list), figsize=(20, 20), sharex='col', sharey='row')

# Go through all the folders and create the subplots
#for folder in folders:
#    plot_data = plot_results(folder)
#    if plot_data is None:
#        continue

for folder in folders:
    monomer_size = read_monomer_size(folder)
    if monomer_size is not None:
        print(f"Folder: {folder}, Monomer Size: {monomer_size}")


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
plt.savefig('grid_plot_layers_vs_size_2.pdf', format='pdf', bbox_inches='tight')
plt.show()


# Define the parameter pairs for creating the grids of subplots
#parameter_pairs = [('Size vs Layers', monomer_sizes, num_layers_list),
#                   ('Size vs Monomers', monomer_sizes, num_monomers_list),
#                   ('Layers vs Monomers', num_layers_list, num_monomers_list)]

#for pair_name, param1_list, param2_list in parameter_pairs:
#    fig, axes = plt.subplots(len(param1_list), len(param2_list), figsize=(20, 20), sharex='col', sharey='row')
#
#    for folder in folders:
#        plot_data = plot_results(folder)
#        if plot_data is None:
#            continue
#        i = param1_list.index(plot_data['monomer_size'])
#        j = param2_list.index(plot_data['num_layers'] if pair_name == 'Size vs Layers' else plot_data['num_monomers'])
#
#        axes[i, j].plot(plot_data['x'], plot_data['y'], '.')
#        axes[i, j].set_title(f'Monomer size: {plot_data["monomer_size"]}, Layers: {plot_data["num_layers"]}, Monomers: {plot_data["num_monomers"]}')
#        axes[0, j].set_xlabel(pair_name.split()[2])
#        axes[i, 0].set_ylabel(pair_name.split()[0])
#
#    plt.suptitle(pair_name, fontsize=24)
#    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
#    plt.savefig(f'{pair_name}_grid_plot.eps', format='eps', bbox_inches='tight')
#    plt.show()


#for pair_name, param1_list, param2_list in parameter_pairs:
#    fig, axes = plt.subplots(len(param1_list), len(param2_list), figsize=(24, 24), sharex='col', sharey='row')
#
#    for folder in folders:
#        plot_data = plot_results(folder)
#
#        if plot_data is None:
#            continue
#
#        i = param1_list.index(plot_data['monomer_size'])
#        j = param2_list.index(plot_data['num_layers'] if pair_name == 'Size vs Layers' else plot_data['num_monomers'])

#        axes[i, j].plot(plot_data['x'], plot_data['y'], '.')
#        axes[0, j].set_xlabel(f"{pair_name.split()[2]}: {param2_list[j]}", fontsize=10)
#        axes[i, 0].set_ylabel(f"{pair_name.split()[0]}: {param1_list[i]}", fontsize=10)
#        axes[i, j].tick_params(axis='both', labelsize=8)

 #   plt.suptitle(pair_name, fontsize=24)
 #   plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2, w_pad=2)
 #   plt.savefig(f'{pair_name}_grid_plot.eps', format='eps', bbox_inches='tight')
 #   plt.show()
