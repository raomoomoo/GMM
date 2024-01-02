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
    parts = directory_path.split('_')
    try:
        num_monomers = int(parts[-2])
        num_layers = int(parts[-1])
    except ValueError:
        print(f"Skipping {directory_path}: Unexpected directory name format.")
        return
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

# Just filter based on number of layers
folders = [
    f for f in folders
    if f.split('_')[-1].isdigit() and int(f.split('_')[-1]) in num_layers_list
]

fig, axes = plt.subplots(len(num_layers_list), len(monomer_sizes_list), figsize=(25, 25), sharex='col', sharey='row')


# Set style for better readability
plt.style.use('seaborn-whitegrid')

# Define a color for the plots
plot_color = 'navy'  # A color that is generally well-received and prints well in black and white

# Adjusting marker style and size
marker_style = 'o'  # This could be any other marker style as per your preference
marker_size = 3     # Adjust size to ensure the plot is not too cluttered

# ... [Your previous code]

# Set style for better readability
plt.style.use('seaborn-whitegrid')

# Define a color for the plots
plot_color = 'navy'  # A color that is generally well-received and prints well in black and white

# Adjusting marker style and size
marker_style = 'o'  # This could be any other marker style as per your preference
marker_size = 3     # Adjust size to ensure the plot is not too cluttered

# Loop through your data to create subplots
for folder in folders:
    plot_data = plot_results(folder)

    if plot_data is None:
        continue
#    # Check if the extracted monomer_size is in the list of sizes you want to plot
    if plot_data['monomer_size'] not in monomer_sizes_list:
        continue
    i = num_layers_list.index(plot_data['num_layers'])
    j = monomer_sizes_list.index(plot_data['monomer_size'])
    ax = axes[i, j]

    # Using a line with marker style here for the plot
    ax.plot(plot_data['x'], plot_data['y'], marker=marker_style, markersize=marker_size, linestyle='-', color=plot_color)

    # Optional: Set individual titles for each subplot if needed
    # ax.set_title(f'Layers: {plot_data['num_layers']} Size: {plot_data['monomer_size']}')

# Adjusting the layout and saving the figure
fig.subplots_adjust(hspace=0.4, wspace=0.4)  # Adjust horizontal and vertical spacing
for ax in axes.flat:
    ax.label_outer()  # Hide x labels and tick labels for top plots and y ticks for right plots.

# Enhancing the x and y axis labels font size
plt.setp(axes[-1, :], xlabel='Wavelength (nm)', fontsize=24)
plt.setp(axes[:, 0], ylabel='Polarization', fontsize=24)

# Set a comprehensive figure title and adjust the font size
plt.suptitle("Effect of Number of Layers and Size of Monomers on Polarisation", fontsize=48)

# Save and show the figure
plt.savefig('grid_plot_layers_vs_size_mon_comp_improved.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.show()

#for folder in folders:
#    plot_data = plot_results(folder)
#    
#    if plot_data is None:
#        continue
#
#    # Check if the extracted monomer_size is in the list of sizes you want to plot
#    if plot_data['monomer_size'] not in monomer_sizes_list:
#        continue
#
#    # The rest of the loop stays unchanged
#    i = num_layers_list.index(plot_data['num_layers'])
#    j = monomer_sizes_list.index(plot_data['monomer_size'])
#    axes[i, j].plot(plot_data['x']
#    ax.plot(plot_data['x'], plot_data['y'], marker=marker_style, markersize=marker_size, linestyle='-', color=plot_color)
#    
# The rest of the script stays unchanged
#for j, monomer_size in enumerate(monomer_sizes_list):
#    axes[-1, j].set_xlabel(f'Monomer size: {monomer_size}')
#
#for i, num_layers in enumerate(num_layers_list):
#    axes[i, 0].set_ylabel(f'Layers: {num_layers}')
#
#plt.suptitle("Effect of Number of Layers and size of Monomers on Polarisation", fontsize=48)
#plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2, w_pad=2)
#plt.savefig('grid_plot_layers_vs_size_mon_comp.pdf', format='pdf', bbox_inches='tight', dpi=300)
#plt.show()

