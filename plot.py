import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

#def read_monomer_size(directory_path):
#    k_file_path = glob.glob(os.path.join(directory_path, "*.k"))
#    if not k_file_path:
#        return None
#    k_file_path = k_file_path[0]
#    with open(k_file_path, 'r') as k_file:
#        for line in k_file:
#            if len(line.strip().split()) == 6:
#                monomer_size = float(line.strip().split()[3])
#                break
#    return monomer_size


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

#    data = pd.read_csv(file_path, skiprows=18, delim_whitespace=True, header=0)

    try:
        data = pd.read_csv(file_path, skiprows=18, delim_whitespace=True, header=0)
    except FileNotFoundError:
        print(f"Skipping {directory_path}: {file_path} not found.")
        return

    x = data['s.a.']
    y = data['pol.']

    plt.figure(figsize=(6,6))
    plt.rcParams.update({'font.size': 20})
 #   plt.plot(x, y, '.', label='GMM CODE')
    plt.xlabel('Scatter angle')
    plt.ylabel('$-F_{12}/F_{11}$')
   # plt.legend()

    plot_filename = os.path.join(directory_path, f'polarisation_{monomer_size}_{num_layers}_{num_monomers}.eps')

    plt.savefig(plot_filename, format='eps', bbox_inches='tight')
    plt.close()
#    plt.show()

agg_models_path = os.path.expanduser('~/runs/Agg_models/NEW_RUNS/')
folders = glob.glob(os.path.join(agg_models_path, "AGG_*_*"))

# Create a plot for each directory
for directory in folders:
    plot_results(directory)

