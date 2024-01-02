#import pandas as pd
#import matplotlib.pyplot as plt
#module load gcc/10.3.0 openmpi/4.1.1
#module load pandas/1.2.4-scipy-bundle-2021.05
#module load matplotlib/3.4.2

import pandas as pd
import matplotlib.pyplot as plt

def plot_results(monomer_size, num_layers, num_monomers):
    # Define the file path based on the parameter values
    file_path = f'./Agg_models/NEW_RUNS/AGG_{monomer_size}_{num_layers}/gmm01s_{num_monomers:05d}.out'

    # Load the data from the file
    data = pd.read_csv(file_path, skiprows=18, delim_whitespace=True, header=0)

    # Plot the results
    x = data['s.a.']
    y = data['pol.']

    plt.figure(figsize=(6,6))
    plt.rcParams.update({'font.size': 20})
    plt.plot(x, y, '.', label='GMM CODE')
    plt.xlabel('Scatter angle')
    plt.ylabel('$-F_{12}/F_{11}$')
    plt.legend()

    # Save the plot with a filename that includes the parameter values
    plot_filename = f'polarisation_{monomer_size}_{num_layers}_{num_monomers}.eps'
    plt.savefig(plot_filename, format='eps', bbox_inches='tight')

    # Show the plot (optional)
    plt.show()

