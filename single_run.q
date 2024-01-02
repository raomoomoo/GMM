#!/bin/bash
#SBATCH --mail-user=rmor0013@student.monash.edu
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
#SBATCH --mem=8G
#SBATCH --array=0-41


# Load required modules
#module load numpy/1.18.3-python-3.8.5
#module load numba/0.51.2-python-3.8.5
#module load scipy/1.6.0-python-3.8.5
#module load pandas/1.0.5-python-3.8.5

module load gcc/11.3.0 openmpi/4.1.4
module load pandas/1.4.2-scipy-bundle-2022.05
module load matplotlib/3.5.2
module load numba/0.56.4
#module load pandas/1.0.5-python-3.8.5

# Define the parameter grid
monomer_sizes=(140 160 180 200 220 240 260)
num_layers_list=(2 3 4 5 6 7)

# Calculate the current parameter combination
#monomer_size_idx=$(( ${SLURM_ARRAY_TASK_ID} / (6 * 3) ))
#num_layers_idx=$(( (${SLURM_ARRAY_TASK_ID} % (6 * 3)) / 3 ))

monomer_size_idx=$(( ${SLURM_ARRAY_TASK_ID} / (6 * 7) ))
num_layers_idx=$(( ${SLURM_ARRAY_TASK_ID} % 7 ))


monomer_size=${monomer_sizes[$monomer_size_idx]}
num_layers=${num_layers_list[$num_layers_idx]}

echo "Running with monomer_size=$monomer_size, num_layers=$num_layers"

# Run the Python script with the current parameter combination
python single_mon_grid.py "$monomer_size" "$num_layers"

