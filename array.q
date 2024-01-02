#!/bin/bash
#SBATCH --job-name=AGG_RUN
#SBATCH --ntasks=1
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=rmor0013@student.monash.edu
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
#SBATCH --mem=8G
#SBATCH --array=0-839

# Load required modules
#module load gcc/10.3.0 openmpi/4.1.1
#module load numpy/1.18.3-python-3.8.5
#module load numba/0.51.2-python-3.8.5
#module load scipy/1.6.0-python-3.8.5
#module load pandas/1.0.5-python-3.8.5

module load gcc/11.3.0 openmpi/4.1.4
module pandas/1.4.2-scipy-bundle-2022.05
module load matplotlib/3.5.2
module load numba/0.56.4
module load pandas/1.0.5-python-3.8.5

# Define the parameter grid
monomer_sizes=(25 50 75 100 125 140 160 180)
num_layers_list=(1 2 3 4 5 6 7) 
num_monomers_list=(100 150 200 250 300 350 400 450 500 550 600 650 700 750 800)

#wavelengths_list=(500 550 600 650 700 750 800 850 900 950 1000)  # Add more wavelengths as needed
# Calculate the current parameter combination
monomer_size_idx=$(( ${SLURM_ARRAY_TASK_ID} / (7 * 15) ))
num_layers_idx=$(( (${SLURM_ARRAY_TASK_ID} % (7 * 15)) / 15 ))
num_monomers_idx=$(( ${SLURM_ARRAY_TASK_ID} % 15 ))

monomer_size=${monomer_sizes[$monomer_size_idx]}
num_layers=${num_layers_list[$num_layers_idx]}
num_monomers=${num_monomers_list[$num_monomers_idx]}

echo "Running with monomer_size=$monomer_size, num_layers=$num_layers, num_monomers=$num_monomers"

# Run the original AGG_RUN.sh script with the current parameter combination
./AGG_RUN.sh "$monomer_size" "$num_layers" "$num_monomers"

# Wait for the job to finish (optional, adjust the sleep time as needed)
sleep 10s

# Call the plot_results function with the current parameter combination
python -c "from plot_results import plot_results; plot_results($monomer_size, $num_layers, $num_monomers)"

