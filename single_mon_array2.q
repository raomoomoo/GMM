#!/bin/bash
#SBATCH --array=0-41
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
#SBATCH --mem=8G
#SBATCH --mail-user=rmor0013@student.monash.edu

# Load required modules
module load matplotlib/3.3.2-python-3.8.5
module load numpy/1.18.3-python-3.8.5
module load numba/0.51.2-python-3.8.5
module load scipy/1.6.0-python-3.8.5
module load pandas/1.0.5-python-3.8.5

# Parent directory to store individual runs
PARENT_DIR=~/runs/Agg_models/NEW_RUNS/single_run

# Define the parameter grid
monomer_sizes=(25 50 75 100 125 150 175)
num_layers_list=(2 3 4 5 6 7)

# Calculate the current parameter combination
monomer_size_idx=$(( SLURM_ARRAY_TASK_ID / 6 ))
num_layers_idx=$(( SLURM_ARRAY_TASK_ID % 6 ))

monomer_size=${monomer_sizes[$monomer_size_idx]}
num_layers=${num_layers_list[$num_layers_idx]}

echo "Running with monomer_size=$monomer_size, num_layers=$num_layers"

# Unique directory for this parameter set
MODEL_DIR="$PARENT_DIR/AGG_${monomer_size}_${num_layers}"

# Create the directory
mkdir -p "$MODEL_DIR"

# Create the run.q file
RUNQ_FILE="$MODEL_DIR/run.q"
cat > "$RUNQ_FILE" << "EOL"
#!/bin/bash
#SBATCH --nodes=1 --ntasks=4
#SBATCH --cpus-per-task=1
# Job name including the parameter set for easy identification
#SBATCH --job-name=GMM_${monomer_size}_${num_layers}
#SBATCH --output=AGG_${num_layers}.in.qout
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --time=0-168:00:00
#SBATCH --mail-user=rmor0013@student.monash.edu
#SBATCH --mem=8G

# Load required modules
module load numpy/1.18.3-python-3.8.5
module load numba/0.51.2-python-3.8.5
module load scipy/1.6.0-python-3.8.5
module load pandas/1.0.5-python-3.8.5

echo "Running with monomer_size=$monomer_size, num_layers=$num_layers"

# Execute your Python script here
python single_mon.py "$monomer_size" "$num_layers"
EOL

# Make the run.q file executable
chmod +x "$RUNQ_FILE"

# Submit the nested job
cd "$MODEL_DIR"
sbatch run.q

