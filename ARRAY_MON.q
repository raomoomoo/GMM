#!/bin/bash
#SBATCH --array=0-13
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
#SBATCH --partition=milan
#SBATCH --mem=4G
#SBATCH --mail-user=rmor0013@student.monash.edu

# Load modules
#module load gcc/10.3.0 openmpi/4.1.1
#module load pandas/1.2.4-scipy-bundle-2021.05
#module load matplotlib/3.4.2
#module load numba/0.51.2-python-3.8.5
#module load scipy/1.6.0-python-3.8.5


module load gcc/11.3.0 openmpi/4.1.4
module load pandas/1.4.2-scipy-bundle-2022.05
module load matplotlib/3.5.2
module load numba/0.56.4

GMM_PATH=~/GMM
export GMM_PATH

# Parent directory to store individual runs
PARENT_DIR=~/runs/Agg_models/NEW_RUNS/single_run

# Parameter grid
monomer_sizes=(25 50 75 100 125 100 140 160 180 200 220 240 260)
#uncomment for equivalent volume
#num_layers_list=(2 3 4 5 6 7)

# Calculate the current parameter combination
monomer_size_idx=$(( SLURM_ARRAY_TASK_ID / 13 ))
#num_layers_idx=$(( SLURM_ARRAY_TASK_ID % 6 ))

monomer_size=${monomer_sizes[$monomer_size_idx]}
#num_layers=${num_layers_list[$num_layers_idx]}

# Generate .k file
#python single_mon.py "$monomer_size" "$num_layers"
python3 single_run.py "$monomer_size" 

# Define directory and filenames
#AGG_NAME="monomer_NIR_${monomer_size}_${num_layers}.k"
AGG_NAME="monomer_NIR_${monomer_size}.k"
AGG_PATH="$HOME/runs/Agg_models/input/single_run/$AGG_NAME"
NEW_RUN_DIR="$PARENT_DIR"

# Create a new directory for the current model
#NEW_DIR="MONO_${monomer_size}_${num_layers}"
NEW_DIR="MONO_${monomer_size}"
MODEL_DIR="$NEW_RUN_DIR/$NEW_DIR"

# Skip if directory already exists
#if [ -d "$MODEL_DIR" ]; then
#  echo "WARNING: Directory $MODEL_DIR already exists. Skipping this run."
#  exit
#fi

mkdir -p "$MODEL_DIR"

# Copy required files
# Copy required files
cp "$AGG_PATH" "$MODEL_DIR" || { echo "Error copying .k file"; exit 1; }
cp "$GMM_PATH/gmm01f.in" "$MODEL_DIR" || { echo "Error copying gmm01f.in"; exit 1; }
cp "$GMM_PATH/gmm01s.f90" "$MODEL_DIR" || { echo "Error copying gmm01s.f90"; exit 1; }
cp "$GMM_PATH/gmm01f.par" "$MODEL_DIR" || { echo "Error copying gmm01f.par"; exit 1; }

#cp "$AGG_PATH" "$MODEL_DIR"
#cp "$GMM_PATH/gmm01f.in" "$MODEL_DIR"
#cp "$GMM_PATH/gmm01s.f90" "$MODEL_DIR"
#cp "$GMM_PATH/gmm01f.par" "$MODEL_DIR"

# Create run.q file
RUNQ_FILE="$MODEL_DIR/run.q"
cat > "$RUNQ_FILE" << EOL
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=GMM_${monomer_size}
#SBATCH --output=MONO_${monomer_size}.in.qout
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --time=0-168:00:00
#SBATCH --partition=milan
#SBATCH --mail-user=rmor0013@student.monash.edu
#SBATCH --mem=8G

echo "HOSTNAME = \$HOSTNAME"
echo "HOSTTYPE = \$HOSTTYPE"
echo Time is \`date\`
echo Directory is \`pwd\`

ulimit -s unlimited
export OMP_SCHEDULE="dynamic"
export OMP_NUM_THREADS=4
export OMP_STACKSIZE=1024m

echo "starting GMM run..."
echo "writing output to \$outfile"

./gmm_${monomer_size}.out
EOL

# Make the run.q file executable
chmod +x "$RUNQ_FILE"

# Update Fortran input file
sed -i "1s|.*|$(basename $AGG_PATH)|" "$MODEL_DIR/gmm01f.in"

# Compile Fortran code
ifort -o "gmm_${monomer_size}.out" "$MODEL_DIR/gmm01s.f90"

# Move compiled output
mv "gmm_${monomer_size}.out" "$MODEL_DIR"

# Submit the run.q file
cd "$MODEL_DIR"
sbatch run.q

