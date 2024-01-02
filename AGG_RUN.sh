#!/bin/bash
#module load matplotlib/3.3.2-python-3.8.5
#module load numpy/1.18.3-python-3.8.5
#module load numba/0.51.2-python-3.8.5
#module load scipy/1.6.0-python-3.8.5 

#module load gcc/10.3.0 openmpi/4.1.1
module load gcc/11.3.0 openmpi/4.1.4
module pandas/1.4.2-scipy-bundle-2022.05
module load matplotlib/3.5.2
module load numba/0.56.4


GMM_PATH=~/GMM
export GMM_PATH
#this is commented out for the grid 
#name the file 
#echo "name the file"
#read -p "Enter the number of monomers: " n_m
#read -p "Enter the number of layers: " r_m

n_m=$1
r_m=$2
r_m=${r_m:-10}
num_monomers=$3

aggregate_output=$(python3 ~/runs/Agg_models/GMM-Aggregates/Multilayer-Agg-code_fast.py "$n_m" "$r_m" "$num_monomers" | tee /dev/stderr)

# Extract the total number of monomers from the output
total_particles=$(echo "${aggregate_output}" | grep "TOTAL_PARTICLES:" | cut -d':' -f2)
echo $total_particles
# Define the input directory and aggregate name
INPUT_DIR=~/runs/Agg_models/input
AGG_NAME="aggregate_NIR${total_particles}_${r_m}.k"
AGG_PATH="$INPUT_DIR/$AGG_NAME"
NEW_RUN_DIR=~/runs/Agg_models/NEW_RUNS

# Create a new directory for the current model
NEW_DIR="AGG_${total_particles}_${r_m}"
MODEL_DIR="$NEW_RUN_DIR/$NEW_DIR"

#if [ -d "$MODEL_DIR" ]; then
#  read -p "WARNING: Overwriting existing directory? (Y/N) " yn
#  case $yn in
#    [Yy]* ) ;;
#    [Nn]* ) exit;;
#    * ) echo "Please answer yes or no.";;
#  esac
#fi
#this is for the grid to automately skip 
if [ -d "$MODEL_DIR" ]; then
  echo "WARNING: Directory $MODEL_DIR already exists. Skipping this run."
  exit
fi

mkdir -p "$MODEL_DIR"

cp "$AGG_PATH" "$MODEL_DIR"
cp "$GMM_PATH/gmm01f.in" "$MODEL_DIR"
cp "$GMM_PATH/gmm01s.f90" "$MODEL_DIR"
cp "$GMM_PATH/gmm01f.par" "$MODEL_DIR"

# Create run.q file
RUNQ_FILE="$MODEL_DIR/run.q"
cat > "$RUNQ_FILE" << EOL
#!/bin/bash
#SBATCH --nodes=1 --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --job-name=GMM_${total_particles}_${r_m}
#SBATCH --output=AGG_${total_particles}_${r_m}.in.qout
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --time=0-168:00:00
#SBATCH --mail-type=END
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

./gmm_${total_particles}.out
EOL

chmod +x "$RUNQ_FILE"

sed -i "1s|.*|$(basename $AGG_PATH)|" "$MODEL_DIR/gmm01f.in"
#sed -i "1s|.*/|$AGG_PATH|" "$INPUT_DIR/$NEW_DIR/gmm01f.in"

ifort -o "gmm_${total_particles}.out" "$MODEL_DIR/gmm01s.f90"

mv "gmm_${total_particles}.out" "$MODEL_DIR"

# Submit the run.q file
cd "$MODEL_DIR"
sbatch run.q
