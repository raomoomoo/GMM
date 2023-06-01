#!/bin/bash

module load numpy/1.18.3-python-3.8.5

GMM_PATH=~/GMM
export GMM_PATH

read -p "Enter the number of monomers: " n_m
read -p "Enter the number of layers: " r_m
r_m=${r_m:-10}

python3 ~/GMM/Multilayer_Agg_Model.py

INPUT_DIR=~/runs/Agg_models/input
AGG_NAME="aggregate_NIR_${n_m}_${r_m}.k"
AGG_PATH="$INPUT_DIR/$AGG_NAME"

NEW_DIR="AGG_${n_m}_${r_m}"

if [ -d "$INPUT_DIR/$NEW_DIR" ]; then
  read -p "WARNING: Overwriting existing directory? (Y/N) " yn
  case $yn in
    [Yy]* ) ;;
    [Nn]* ) exit;;
    * ) echo "Please answer yes or no.";;
  esac
fi

mkdir -p "$INPUT_DIR/$NEW_DIR"

cp "$AGG_PATH" "$INPUT_DIR/$NEW_DIR"
cp "$GMM_PATH/gmm01f.in" "$INPUT_DIR/$NEW_DIR"
cp "$GMM_PATH/gmm01s.f90" "$INPUT_DIR/$NEW_DIR"
cp "$GMM_PATH/gmm01f.par" "$INPUT_DIR/$NEW_DIR"

# Create run.q file
RUNQ_FILE="$INPUT_DIR/$NEW_DIR/run.q"
cat > "$RUNQ_FILE" << EOL
#!/bin/bash
#SBATCH --nodes=1 --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --job-name=GMM_${n_m}_${r_m}
#SBATCH --output=AGG_${n_m}_${r_m}.in.qout
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --time=0-168:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=rmor0013@student.monash.edu
#SBATCH --mem=16G
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

./gmm_${n_m}.out
EOL

chmod +x "$RUNQ_FILE"

sed -i "1s|.*/|$AGG_PATH|" "$INPUT_DIR/$NEW_DIR/gmm01f.in"

ifort -o "gmm_${n_m}.out" "$INPUT_DIR/$NEW_DIR/gmm01s.f90"

mv "gmm_${n_m}.out" "$INPUT_DIR/$NEW_DIR"

# Submit the run.q file
cd "$INPUT_DIR/$NEW_DIR"
sbatch run.q
