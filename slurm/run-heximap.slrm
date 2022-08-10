#!/bin/bash
#SBATCH --time=06:00:00 # walltime, abbreviated by -t
#SBATCH --nodes=1      # number of cluster nodes, abbreviated by -N
#SBATCH -o slurm-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH --ntasks=1    # number of MPI tasks, abbreviated by -n
#SBATCH --cpus-per-task=8
# #SBATCH --constraint="c20"
#SBATCH --mem=256G
# additional information for allocated clusters
#SBATCH --account=rupper     # account - abbreviated by -A
#SBATCH --partition=lonepeak  # partition, abbreviated by -p
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=durban.keeler@gmail.com

###### USER-DEFINED INPUTS ######
# Define paths to source, data, and output directories
StartDir=$(pwd)
HexDir="/uufs/chpc.utah.edu/common/home/u1046484/Codebase/heximap/"
SourceDir=$HexDir"src/"
ScratchDir="/scratch/general/vast/u1046484/heximap/"
DataDir="/uufs/chpc.utah.edu/common/home/u1046484/Documents/Research/GSLR/data/"
ImageDir=$DataDir"hexagon/declass-ii/imagery/"
GeoDir=$DataDir"refDEMs/"
OutDir=$StartDir/"outputs/"
PyScriptPath=$HexDir"scripts/chpc-main.py"
MatScriptPath=$HexDir"scripts/HPC_HEX.m"
###### USER-DEFINED INPUTS ######

echo "LOADING MODULES..."
# Purge old modules, load required modules, and activate conda env
module purge
module load gcc/8.5.0 matlab/R2022a opencv/3.4.1-mex-nomkl mexopencv/R2022a-nomkl
module use $HOME/MyModules
module load miniconda3/latest
conda activate heximap

# Set up scratch directory with necessary contents
mkdir -p $ScratchDir
cp $PyScriptPath $ScratchDir
cp $MatScriptPath $ScratchDir
cd $ScratchDir

echo "RUNNING PYTHON SCRIPT"
# Run Python main script
python chpc-main.py $SourceDir $DataDir $ImageDir $GeoDir $OutDir

echo "RUNNING MATLAB SCRIPT"
# Run MATLAB main script
matlab -nodisplay -batch 'HPC_HEX' -logfile matlab_stdout.log
