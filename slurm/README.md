# Batch processing of HEXIMAP on CHPC

There are a number of specific requirements to get set up to run the included slurm scripts on CHPC resources.
This document details those steps so that someone who is reasonably familiar with CHPC resources and processing should be able to reproduce the methods.

A number of additional modules are required for processing hexagon imagery.
Some of these resources are available across all Utah CHPC environments (e.g. MATLAB, customized mexopencv library, etc.) while others will need to be installed on each individual user's home partition that wishes to run the scripts (e.g. Python conda environment).

The first step is therefore to get these modules set up.
This code assumes a specific conda environment called `heximap` is in the reachable conda paths.
To get conda set up on your individual home partition, follow the instructions for installing and using miniconda found on the [CHPC website](https://www.chpc.utah.edu/documentation/software/python-anaconda.php#mi).

Once you have `conda` set up, you will need to create the `heximap` environment.
Do this by navigating to the heximap root directory (the root of the repository cloned from GitHub) and run the following command in the commandline terminal:
```{bash}
conda env create -f environment.yml
```
This will install all the necessary Python modules/versions to run `heximap` and will install them into a conda environment of the same name.

You will then need to edit two scripts for your particular use case, `run-heximap.slrm` and `chpc-main.py`.
Modify the `#SBATCH` lines of `run-heximap.slrm` to match your needs for your specific processing task (e.g. account name, nodes to use, walltime, etc.).
You will also need to edit the input paths for the source, data, and output directories.
The source directory is the path to the `heximap` repository, the data directory is the path to the hexagon data (both raw imagery and metadata),
<!-- NOTE: Need to explain the proper structure for the data directory -->
and the output directory is the path to save the final processed results (this folder does not need to exist, as it will be created if necessary during processing).
