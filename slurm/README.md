# Batch processing of HEXIMAP on CHPC

There are a number of specific requirements to get set up to run the included batch processing script on CHPC resources.
This document details those steps so that someone who is reasonably familiar with CHPC resources and processing should be able to reproduce the methods.

A number of additional modules are required for processing hexagon imagery.
Some of these resources are available across all Utah CHPC environments (e.g. MATLAB, customized mexopencv library, etc.) while others will need to be installed on each individual user's home partition that wishes to run the scripts (e.g. the appropriate Python conda environment).

The first step is therefore to get these modules set up.
This code assumes a specific conda environment called "heximap" is at a reachable conda path.
To get conda set up on your individual home partition, follow the instructions for installing and using miniconda found on the [CHPC website](https://www.chpc.utah.edu/documentation/software/python-anaconda.php#mi).

Once you have conda set up, you will need to create the "heximap" environment.
Do this by navigating to the `heximap` root directory (the root of the repository cloned from GitHub) and run the following command in the commandline terminal:
```{bash}
conda env create -f environment.yml
```
This command creates a new virtual environment based on the necessary packages listed in the `environment.yml` file (this file is included as part of the `heximap` repository).
It therefore installs all the necessary Python modules/versions to run `heximap` and places them into a conda environment of the same name.

<!-- To complete the conda environment setup, you will also need to install the MATLAB Python engine to your environment.
Instructions for how to do this are found in [this answer](https://www.mathworks.com/matlabcentral/answers/346068-how-do-i-properly-install-matlab-engine-using-the-anaconda-package-manager-for-python).
Just make sure that you run the commands from within the `heximap` environment (e.g. run `conda activate heximap` prior to running the steps in that answer). -->

You will then need to edit two scripts for your particular use case, `run-heximap.slrm` (found in the `slurm` directory of the `heximap` repository) and `chpc-main.py` (found in the `scripts` subdirectory).
Modify the `#SBATCH` lines of `run-heximap.slrm` to match your needs for your specific processing task (e.g. account name, nodes to use, walltime, etc.).
You will also need to edit the input paths for `HexDir` (where you have the `heximap` repository stored on your home partition), `ScratchDir` (where should `heximap` create the temporary processing directory---this will usually be in `/scratch/general/vast/` under a subdirectory of your CHPC username e.g. u1234567), `DataDir`[^datadir] (the path to where you have hexagon imagery and glacier data saved), and `OutDir` (the location you wish to save the final processed results).
The `ScratchDir` and the `OutDir` do not need to exist prior to running the script.
They both will be created as necessary during `heximap` processing.

[^datadir]: This directory should follow a particular format in order for everything to run properly. The directory layout should exactly match the following form:
    ```{nohighlight}
    {/Path/to/Data/Dir/}
      ├──hexagon
      │   ├──declass-ii
      │   |   ├──imagery
      │   |   |   ├──hex-ii-file-downloaded-from-EarthExplorer_a.tif
      │   |   |   ├──hex-ii-file-downloaded-from-EarthExplorer_a.tif
      │   |   |   ├── ...
      │   |   |   ├──more-files-from-EarthExplorer_a.tif
      │   |   |   └──more-files-from-EarthExplorer_b.tif
      │   |   └── metadata
      │   |       └──declass-ii-shapefiles-from-EarthExplorer.shp
      │   └──declass-iii
      │       └── ...{same structure as declass-ii}
      └── RGI-data
          ├──RGI-region-folder
          |   └──RGI-region-shapefile.shp
          ├── ...
          └──More-RGI-regions-folders
              └──More-RGI-regions-shapefile.shp
    ```
    It may also be necessary to modify some of the Path assignments in `chpc-main.py` to match the specific folder names on your system.

Some of the details of what and how you are specifically processing Hexagon imagery will also impact the precise makeup of the `chpc-main.py` script.
You will therefore likely need to make some minor modifications to this as well.
As exact specifications will be case-dependent, we cannot state precisely what all will need to change and be updated in this file.
We do provide some general guidance as to which aspects will typically change from task to task.
Most of the required edits relate to defining what specific regions, glaciers, and attributes are of interest in your specific task.
