% Experimental workspace

% Load environment settings
addpath(genpath('~/Codebase/heximap'))
% addpath(genpath('/uufs/chcpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2018a-nomkl'))
% addpath(genpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2018a-nomkl/opencv_contrib/'))
addpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2022a-nomkl/')
addpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2022a-nomkl/opencv_contrib/')

% Define path directory containing hexagon imagery and define names of
% (stiched) images in hex pair
IM_dir = '~/Documents/Research/GSLR/data/hexagon/declass-ii/imagery/';
IM1_name = 'DZB1212-500105L005001';
IM2_name = 'DZB1212-500105L006001';

% Define path of directory to save outputs
ex_dir = '~/Documents/Research/GSLR/data/interim/';

% Define image halves for first image in pair
sInfoL1 = imfinfo(char(fullfile(IM_dir, strcat(IM1_name, "_a.tif"))));
sInfoR1 = imfinfo(char(fullfile(IM_dir, strcat(IM1_name, "_b.tif"))));

% Stitch together the two image halves of first image and save to ex_dir
objIM1 = StitchAuto(sInfoL1, sInfoR1, ex_dir);

% Define image halves for second image in pair
sInfoL2 = imfinfo(char(fullfile(IM_dir, strcat(IM2_name, "_a.tif"))));
sInfoR2 = imfinfo(char(fullfile(IM_dir, strcat(IM2_name, "_b.tif"))));

% Stitch together the two image halves of second image and save to ex_dir
objIM2 = StitchAuto(sInfoL2, sInfoR2, ex_dir);

% Import image metadata (generated with Python script)
IM1_meta = load(fullfile("data/tmp/", strcat(IM1_name, "_meta.mat")));
IM2_meta = load(fullfile("data/tmp/", strcat(IM2_name, "_meta.mat")));

% Extract DEM for overlap between hexagon pair
objDEM = ExtractAuto(objIM1, IM1_meta, objIM2, IM2_meta);