% Experimental workspace

% Load environment settings
addpath(genpath('~/Codebase/heximap'))
% addpath(genpath('/uufs/chcpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2018a-nomkl'))
% addpath(genpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2018a-nomkl/opencv_contrib/'))
addpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2022a-nomkl/')
addpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2022a-nomkl/opencv_contrib/')

%%%%%% USER-DEFINED INPUTS %%%%%%
% Define path directory containing hexagon imagery and define names of
% (stiched) images in hex pair (this will later be provided programatically
% from Python call)
IM_dir = '~/Documents/Research/GSLR/data/hexagon/declass-ii/imagery/';
IM1_name = 'DZB1212-500105L005001';
IM2_name = 'DZB1212-500105L006001';

% Define path for DEM to be used for georeferencing
GeorefPath = '';

% Define path of directory to save interim results (create dir if needed)
ScratchPath = '~/Codebase/heximap/scratch/';
if ~isfolder(ScratchPath)
    mkdir(ScratchPath)
end

% Define path of directory to save final results (create dir if needed)
ExportPath = '~/Documents/Research/GSLR/data/';
if ~isfolder(ExportPath)
    mkdir(ExportPath)
end

% % Define heximap defaults (will later re-implement so that these can be
% % defined outside the function and the defaults overridden)
% strRes = '1/2'; %Desired image resolution for stereo matching
% iBlkSz = 3; %The block size (in pixels) used for matching in disparity mapping (should be odd 3-11)
% verbose = false; %Should iterim data be kept after processing?

%%%%%% USER-DEFINED INPUTS %%%%%%


% Define image halves for first image in pair
sInfoL1 = imfinfo(char(fullfile(IM_dir, strcat(IM1_name, "_a.tif"))));
sInfoR1 = imfinfo(char(fullfile(IM_dir, strcat(IM1_name, "_b.tif"))));

% Define image halves for second image in pair
sInfoL2 = imfinfo(char(fullfile(IM_dir, strcat(IM2_name, "_a.tif"))));
sInfoR2 = imfinfo(char(fullfile(IM_dir, strcat(IM2_name, "_b.tif"))));

% Stitch together the two image halves of each image and save to scratch
StitchPath = fullfile(ScratchPath, "stitched/");
if ~isfolder(StitchPath)
    mkdir(StitchPath)
end
% objIM1 = StitchAuto(sInfoL1, sInfoR1, StitchPath);
% objIM2 = StitchAuto(sInfoL2, sInfoR2, StitchPath);
% For development purposes, skip stitching and load previously stitched
% images
objIM1 = matfile(fullfile(StitchPath, strcat(IM1_name,".mat")),...
    'Writable',true);
objIM2 = matfile(fullfile(StitchPath, strcat(IM2_name,".mat")),...
    'Writable',true);

% Import image metadata (generated with Python script)
IM1_meta = load(fullfile(ScratchPath, "metadata", ...
    strcat(IM1_name, "_meta.mat")));
IM2_meta = load(fullfile(ScratchPath, "metadata", ...
    strcat(IM2_name, "_meta.mat")));
% IM1_meta = load(fullfile("data/tmp/", strcat(IM1_name, "_meta.mat")));
% IM2_meta = load(fullfile("data/tmp/", strcat(IM2_name, "_meta.mat")));

% Load ROI locations from Python-generated file
ROIs = load(fullfile(ScratchPath, "hexROIs.mat"));

% Extract DEM for overlap between hexagon pair
ExtractAuto(objIM1, IM1_meta, objIM2, IM2_meta, ROIs, ScratchPath);
% ExtractAuto(objIM1, IM1_meta, objIM2, IM2_meta, ScratchPath,...
%     strRes, iBlkSz);