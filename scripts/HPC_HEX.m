%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

now0 = tic();

% Load script parameters from Python-generated file
sParams = load("tmp/sParams.mat");

% Load environment settings
addpath(genpath(sParams.strSourcePath))
addpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2022a-nomkl/')
addpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2022a-nomkl/opencv_contrib/')

%%%%%% USER-DEFINED INPUTS %%%%%%
% Define path directory containing hexagon imagery and define names of
% (stiched) images in hex pair (this will later be provided programatically
% from Python call)
% IM_dir = '~/Documents/Research/GSLR/data/hexagon/declass-ii/imagery/';
% IM1_name = 'DZB1212-500105L005001';
% IM2_name = 'DZB1212-500105L006001';

% % Define heximap defaults (will later re-implement so that these can be
% % defined outside the function and the defaults overridden)
% strRes = '1/2'; %Desired image resolution for stereo matching
% iBlkSz = 3; %The block size (in pixels) used for matching in disparity mapping (should be odd 3-11)
% verbose = false; %Should iterim data be kept after processing?

%%%%%% USER-DEFINED INPUTS %%%%%%

% Define tmp processing directory (create if necessary)
strTmpPath = fullfile(sParams.strRootPath, 'tmp/');
if ~isfolder(strTmpPath)
    mkdir(strTmpPath)
end

% Define image halves for first image in pair
sInfoL1 = imfinfo(char(fullfile(sParams.strImageDir, ...
    strcat(sParams.strIM1Name, "_a.tif"))));
sInfoR1 = imfinfo(char(fullfile(sParams.strImageDir, ...
    strcat(sParams.strIM1Name, "_b.tif"))));

% Define image halves for second image in pair
sInfoL2 = imfinfo(char(fullfile(sParams.strImageDir, ...
    strcat(sParams.strIM2Name, "_a.tif"))));
sInfoR2 = imfinfo(char(fullfile(sParams.strImageDir, ...
    strcat(sParams.strIM2Name, "_b.tif"))));

disp('Stitching Hexagon images...\n')
now1 = tic();

% Stitch together the two image halves of each image and save to scratch
strStitchPath = fullfile(strTmpPath, "stitched/");
if ~isfolder(strStitchPath)
    mkdir(strStitchPath)
end

fprintf('Stitching halves for Image %s...', sParams.strIM1Name)
objIM1 = StitchAuto(sInfoL1, sInfoR1, strStitchPath);
fprintf('Stitching halves for Image %s...', sParams.strIM2Name)
objIM2 = StitchAuto(sInfoL2, sInfoR2, strStitchPath);

% % For development purposes, skip stitching and load previously stitched
% % images
% objIM1 = matfile(fullfile(strStitchPath, strcat(sParams.strIM1Name,".mat")),...
%     'Writable',true);
% objIM2 = matfile(fullfile(strStitchPath, strcat(sParams.strIM2Name,".mat")),...
%     'Writable',true);

T_stitch = toc(now1);
fprintf('Total stitching time: %.0f seconds', T_stitch)

% Import image metadata (generated with Python script)
IM1_meta = load(fullfile(strTmpPath, "metadata", ...
    strcat(sParams.strIM1Name, "_meta.mat")));
IM2_meta = load(fullfile(strTmpPath, "metadata", ...
    strcat(sParams.strIM2Name, "_meta.mat")));

% Load ROI locations from Python-generated file
ROIs = load(fullfile(strTmpPath, "hexROIs.mat"));

disp('Begin extraction of DEM from Hexagon imagery pair...\n')
now2 = tic();
% Extract DEM for overlap between hexagon pair
ExtractAuto(objIM1, IM1_meta, objIM2, IM2_meta, ROIs, strTmpPath);
% ExtractAuto(objIM1, IM1_meta, objIM2, IM2_meta, strTmpPath,...
%     strRes, iBlkSz);

T_extract = toc(now2);
fprintf('Extraction time: %.0f seconds', T_extract)

disp('Georeferencing hexagon DEM...')
% now3 = tic()
% Georeference DEMs using reference DEM
tmp = sParams.strGeoRefPath;

% T_ref = toc(now3);
% fprintf('Georeferencing time: %.0f seconds', T_ref)

T_total = toc(now0);
fprintf('Total time: %.0f seconds', T_total)
