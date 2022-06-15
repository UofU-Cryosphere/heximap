% Experimental workspace

% Load environment settings
addpath(genpath('~/Codebase/heximap'))
% addpath(genpath('/uufs/chcpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2018a-nomkl'))
% addpath(genpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2018a-nomkl/opencv_contrib/'))
addpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2022a-nomkl/')
addpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2022a-nomkl/opencv_contrib/')

% Define path to image directory and (stitched) image name
im_dir = ['~/Documents/Research/GSLR/tutorial/npi/'...
    'heximap/DZB1212-500105L005001_6001/'];
im_names = {'DZB1212-500105L005001', 'DZB1212-500105L006001'};

% Define path of directory to save outputs
ex_dir = '~/Documents/Research/GSLR/data/interim/';

for i=1:length(im_names)
    
    % Define image halves
    sInfoL = imfinfo(char(fullfile(...
        im_dir, strcat(im_names{i}, "_a.tif"))));
    sInfoR = imfinfo(char(fullfile(...
        im_dir, strcat(im_names{i}, "_b.tif"))));
    
    % Stitch together the two image halves and save to ex_dir
    StitchAuto(sInfoL, sInfoR, ex_dir);
    
end
