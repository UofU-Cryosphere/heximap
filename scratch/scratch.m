% Experimental workspace

% Load environment settings
addpath(genpath('/uufs/chcpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2018a-nomkl'))
% addpath(genpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2018a-nomkl/opencv_contrib/'))
addpath(genpath('~/Codebase/heximap'))

% Define path to image directory and (stitched) image name
im_dir = '~/Documents/Research/GSLR/tutorial/npi/heximap/DZB1212-500105L005001_6001/';
im_name = 'DZB1212-500105L006001';

% Define image halves
sInfoL = imfinfo(char(fullfile(im_dir, strcat(im_name, "_a.tif"))));
sInfoR = imfinfo(char(fullfile(im_dir, strcat(im_name, "_b.tif"))));

test = StitchAuto(sInfoL, sInfoR);


% corners = STI_FUNC.GetCorners(im_dir, im_name);