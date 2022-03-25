% Experimental workspace

% Load environment settings
addpath(genpath('/uufs/chcpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2018a-nomkl'))
% addpath(genpath('/uufs/chpc.utah.edu/sys/srcdir/mexopencv/3.4.1-R2018a-nomkl/opencv_contrib/'))
addpath(genpath('~/Codebase/heximap'))

% Define path to image directory and (stitched) image name
im_dir = '../tutorial/npi/heximap/DZB1212-500105L005001_6001/';
im_name = 'DZB1212-500105L006001';

corners = STI_FUNC.GetCorners(im_dir, im_name);