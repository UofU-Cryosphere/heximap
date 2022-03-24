% Function to get image corners from HEX image when given directory and
% filename. This will automatically select left and right image halves and
% extract out corners for the stitched image.

function corners = GetCorners(im_dir, im_name)

% Define corner region subset size
nrow = 5000;
ncol = 5000;

% Define image halves
sInfoL = imfinfo(char(fullfile(im_dir, strcat(im_name, "_a.tif"))));
sInfoR = imfinfo(char(fullfile(im_dir, strcat(im_name, "_b.tif"))));


NW_region = imread(sInfoL.Filename,...
    'PixelRegion',{[1 nrow] [1 ncol]});
region_corner = findCorner(NW_region);
NW_corner = region_corner;

SW_region = imread(sInfoL.Filename,...
    'PixelRegion',{[sInfoL.Height-nrow sInfoL.Height] [1 ncol]});
region_corner = findCorner(SW_region);
SW_corner = [(sInfoL.Height-nrow+region_corner(1)) region_corner(2)];

NE_region = imread(sInfoR.Filename,...
    'PixelRegion',{[1 nrow] [sInfoR.Width-ncol sInfoR.Width]});
region_corner = findCorner(NE_region);
NE_corner = [region_corner(1) (sInfoR.Width-ncol+region_corner(2))];

SE_region = imread(sInfoR.Filename,...
    'PixelRegion', {[sInfoR.Height-nrow sInfoR.Height] [sInfoR.Width-ncol sInfoR.Width]});
region_corner = findCorner(SE_region);
SE_corner = [(sInfoR.Height-nrow+region_corner(1)) (sInfoR.Width-ncol+region_corner(2))];


% W_idx = round(mean([NW_corner(2) SW_corner(2)]));
N_diff = abs(NW_corner(1) - NE_corner(1));
S_diff = abs(SW_corner(1) - SE_corner(1));
E_diff = abs(NE_corner(2) - SE_corner(2));
W_diff = abs(NW_corner(2) - SW_corner(2));


corners = [];


function corner_loc = findCorner(IM_subset)

% Define defaults
gauss_sigma = 50;
text_threshold = 250;


% Set pixels above threshold (should represent text) to zero
IM_subset(IM_subset > text_threshold) = 0;

% Convert to logical and apply Gaussian convolution
dat_log = double(logical(IM_subset));
dat_blur = imgaussfilt(dat_log, gauss_sigma);

% Find max gradient in columns
[~, c_idx] = max(abs(diff(sum(dat_blur,1))));

% Find max gradient in rows
[~, r_idx] = max(abs(diff(sum(dat_blur,2))));

% % Diagnostic plot
% [nrow, ncol] = size(IM_subset);
% figure
% imshow(IM_subset)
% hold on
% plot(1:ncol, repmat(r_idx,ncol,1), 'r', 'LineWidth', 1)
% plot(repmat(c_idx, nrow,1), 1:nrow,'r', 'LineWidth', 1)
% hold off
% pause(1) % Issues with plotting require pausing

corner_loc = [r_idx, c_idx];

end
end