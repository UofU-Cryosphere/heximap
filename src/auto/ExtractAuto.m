% Function to automatically extract DEM values from a pair of Hexagon
% images.

function [objDEM] = ExtractAuto(objM1, IM1_meta, objM2, IM2_meta)

% Define rough corner positions for initial georeferencing and alignment
EXT_FUNC.ControlPoints(objM1, IM1_meta);
EXT_FUNC.ControlPoints(objM2, IM2_meta);

% Sort the Hexagon images to establish left and right images when computing
% stereo disparity maps
[vOrder] = EXT_FUNC.SortImages({objM1, objM2});

% Compute fundamental matrix, relative pose matrices, and rough
% homographies using full images. More accurate matrices are computed
% later for each processing window.
EXT_FUNC.InitTrans(cM,strHexPath);

end