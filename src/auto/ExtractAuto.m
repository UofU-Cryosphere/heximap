% Function to automatically extract DEM values from a pair of Hexagon
% images.

function [objDEM] = ExtractAuto(objM1, IM1_meta, objM2, IM2_meta)

% Define rough corner positions for initial georeferencing and alignment
EXT_FUNC.ControlPoints(objM1, IM1_meta);
EXT_FUNC.ControlPoints(objM2, IM2_meta);

% Sort the Hexagon images to establish left and right images when computing
% stereo disparity maps
[cHexFile,cM] = EXT_FUNC.SortImages(cHexFile,cM,hW);

end