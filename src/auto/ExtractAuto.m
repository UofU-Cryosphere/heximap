% Function to automatically extract DEM values from a pair of Hexagon
% images.

function [dem] = ExtractAuto(strHexPath1, strPtsPath1, strHexPath2, strPtsPath2)

% Define Hexagon pair mat files
ObjM1 = matfile(strHexPath1, 'Writable',true);
ObjM2 = matfile(strHexPath2, 'Writable',true);

% Define rough corner positions for initial georeferencing and alignment
EXT_FUNC.ControlPoints(ObjM1, strPtsPath1);
EXT_FUNC.ControlPoints(ObjM2, strPtsPath2);

end