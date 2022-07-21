% Function to perform georeferencing of Hexagon DEM using a reference DEM
function [] = GeorefAuto(strHexPath, strDEMrefPath, strShpPath)

% Check input data for errors
checkInput(strDEMrefPath)
checkInput(strShpPath)

% Get mat files containing data for image windows
[cL,cR] = geoGetMatFiles(strHexPath);

% Select window for initial georeferencing
vRankIdx = GEO_FUNC.RankWindows(cL,cR);
