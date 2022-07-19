% Function to perform georeferencing of Hexagon DEM using a reference DEM
function [] = GeorefAuto(strHexPath, strDEMrefPath, strShpPath)

% Check input data for errors
geoCheckInput(strDEMrefPath,strShpPath);

% Get mat files containing data for image windows
[cL,cR] = GEO_FUNC.GetMatFiles(strHexPath);

