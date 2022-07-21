% Function to perform georeferencing of Hexagon DEM using a reference DEM
function [] = GeorefAuto(strHexPath, strDEMrefPath, strShpPath)

% Turn off lVis (this is only used in manual processing to visualize
% georeferencing
lVis = false;

% Check input data for errors
checkInput(strDEMrefPath)
checkShpPath(strShpPath);

% Get mat files containing data for image windows
[cL,cR] = geoGetMatFiles(strHexPath);

% Select window for initial georeferencing
vRankIdx = GEO_FUNC.RankWindows(cL,cR);

% Define control points
geoControlPoints(cL{vRankIdx(1)}, strDEMrefPath, 'full-auto');

% Compute transformations for initial georeferencing (using chosen window)
geoInitTrans(cL{vRankIdx(1)},strDEMrefPath,strShpPath,lVis);

% Apply georeferencing transformations (from chosen window to all windows)
geoApplyTrans(cL,cR,vRankIdx(1));

% Optimize orientation of all regions
geoOptimize(cL,vRankIdx(1),strRef,strShpPath,lVis,lPoly,lWin,dMaxDelZ,hW);

end
