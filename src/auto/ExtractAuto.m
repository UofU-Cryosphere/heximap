% Function to automatically extract DEM values from a pair of Hexagon
% images.

function [] = ExtractAuto(objM1, IM1_meta, objM2, IM2_meta, ROIs, ...
    strExportPath, strRes, iBlkSz)

arguments
    objM1
    IM1_meta
    objM2
    IM2_meta
    ROIs
    strExportPath {mustBeTextScalar}
    strRes {mustBeTextScalar} = '1/2';
    iBlkSz {mustBeInteger} = 3;
end

% Define export directory (create if needed)
strExPath = fullfile(strExportPath, "extraction/");
if ~isfolder(strExPath)
    mkdir(strExPath)
end

% Define rough corner positions for initial georeferencing and alignment
EXT_FUNC.ControlPoints(objM1, IM1_meta);
EXT_FUNC.ControlPoints(objM2, IM2_meta);

% Sort the Hexagon images to establish left and right images when computing
% stereo disparity maps
[vOrder] = EXT_FUNC.SortImages({objM1, objM2});

% Define left-right ordering for image pair
cM = {objM1, objM2};
cM = cM(vOrder);

% Compute fundamental matrix, relative pose matrices, and rough
% homographies using full images. More accurate matrices are computed
% later for each processing window.
EXT_FUNC.InitTrans(cM);

% Define regions of interest within overlap between pairs
% This will eventually subset to Python-generated ROIs based on glacier
% locations
cWindows = EXT_FUNC.GetRegions(cM, ROIs);

% % Until that is implemented, here is a manually-curated variable for
% % development and testing (a 2x2 window region near the center)
% sWin1 = struct(...
%     'left', [32251, 9388; 37146, 14350], ...
%     'right', [13243, 9288; 18138, 14250], ...
%     'region', 1);
% sWin2 = struct(...
%     'left', [32251, 13787; 37146, 18751], ...
%     'right', [13228, 13653; 18123, 18617], ...
%     'region', 1);
% sWin3 = struct(...
%     'left', [36646, 9386; 41549, 14353], ...
%     'right', [17565, 9286; 22468, 14253], ...
%     'region', 1);
% sWin4 = struct(...
%     'left', [36647, 13784; 41549, 18753], ...
%     'right', [17551, 13655; 22453, 18624], ...
%     'region', 1);
% cWindow = {sWin1, sWin2, sWin3, sWin4};

now0 = tic();
% Rectify the stereo images, compute disparity maps
EXT_FUNC.DisparityLoop(cM,strExPath,cWindows(1:4),strRes,iBlkSz);
% EXT_FUNC.DisparityLoop(cM,strExPath,cWindows,strRes,iBlkSz);
T_disp = toc(now0);
fprintf('Disparity loop time: %.0f seconds', T_disp)

now1 = tic();
% Refine camera orientations using bundle adjustment
EXT_FUNC.BundleAdjustLoop(strExPath);

T_adjust = toc(now0);
fprintf('Bundle adjustment time: %.0f seconds', T_adjust)

now2 = tic();
% Triangulate the points
EXT_FUNC.TriangulateLoop(strExPath);

T_tri = toc(now2);
fprintf('Bundle adjustment time: %.0f seconds', T_tri)

end