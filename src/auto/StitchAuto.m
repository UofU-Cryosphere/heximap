% Function to perform image stitching in a fully automated manner
function [objM] = StitchAuto(sInfoL, sInfoR, export_path)

% Find corners of image halves
corners = STI_FUNC.GetCorners(sInfoL, sInfoR);

% Stitch together the image halves
objM = STI_FUNC.Stitch(sInfoL, sInfoR, corners, export_path);

% % Save lower resolution copies of the images
% stiResize(objM, []);

% Save camera scan resolution, focal length, and estimated principal
% point in mat file
iScanRes = 7;
objM.ScanResMicrometers = iScanRes;
objM.FocalLengthMicrometers = 304.8E3;
objM.FocalLengthPixels = 304.8E3 / iScanRes;
objM.PrincipalPointPixels = fliplr(size(objM,'Image')/2);

end