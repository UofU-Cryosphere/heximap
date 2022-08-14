% Function to perform image stitching in a fully automated manner
function [objM] = StitchAuto(sInfoL, sInfoR, export_path)

% Find corners of image halves
corners = STI_FUNC.GetCorners(sInfoL, sInfoR);

% Stitch together the image halves
objM = STI_FUNC.Stitch(sInfoL, sInfoR, corners, export_path);

fprintf('Upscaling stitched Hexagon Image...\n')
now = tic();
% Save lower resolution copies of the images
STI_FUNC.Resize(objM);

T_ref = toc(now);
fprintf('Upscaling time: %.0f seconds\n', T_ref)

% Save camera scan resolution, focal length, and estimated principal
% point in mat file
iScanRes = 7;
objM.ScanResMicrometers = iScanRes;
objM.FocalLengthMicrometers = 304.8E3;
objM.FocalLengthPixels = 304.8E3 / iScanRes;
objM.PrincipalPointPixels = fliplr(size(objM,'Image')/2);

end