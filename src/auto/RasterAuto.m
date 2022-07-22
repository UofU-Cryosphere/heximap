function [] = RasterAuto(strHexPath, strSavePath, sRasterParams)
%RASTERAUTO Summary of this function goes here
%   Detailed explanation goes here


% Get mat files containing data for image windows
cL = rasGetMatFiles(strHexPath);

% Create directories for completed dems and images
mkdir(strSavePath,'dems')
mkdir(strSavePath,'images')

% Loop through each window
for iW = 1:numel(cL)

    try

        % Rasterize the DEM
        rasDem(cL{iW}, strHexPath, sRasterParams.lClean, sRasterParams.lMed,...
            sRasterParams.lDen, sRasterParams.iGap, sRasterParams.iSpec, ...
            sRasterParams.iMed, sRasterParams.iDenT, sRasterParams.iDenN);

        % Rasterize the orthoimage
        rasOrtho(cL{iW}, strHexPath);

    catch objExc

        warning(objExc.message)
        warning(['An error occurred. Skipping window ' num2str(iW) '...'])
        cL{iW}.Error = objExc.message;

    end
end

end

