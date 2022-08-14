function [mArray] = InpaintNaN_chunks(mArray, iMaxSize, iInpaintMethod)
%INPAINTNAN_CHUNK Helper function for inpaint_nan that will process the
%data in chunks rather than the full array if the array is larger than a
%given threshold. This helps avoid out-of-memory issues when dealing with
%large matrices with a high number of NaN values.

%   Detailed explanation goes here

% If mArray is overly large, process nan-filling in chunks
if any(size(mArray) > iMaxSize)

    % Determine chunking size for each dimension
    [iY, iX] = size(mArray);
    iYstep = min([iY iMaxSize]);
    iXstep = min([iX iMaxSize]);

    % Get indices for chunk edges in rows
    vY = 1:iYstep:iY;
    if iY ~= vY(end)
        vY = [vY iY];
    end

    % Get indices for chunk edges in columns
    vX = 1:iXstep:iX;
    if iX ~= vX(end)
        vX = [vX iX];
    end
    
    % Inpaint nan values for each chunk
    for i=1:length(vY)-1
        for j=1:length(vX)-1
            mArray(vY(i):vY(i+1),vX(j):vX(j+1)) = inpaint_nans(...
                mArray(vY(i):vY(i+1),vX(j):vX(j+1)), iInpaintMethod);
        end
    end

else
    % When array is sufficiently small, inpaint nan values on full array
    mArray = inpaint_nans(mArray);
end
end

