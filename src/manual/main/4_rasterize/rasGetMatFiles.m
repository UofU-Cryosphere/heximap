function cL = rasGetMatFiles(strWinPath,varargin)
% Get mat files containing image windows

% Retain interactive plotting when in manual mode
if nargin > 1
    hW = varargin{1};
    try
        % Update waitbar
        waitbar(0/1,hW,'checking .mat files...')
    catch
    end
end

% Get paths to mat files
cL = cellfun(@(x) matfile(x,'Writable',true),getFiles(strWinPath, ...
    'Left.mat'),'Uni',0);

% Remove any windows where georeferencing was not completed
lIn = true(size(cL));
for i = 1:numel(cL)
    sAccuracy = cL{i}.Accuracy;
    if ~isfield(sAccuracy,'VerticalRMSError')
        lIn(i) = false;
        warning(['Georeferencing was not completed for window ' ...
            num2str(cL{i}.WindowID) '. Skipping.'])
    end
end
cL = cL(lIn);

% Make sure at least one valid window exists
if isempty(cL)
    error('No valid windows found in the specified folder.')
end

% Sort the windows
[~,vIdx] = sort(cell2mat(cellfun(@(x) x.WindowID,cL,'Uni',0)));
cL = cL(vIdx);
