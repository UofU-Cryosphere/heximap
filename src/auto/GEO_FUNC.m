classdef GEO_FUNC
    methods(Static)

        function [cL, cR] = GetMatFiles(strWinPath)

            % Get paths to mat files
            cL = cellfun(@(x) matfile(x,'Writable',true),getFiles(strWinPath, ...
                'Left.mat'),'Uni',0);
            cR = cellfun(@(x) matfile(x,'Writable',true),getFiles(strWinPath, ...
                'Right.mat'),'Uni',0);
            
            % Make sure corresponding left and right images exist for each window
            if numel(cL) ~= numel(cR)
                error('Number of left and right .mat files is not equal.')
            end


        end

    end
end