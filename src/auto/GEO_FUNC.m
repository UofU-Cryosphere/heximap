classdef GEO_FUNC
    methods(Static)

        function [vRankIdx] = RankWindows(cL, cR)

           % Function to determine best window to perform initial
           % georeferencing
           
           % Determine which window has the least gaps in its disparity
           % map
           vValid = cellfun(@(x) ...
               1 - sum(isnan(x.Disparity), 'all')/numel(x.Disparity), cL);

           % Calcuate contrast within windows (using RMS contrast)
           vConL = cellfun(@(x) ...
               std(im2single(x.Image), 1, 'all', 'omitnan'), cL);
           vConR = cellfun(@(x) ...
               std(im2single(x.Image), 1, 'all', 'omitnan'), cR);
           vCon = mean([vConL, vConR], 2);

           % Determine mean luminosity for each window (very rough
           % approximation of impact of clouds/snow/ice in image)
           vLumL = cellfun(@(x) ...
               mean(im2single(x.Image), 'all', 'omitnan'), cL);
           vLumR = cellfun(@(x) ...
               mean(im2single(x.Image), 'all', 'omitnan'), cR);
           vLum = mean([vLumL, vLumR], 2);

           % Determine reasonable range in elevation for each window
           vQuant10 = cellfun(@(x) ...
               quantile(x.Disparity, 0.1, 'all'), cL);
           vQuant90 = cellfun(@(x) ...
               quantile(x.Disparity, 0.9, 'all'), cL);
           vRange = vQuant90 - vQuant10;
           vRange = (vRange-min(vRange))/(max(vRange)-min(vRange));

           % Combine metrics to weighted average (with lower weight on
           % luminosity (higher uncertainty) and range (artificially
           % inflated due to normalization)
           mMetrics = [vValid, vCon, vLum, vRange];
           vW = [1, 1, 0.5, 0.5];
           vScore = sum(vW.*mMetrics,2)/sum(vW);

           % Get ranked indices of best scores
           [~,vRankIdx] = sort(vScore, 'descend');
           
        end


    end
end