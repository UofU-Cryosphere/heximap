classdef STI_FUNC
    methods(Static)
        
        function corners = GetCorners(im_dir, im_name)
            
            % Define corner region subset size
            nrow = 5000;
            ncol = 5000;
            
            % Define image halves
            sInfoL = imfinfo(char(fullfile(im_dir, strcat(im_name, "_a.tif"))));
            sInfoR = imfinfo(char(fullfile(im_dir, strcat(im_name, "_b.tif"))));
            
            
            NW_region = imread(sInfoL.Filename,...
                'PixelRegion',{[1 nrow] [1 ncol]});
            region_corner = STI_FUNC.findCorner(NW_region, sInfoL);
            NW_corner = region_corner;
            
            SW_region = imread(sInfoL.Filename,...
                'PixelRegion',{[sInfoL.Height-nrow sInfoL.Height] [1 ncol]});
            region_corner = STI_FUNC.findCorner(SW_region, sInfoL);
            SW_corner = [(sInfoL.Height-nrow+region_corner(1)); region_corner(2)];
            
            NE_region = imread(sInfoR.Filename,...
                'PixelRegion',{[1 nrow] [sInfoR.Width-ncol sInfoR.Width]});
            region_corner = STI_FUNC.findCorner(NE_region, sInfoR);
            NE_corner = [region_corner(1); (sInfoR.Width-ncol+region_corner(2))];
            
            SE_region = imread(sInfoR.Filename,...
                'PixelRegion', {[sInfoR.Height-nrow sInfoR.Height] [sInfoR.Width-ncol sInfoR.Width]});
            region_corner = STI_FUNC.findCorner(SE_region, sInfoR);
            SE_corner = [(sInfoR.Height-nrow+region_corner(1)); (sInfoR.Width-ncol+region_corner(2))];
            
            
            % W_idx = round(mean([NW_corner(2) SW_corner(2)]));
            N_diff = abs(NW_corner(1) - NE_corner(1));
            S_diff = abs(SW_corner(1) - SE_corner(1));
            E_diff = abs(NE_corner(2) - SE_corner(2));
            W_diff = abs(NW_corner(2) - SW_corner(2));
            
            % Add warning for edge values that are too different
            px_threshold = 25;
            if N_diff > px_threshold
                warning(['Detected upper image edge differs between halves '...
                'by more than %i pixels for %s'],...
                px_threshold, im_name);
            end
            if S_diff > px_threshold
                warning(['Detected lower image edge differs between halves '...
                'by more than %i pixels for %s'],...
                px_threshold, im_name);
            end
            if W_diff > px_threshold
                warning(['Detected left image edge differs '...
                'by more than %i pixels for %s'],...
                px_threshold, im_name);
            end
            if E_diff > px_threshold
                warning(['Detected right image edge differs '...
                'by more than %i pixels for %s'],...
                px_threshold, im_name);
            end
            
            corners = table(...
                NW_corner, NE_corner, SW_corner, SE_corner,...
                'VariableNames', {'UL', 'UR', 'LL', 'LR'},...
                'RowNames', {'row_idx', 'col_idx'});
            
        end
        
        function corner_loc = findCorner(IM_subset, sInfo)
            
            % Define defaults
            gauss_sigma = 100;
            text_threshold = 250;
            
            
            % Set pixels above threshold (should represent text) to zero
            IM_subset(IM_subset > text_threshold) = 0;
            
            % Convert to logical and apply Gaussian convolution
            dat_log = double(logical(IM_subset));
            dat_blur = imgaussfilt(dat_log, gauss_sigma);
            
            % Find max gradient in columns
            [c_pks,c_locs,~,c_P] = findpeaks(abs(diff(sum(dat_blur,1))), ...
                'SortStr', 'descend', 'NPeaks', 10, ...
                'MinPeakDistance', gauss_sigma/2);
%             [~, c_idx] = max(abs(diff(sum(dat_blur,1))));
            c_idx = c_locs(1);
            
            % Find max gradient in rows
            [r_pks,r_locs,~,r_P] = findpeaks(abs(diff(sum(dat_blur,2))), ...
                'SortStr', 'descend', 'NPeaks', 10, ...
                'MinPeakDistance', gauss_sigma/2);
%             [~, r_idx] = max(abs(diff(sum(dat_blur,2))));
            r_idx = r_locs(1);
            
            
            
            % Add warning if another peak is of similar magnitude
            if c_pks(2) > 0.90*c_pks(1)
                [~,fn] = fileparts(sInfo.Filename);
                warning('Similar gradient values for columns in %s', fn)
            end
            if r_pks(2) > 0.90*r_pks(1)
                [~,fn] = fileparts(sInfo.Filename);
                warning('Similar gradient values for rows in %s', fn)
            end
            
%             % Diagnostic plot
%             [nrow, ncol] = size(IM_subset);
%             figure
%             imshow(IM_subset)
%             hold on
%             plot(1:ncol, repmat(r_idx,ncol,1), 'r', 'LineWidth', 1)
%             plot(repmat(c_idx, nrow,1), 1:nrow,'r', 'LineWidth', 1)
%             hold off
%             pause(1) % Issues with plotting require pausing
            
            corner_loc = [r_idx; c_idx];
            
        end
        
    end
    
end