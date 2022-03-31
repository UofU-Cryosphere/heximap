classdef STI_FUNC
    methods(Static)
        
        function corners = GetCorners(sInfoL, sInfoR)
            
            % Define corner region subset size
            nrow = 5000;
            ncol = 5000;
            
%             % Define image halves
%             sInfoL = imfinfo(char(fullfile(im_dir, strcat(im_name, "_a.tif"))));
%             sInfoR = imfinfo(char(fullfile(im_dir, strcat(im_name, "_b.tif"))));
            
            [~,fn] = fileparts(sInfoL.Filename);
            im_name = fn(1:end-2);
            
            NW_region = imread(sInfoL.Filename,...
                'PixelRegion',{[1 nrow] [1 ncol]});
            region_corner = STI_FUNC.findCorner(NW_region, sInfoL, 'upper left');
            NW_corner = region_corner;
            
            SW_region = imread(sInfoL.Filename,...
                'PixelRegion',{[sInfoL.Height-nrow sInfoL.Height] [1 ncol]});
            region_corner = STI_FUNC.findCorner(SW_region, sInfoL, 'lower left');
            SW_corner = [(sInfoL.Height-nrow+region_corner(1)); region_corner(2)];
            
            NE_region = imread(sInfoR.Filename,...
                'PixelRegion',{[1 nrow] [sInfoR.Width-ncol sInfoR.Width]});
            region_corner = STI_FUNC.findCorner(NE_region, sInfoR, 'upper right');
            NE_corner = [region_corner(1); (sInfoR.Width-ncol+region_corner(2))];
            
            SE_region = imread(sInfoR.Filename,...
                'PixelRegion', {[sInfoR.Height-nrow sInfoR.Height] [sInfoR.Width-ncol sInfoR.Width]});
            region_corner = STI_FUNC.findCorner(SE_region, sInfoR, 'lower right');
            SE_corner = [(sInfoR.Height-nrow+region_corner(1)); (sInfoR.Width-ncol+region_corner(2))];
            
            
            % W_idx = round(mean([NW_corner(2) SW_corner(2)]));
            N_diff = abs(NW_corner(1) - NE_corner(1));
            S_diff = abs(SW_corner(1) - SE_corner(1));
            E_diff = abs(NE_corner(2) - SE_corner(2));
            W_diff = abs(NW_corner(2) - SW_corner(2));
            
            % Add warning for edge values that are too different
            px_threshold = 30;
            if N_diff > px_threshold
                warning(['Detected upper image edge differs between halves '...
                'by more than %i pixels for %s. '... 
                'This is likely caused by the scannned image edge not '...
                'properly aligned to true horiztonal, but could also '
                'indicate an issue with corner selection.'],...
                px_threshold, im_name);
            end
            if S_diff > px_threshold
                warning(['Detected lower image edge differs between halves '...
                'by more than %i pixels for %s. '... 
                'This is likely caused by the scannned image edge not '...
                'properly aligned to true horiztonal, but could also '
                'indicate an issue with corner selection.'],...
                px_threshold, im_name);
            end
            if W_diff > px_threshold
                warning(['Detected left image edge differs '...
                'by more than %i pixels for %s. '... 
                'This is likely caused by the scannned image edge not '...
                'properly aligned to true vertical, but could also '
                'indicate an issue with corner selection.'],...
                px_threshold, im_name);
            end
            if E_diff > px_threshold
                warning(['Detected right image edge differs '...
                'by more than %i pixels for %s. '... 
                'This is likely caused by the scannned image edge not '...
                'properly aligned to true vertical, but could also '
                'indicate an issue with corner selection.'],...
                px_threshold, im_name);
            end
            
            corners = table(...
                NW_corner, NE_corner, SW_corner, SE_corner,...
                'VariableNames', {'UL', 'UR', 'LL', 'LR'},...
                'RowNames', {'row_idx', 'col_idx'});
            
        end
        
        function corner_loc = findCorner(IM_subset, sInfo, corner_name)
            
            % Define defaults
            gauss_sigma = 30;
            text_threshold = 250;
            
            
            % Set pixels above threshold (should represent text) to zero
            IM_subset(IM_subset > text_threshold) = 0;
            
            % Convert to logical and apply Gaussian convolution
            dat_log = double(logical(IM_subset));
            dat_blur = imgaussfilt(dat_log, gauss_sigma);
            
            [Gmag,~] = imgradient(dat_blur);
%             [R_mag,Rm_idx] = max(Gmag);
%             R_weights = R_mag.^2/sum(R_mag.^2);
%             r_idx = round(sum(R_weights .* Rm_idx));
            
            % Sum row-wise and col-wise gradients for each pixel
            px_sum = sum(Gmag,1) + sum(Gmag,2);
            
%             % Diagnostic plot
%             figure
%             imagesc(px_sum)
%             colorbar
%             title('Row-wise and col-wise gradient summation for pixels')
%             xlabel('Column number')
%             ylabel('Row number')
%             pause(1) % Issues with plotting require pausing
            
            % Find global maximum position (in the future could consider
            % this as an a optimaztion to return local maxima as well)
            [px_max,px_idx] = max(px_sum(:));
            [r_idx, c_idx] = ind2sub(size(px_sum), px_idx);
            
            % Check if other distant positions have a large value as well
            idx_check = find(px_sum > 0.90*px_max);
            [r_check, c_check] = ind2sub(size(px_sum), idx_check);
            px_dist = sqrt((r_check-r_idx).^2 + (c_check-c_idx).^2);
            px_power = (px_sum(idx_check)/px_max) .* px_dist;
            
            % Get filename (for warning messages)
            [~,fn] = fileparts(sInfo.Filename);

            if any(px_power > 3*gauss_sigma)
                warning(['Additional possible positions for the %s corner '...
                    'were detected in image %s. '...
                    'Consider manually verifying the %s corner position.'], ...
                    corner_name, fn, corner_name)
            end
            
            % Refine search based on detected corner using unblurred image
            dat_fine = IM_subset(r_idx-2*gauss_sigma:r_idx+2*gauss_sigma,...
                c_idx-2*gauss_sigma:c_idx+2*gauss_sigma);
            [Gmag,~] = imgradient(dat_fine);
            px_sum = sum(Gmag,1) + sum(Gmag,2);
            [~,px_idx] = max(px_sum(:));
            [r_fine, c_fine] = ind2sub(size(px_sum), px_idx);
            
            r_idx = r_idx-2*gauss_sigma + r_fine-1;
            c_idx = c_idx-2*gauss_sigma + c_fine-1;
            
            % Output corner position
            corner_loc = [r_idx; c_idx];
            
            % Diagnostic plot            
            [nrow, ncol] = size(IM_subset);
            figure
            imshow(IM_subset)
            hold on
            plot(1:ncol, repmat(r_idx,ncol,1), 'r', 'LineWidth', 1)
            plot(repmat(c_idx, nrow,1), 1:nrow,'r', 'LineWidth', 1)
            hold off
            pause(1) % Issues with plotting require pausing
            
            
%             % Find max gradient in columns
%             [c_pks,c_locs,~,~] = findpeaks(sum(Gmag,1), ...
%                 'SortStr', 'descend', 'NPeaks', 10, ...
%                 'MinPeakDistance', gauss_sigma);
% %             [~, c_idx] = max(abs(diff(sum(dat_blur,1))));
%             c_idx = c_locs(1);
%             
%             % Estimate dis-similarity in col gradient peaks
%             c_pk_range = max(c_pks) - min(c_pks);
%             c_dissim = (c_pks-min(c_pks))/c_pk_range .* abs(c_locs-c_idx);
%             
%             
%             % Find max gradient in rows
%             [r_pks,r_locs,~,~] = findpeaks(sum(Gmag,2), ...
%                 'SortStr', 'descend', 'NPeaks', 10, ...
%                 'MinPeakDistance', gauss_sigma);
% %             [~, r_idx] = max(abs(diff(sum(dat_blur,2))));
%             r_idx = r_locs(1);
%             
%             % Estimate dis-similarity in row gradient peaks
%             r_pk_range = max(r_pks) - min(r_pks);
%             r_dissim = (r_pks-min(r_pks))/r_pk_range .* abs(r_locs-r_idx);
%             
%             % Get filename (for warning messages)
%             [~,fn] = fileparts(sInfo.Filename);
%             
%             % Add warning if errors are greater than scaled gaussian blur
%             if any(c_dissim > gauss_sigma)
%                 warning(['Additional possible values for %s corner column '...
%                     'position were detected for image %s. '...
%                     'Consider manually verifying %s corner position.'], ...
%                     corner_name, fn, corner_name)
%             end
%             if any(r_dissim > gauss_sigma)
%                 warning(['Additional possible values for %s corner row '...
%                     'position were detected for image %s. '...
%                     'Consider manually verifying %s corner position.'], ...
%                     corner_name, fn, corner_name)
%             end
%             
%             % Add warning if range in peak values is too small
%             if c_pks(end) > 0.90*c_pks(1)
%                 warning(['Limited range in column gradient peaks in %s. '...
%                     'Detected %s corner could be unreliable.']', ...
%                     fn, corner_name)
%             end
%             if r_pks(end) > 0.90*r_pks(1)
%                 warning(['Limited range in row gradient peaks in %s. '...
%                     'Detected %s corner could be unreliable.']', ...
%                     fn, corner_name)
%             end
            
            
        end
        
    end
    
end