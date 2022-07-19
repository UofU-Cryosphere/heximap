classdef STI_FUNC
    methods(Static)
        
        function corners = GetCorners(sInfoL, sInfoR)
            
            % Define corner region subset size
            nrow = 5000;
            ncol = 5000;
            
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
            
            % %             Diagnostic plot
            %             [nrow, ncol] = size(IM_subset);
            %             figure
            %             imshow(IM_subset)
            %             hold on
            %             plot(1:ncol, repmat(r_idx,ncol,1), 'r', 'LineWidth', 1)
            %             plot(repmat(c_idx, nrow,1), 1:nrow,'r', 'LineWidth', 1)
            %             hold off
            %             pause(1) % Issues with plotting require pausing
            
            
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
        
        function [objM] = Stitch(sInfoL, sInfoR, tblCorners, strExpPath)
            
            % Create mat file to store images
            [~,strN,~] = fileparts(sInfoL.Filename);
            strFile = strcat(strExpPath, strN(1:end-2), '.mat');
            objM = matfile(strFile,'Writable',true);
            
            % Deterime if figure was previously created
            if exist(strFile,'file')
                [~,strN,~] = fileparts(strFile);
                warning(['Stitched image %s already exists in %s. '...
                    'Moving on to next image.'], ...
                    strN, strExpPath);
                
            else
                
                % Get averaged edge positions
                edgeL = max([tblCorners.UL(2) tblCorners.LL(2)]);
                edgeR = min([tblCorners.UR(2) tblCorners.LR(2)]);
                
                % Get matrix of left image start/stop row/col indices
                % This has the following setup:
                % | starting row_idx | starting col_idx |
                % |  ending row_idx  |  ending col_idx  |
                mCornL = [edgeL, tblCorners.UL(1); sInfoL.Width, tblCorners.LL(1)];
                mCornR = [1, tblCorners.UR(1); edgeR, tblCorners.LR(1)];
                
                % Save left image half in the mat file
                STI_FUNC.saveLeftImageHalf(sInfoL,mCornL,objM)
                
                % Estimate transformation for right image half to match left image half
                mT = STI_FUNC.estimateTransform(sInfoR,mCornR,objM);
                
                % Save right image half in the mat file
                STI_FUNC.saveRightImageHalf(sInfoR,mCornR,objM,mT)
            end
            
            
        end
        
        function [] = saveLeftImageHalf(sInfoL,mCornL,objM)
            % Save left image half in the mat file
            
            % Make index vectors
            dNumSections = 10;
            vR = mCornL(:,2)';
            vRi = vR-mCornL(3)+1;
            vC = round(linspace(mCornL(1),mCornL(2),dNumSections));
            
            % Loop for each image section
            for i = 1:length(vC)-1
                
                % Save image section in mat file
                vCi = [vC(i) vC(i+1)]-mCornL(1)+1;
                objM.Image(vRi(1):vRi(2),vCi(1):vCi(2)) = imread(sInfoL.Filename, ...
                    'PixelRegion',{vR,[vC(i) vC(i+1)]});
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function mT = estimateTransform(sInfoR,mCornR,objM)
            
            % Get size of left image half
            [iH,iW] = size(objM,'Image');
            
            % Define pixel region boundaries
            iX = round(iW/8);
            vY = round(linspace(1,iH,8));
            iNumWin = length(vY)-1;
            
            % Initialize
            mL = []; mR = [];
            
            % Loop through each pixel region
            for i = 1:iNumWin
                
                try
                    
                    % Construct classes
                    clDetector = cv.FeatureDetector('ORB','MaxFeatures',10000);
                    clExtractor = cv.DescriptorExtractor('ORB');
                    clMatcher = cv.DescriptorMatcher('BruteForce-Hamming');
                    
                    % Read the left image, detect keypoints and descriptors, make
                    % keypoint matrix
                    mPixRegL = [iW-iX vY(i); iW vY(i+1)];
                    mI = objM.Image(mPixRegL(3):mPixRegL(4),mPixRegL(1):mPixRegL(2));
                    sKeyPoints = clDetector.detect(mI);
                    sDescriptorsL = clExtractor.compute(mI,sKeyPoints);
                    mKeyPoints = [sKeyPoints.pt];
                    mKeyPointsL = [mKeyPoints(1:2:end); mKeyPoints(2:2:end)];
                    
                    % Read the right image, detect keypoints and descriptors, make
                    % keypoint matrix
                    cPixRegR = {[vY(i) vY(i+1)]+mCornR(3) [1 iX]+mCornR(1)};
                    mI = imread(sInfoR.Filename,'PixelRegion',cPixRegR);
                    sKeyPoints = clDetector.detect(mI);
                    sDescriptorsR = clExtractor.compute(mI,sKeyPoints);
                    mKeyPoints = [sKeyPoints.pt];
                    mKeyPointsR = [mKeyPoints(1:2:end); mKeyPoints(2:2:end)];
                    clear mI sKeyPoints mKeyPoints
                    
                    % Match the keypoints
                    sMatches = clMatcher.match(sDescriptorsL,sDescriptorsR);
                    
                    % Make vectors of matched points
                    vMatchIdxL = [sMatches.queryIdx]+1;
                    vMatchIdxR = [sMatches.trainIdx]+1;
                    mPtsL = mKeyPointsL(:,vMatchIdxL);
                    mPtsR = mKeyPointsR(:,vMatchIdxR);
                    
                    % Use ransac to estimate inliers
                    sRansac.numIter = 2000;
                    sRansac.inlierDist = 1;
                    [~,lIn] = estimateTransformRansac([mPtsL;ones(1,size(mPtsL,2))], ...
                        [mPtsR;ones(1,size(mPtsR,2))], ...
                        sRansac);
                    mPtsL = mPtsL(:,lIn);
                    mPtsR = mPtsR(:,lIn);
                    
                    % Convert points to full image coordinates
                    mPtsL(1,:) = mPtsL(1,:) + mPixRegL(1) - 1;
                    mPtsL(2,:) = mPtsL(2,:) + mPixRegL(3) - 1;
                    mPtsR(1,:) = mPtsR(1,:) + cPixRegR{2}(1) - 1;
                    mPtsR(2,:) = mPtsR(2,:) + cPixRegR{1}(1) - 1;
                    
                    % Make points relative to axis with origin at lower left corner
                    mPtsL(2,:) = iH - mPtsL(2,:) + 1;
                    mPtsR(2,:) = sInfoR.Height - mPtsR(2,:) + 1;
                    
                    % Save points
                    mL = [mL mPtsL];
                    mR = [mR mPtsR];
                    
                catch
                end
                
            end
            
            % Estimate geometric transform
            sRansac.NumIter = 2000;
            sRansac.InlierDist = 1;
            [mT,~] = estimateTransformRansac([mL;ones(1,size(mL,2))], ...
                [mR;ones(1,size(mR,2))], ...
                sRansac);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [] = saveRightImageHalf(sInfoR,mCornR,objM,mT)
            
            % Get image size
            [iH,iW] = size(objM,'Image');
            
            % Spatial referencing info for image to be transformed (right half)
            sR.DeltaLon = 1;
            sR.DeltaLat = -1;
            sR.Lonlim = [1 sInfoR.Width] + [-0.5 0.5];
            sR.Latlim = [1 sInfoR.Height] + [-0.5 0.5];
            sInfoR.SpatialRef = sR;
            
            % Spatial referencing vectors for saving in mat file
            mWin = round(mT * [mCornR(2,:) 1 1; mCornR(2) 1 1 1]');
            vX = iW:max(mWin(1,:));
            vY = fliplr(1:iH);
            
            % Transform the right image half and save it in the mat file
            sParams.blockSize = 1000;
            sParams.matSource = objM.Properties.Source;
            sParams.matField = 'Image';
            sParams.matStartIndex = [min(vY) min(vX)];
            sParams.transform = mT;
            warning('off','block_process:class_match')
            grid2grid(sInfoR,vX,vY,sParams);
            warning('on','block_process:class_match')
            
        end

        function [] = Resize(objM)
            % Scales for downsampling. Note that each successive scale is relative to
            % the preceding one.
            vScale = [1/2 1/5];

            % Field names for saving downscaled images
            cFields = {'Image','Image2','Image10'};

            % Loop through each scale
            for i = 1:length(vScale)

                % Get size of higher resolution image
                [iH,iW] = size(objM,cFields{i});

                % Spatial referencing info for higher resolution image
                sR.DeltaLon = 1;
                sR.DeltaLat = -1;
                sR.Lonlim = [1 iW] + [-0.5 0.5];
                sR.Latlim = [1 iH] + [-0.5 0.5];
                sInfo.SpatialRef = sR;
                sInfo.Width = iW;
                sInfo.Height = iH;
                sInfo.matSource = objM.Properties.Source;
                sInfo.matField = cFields{i};

                % Spatial referencing vectors for downscaled image
                dS = vScale(i);
                vSz = round(dS*[iH iW]);
                vX = 1:vSz(2);
                vY = fliplr(1:vSz(1));

                % Resampling parameters
                sParams.blockSize = round(1000*dS);
                sParams.matSource = objM.Properties.Source;
                sParams.matField = cFields{i+1};
                sParams.transform = [dS 0 0 0; 0 dS 0 0; 0 0 1 0; 0 0 0 1];

                % Initialize field with uint8 data type
                objM.(cFields{i+1}) = uint8(0);

                % Resample
                warning('off','block_process:class_match')
                grid2grid(sInfo,vX,vY,sParams);
                warning('on','block_process:class_match')

            end
        end
    end
end