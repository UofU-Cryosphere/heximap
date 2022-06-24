classdef EXT_FUNC
    methods(Static)
        function [] = ControlPoints(objM, sMeta)
            % Put corner coordinates in correct order. This assumes the Hexagon
            % image was scannned with left side = north, right side = south
            % i.e. with the following layout:
            % NE Lon, NE Lat; SW Lon, SW Lat; SE Lon, SE Lat; NW Lon, NW Lat
            mPtsWld = str2double(...
                {sMeta.NE_Corne_3, sMeta.NE_Corne_2; ...
                sMeta.SW_Corne_3, sMeta.SW_Corne_2; ...
                sMeta.SE_Corne_3, sMeta.SE_Corne_2; ...
                sMeta.NW_Corne_3, sMeta.NW_Corne_2});

            % Define image size and corners. Also note that the y-coordinates are
            % made negative to correspond with right-handed world coordinate system
            [iH,iW] = size(objM,'Image');
            mPtsImg = [1 1; iW iH; iW 1; 1 iH];
            mPtsImg(:,2) = -mPtsImg(:,2);

            % Compute spatial transformation from control point pairs. Can ignore
            % condition number warning.
            warning('off','images:maketform:conditionNumberofAIsHigh')
            sT = cp2tform(mPtsImg,mPtsWld,'nonreflective similarity');
            warning('on','images:maketform:conditionNumberofAIsHigh')

            % Save spatial transformation structure and corner coordinates to mat
            % file
            objM.SpatialTrans = affine2d(sT.tdata.T);
            objM.CornerGCPs = mPtsWld;
        end

        function [vOrder] = SortImages(cM)
            % Get height and width of images
            [cH,cW] = cellfun(@(x) size(x,'Image'),cM,'UniformOutput',0);

            % Get spatial transformation structures
            cT = cellfun(@(x) x.SpatialTrans,cM,'UniformOutput',0);

            % Loop for each image
            mCenters = zeros(numel(cM),2);
            mDir = zeros(numel(cM),2);
            for i = 1:numel(cM)

                % Define corners. Note y-axis sign is reversed
                mPtsImg = [1 1; cW{i} cH{i}; cW{i} 1; 1 cH{i}];
                mPtsImg(:,2) = -mPtsImg(:,2);

                % Transform corners to world coordinates
                mPtsWld = transformPointsForward(cT{i},mPtsImg);

                % Compute image center
                mCenters(i,:) = mean(mPtsWld);

                % Compute vector from corner 1 to corner 3
                mDir(i,:) = mPtsWld(3,:) - mPtsWld(1,:);

            end

            % Determine image order by rotating the images so image x-axes are
            % horizontal in the world coordinate system
            vPosX = mean([mDir zeros(size(mDir,1),1)]) / ...
                norm(mean([mDir zeros(size(mDir,1),1)]));
            mR = vecRotMat(vPosX,[1 0 0]);
            mCentersR = (mR * [mCenters zeros(size(mCenters,1),1)]')';
            [~,vOrder] = sortrows(mCentersR,1);

        end

        function [] = InitTrans(cM)
            % Find and match corresponding points (using ORB descriptors), use
            % these points to estimate the fundamental matrix (using 8-point
            % algorithm), use the fundamental and intrinsic matrices to compute the
            % essential matrix, estimate the relative camera pose matrices
            % directly from the essential matrix (using SVD decomposition and choosing
            % the correct solution out of 4 possible), then estimate rough homography
            % transformations for the full Hexagon images (using nonlinear least
            % squares solver).  These solutions are quite rough, and are later refined
            % seperately for each processing window (using more accurate point
            % correspondences during stereo rectification and bundle adjustment).

            % Get mat file objects
            objML = cM{1};
            objMR = cM{2};

            % Find where images overlap
            vSz1 = size(objML,'Image') - 10;
            vSz2 = size(objMR,'Image') - 10;
            mWinL = [1 vSz1(2);1 vSz1(1)]';
            mWinL = flipud(mWinL([1 3; 1 4; 2 4; 2 3]));
            mWinR = [1 vSz2(2);1 vSz2(1)]';
            mWinR = flipud(mWinR([1 3; 1 4; 2 4; 2 3]));
            mWinL = transformPointsForward(objML.SpatialTrans, ...
                [mWinL(:,1) mWinL(:,2)*-1]);
            mWinR = transformPointsForward(objMR.SpatialTrans, ...
                [mWinR(:,1) mWinR(:,2)*-1]);
            [vX,vY] = polybool('Intersection',mWinL(:,1),mWinL(:,2), ...
                mWinR(:,1),mWinR(:,2));
            mWinL = transformPointsInverse(objML.SpatialTrans,[vX vY]);
            mWinR = transformPointsInverse(objMR.SpatialTrans,[vX vY]);
            mWinL(:,2) = -mWinL(:,2);
            mWinR(:,2) = -mWinR(:,2);
            mWinL = round([min(mWinL);max(mWinL)]);
            mWinR = round([min(mWinR);max(mWinR)]);
            mWinL(3:4) = [1 vSz1(1)]';
            mWinR(3:4) = [1 vSz2(1)]';

            % Make sure windows don't extend past image boundaries
            mWinL(mWinL < 10) = 10;
            mWinR(mWinR < 10) = 10;
            mB = fliplr([vSz1;vSz1]);
            mWinL(mWinL > mB) = mB(mWinL > mB);
            mB = fliplr([vSz2;vSz2]);
            mWinR(mWinR > mB) = mB(mWinR > mB);

            % Compute point correspondences
            [mPts1,mPts2] = EXT_FUNC.getPointCorrespondences(objML,objMR,mWinL,mWinR);

            % Initialize
            iNumTries = 50;
            cPose2 = cell(1,iNumTries);

            % Try computing pose multiple times, then save best solution
            for i = 1:iNumTries

                % Compute intrinsic and pose matrices
                [mPts1s,mPts2s,mK1,mK2,mPose1,mPose2] = EXT_FUNC.computeMatrices(objML, ...
                    objMR,mPts1,mPts2);

                % Save right camera pose output
                cPose2{i} = mPose2;

            end

            % Save right camera pose solution with z-coordinate closest to zero
            lE = cell2mat(cellfun(@(x) all(x(:) == 0),cPose2,'Uni',0));
            cPose2 = cPose2(~lE);
            [~,iIdx] = min(cell2mat(cellfun(@(x) abs(x(3,4)),cPose2,'Uni',0)));
            mPose2 = cPose2{iIdx};

%             % Save ORB point figures
%             hF = figure;
%             imagesc(objML.Image10),colormap(bone(256)),axis equal,axis off,hold on
%             scatter(mPts1s(:,1)/10,mPts1s(:,2)/10,[],[1,0.5,0])
%             title('Left Image ORB Matches')
%             saveas(hF,[strHexPath 'LeftFullImageMatches'],'fig')
%             close(hF)
%             hF = figure;
%             imagesc(objMR.Image10),colormap(bone(256)),axis equal,axis off,hold on
%             scatter(mPts2s(:,1)/10,mPts2s(:,2)/10,[],[1,0.5,0])
%             title('Right Image ORB Matches')
%             saveas(hF,[strHexPath 'RightFullImageMatches'],'fig')
%             close(hF)

%             % Plot camera poses and save figure
%             hF = plotCameraPoses(inv([mPose1; 0 0 0 1]),inv([mPose2; 0 0 0 1]));
%             saveas(hF,[strHexPath 'RelativeCameraPoses'],'fig');
%             close(hF)

            % Compute rough homographies. More refined solutions will be computed later
            % for each processing window.
            [mH1,mH2] = EXT_FUNC.homographies(mK1,mK2,mPts1s,mPts2s);

            % Save intrinsic and pose matrices in mat file
            objML.IntrinsicMatrix = mK1;
            objML.LeftPoseMatrix = mPose1;
            objML.RightPoseMatrix = mPose2;

            % Save homography matrices in mat file
            objML.LeftHomography = mH1;
            objML.RightHomography = mH2;
        end
        
        function [mPts1,mPts2] = getPointCorrespondences(objML,objMR,mWinL,mWinR)
            % Compute ORB correspondences

            % Read overlapping regions of downsized images
            mWinL = round(mWinL / 10);
            mWinR = round(mWinR / 10);
            mI1 = objML.Image10(mWinL(3):mWinL(4),mWinL(1):mWinL(2));
            mI2 = objMR.Image10(mWinR(3):mWinR(4),mWinR(1):mWinR(2));

            % Filter the images
            sOpt.histmatch = true;
            sOpt.adapthisteq = true;
            sOpt.wiener2 = true;
            [mI1,mI2] = EXT_FUNC.FilterImages(mI1,mI2,0,sOpt);

            % Construct classes
            clDetector = cv.FeatureDetector('ORB','MaxFeatures',30000);
            clExtractor = cv.DescriptorExtractor('ORB');
            clMatcher = cv.DescriptorMatcher('BruteForce-Hamming');

            % Detect keypoints and descriptors
            sKeyPoints1 = clDetector.detect(mI1);
            sKeyPoints2 = clDetector.detect(mI2);
            sDescriptors1 = clExtractor.compute(mI1,sKeyPoints1);
            sDescriptors2 = clExtractor.compute(mI2,sKeyPoints2);

            % Make matrices of keypoints
            mKeyPoints1 = [sKeyPoints1.pt];
            mKeyPoints2 = [sKeyPoints2.pt];
            mKeyPoints1 = [mKeyPoints1(1:2:end); mKeyPoints1(2:2:end)]';
            mKeyPoints2 = [mKeyPoints2(1:2:end); mKeyPoints2(2:2:end)]';

            % Make index vectors
            vSz = size(mI1);
            vIdxX = round(linspace(1,vSz(2),5));
            vIdxY = round(linspace(1,vSz(1),5));

            % Initialize point matrices
            mPts1 = [];
            mPts2 = [];

            % Loop for each sub-window
            for iY = 1:length(vIdxY)-1
                for iX = 1:length(vIdxX)-1

                    try

                        % Make sub-array of left image points within current window
                        lIn1 = mKeyPoints1(:,1) > vIdxX(iX) & ...
                            mKeyPoints1(:,1) < vIdxX(iX+1) &...
                            mKeyPoints1(:,2) > vIdxY(iY) & ...
                            mKeyPoints1(:,2) < vIdxX(iY+1);
                        sDescriptors1s = sDescriptors1(lIn1,:);

                        % Make sub-array of right image points within current window
                        lIn2 = mKeyPoints2(:,1) > vIdxX(iX) & ...
                            mKeyPoints2(:,1) < vIdxX(iX+1) & ...
                            mKeyPoints2(:,2) > vIdxY(iY) & ...
                            mKeyPoints2(:,2) < vIdxX(iY+1);
                        sDescriptors2s = sDescriptors2(lIn2,:);

                        % Match the keypoints. Save 100 strongest matches
                        sMatches = clMatcher.match(sDescriptors1s,sDescriptors2s);
                        [~,vIdx] = sort([sMatches.distance],'ascend');
                        sMatches = sMatches(vIdx(1:100));

                        % Make vectors of matched points
                        vMatchIdx1 = [sMatches.queryIdx]+1;
                        vMatchIdx2 = [sMatches.trainIdx]+1;

                        % Convert left image matched points to full array indices
                        vIdx1 = 1:length(lIn1);
                        vIdx1 = vIdx1(lIn1);
                        vMatchIdx1 = vIdx1(vMatchIdx1);

                        % Convert right image matched points to full array indices
                        vIdx2 = 1:length(lIn2);
                        vIdx2 = vIdx2(lIn2);
                        vMatchIdx2 = vIdx2(vMatchIdx2);

                        % Save matched keypoints
                        mPts1 = [mPts1; mKeyPoints1(vMatchIdx1,:)];
                        mPts2 = [mPts2; mKeyPoints2(vMatchIdx2,:)];

                    catch
                    end
                end
            end

            % Convert points to full image coordinates
            mPts1(:,1) = mWinL(1) + mPts1(:,1) - 1;
            mPts2(:,1) = mWinR(1) + mPts2(:,1) - 1;
            mPts1(:,2) = mWinL(3) + mPts1(:,2) - 1;
            mPts2(:,2) = mWinR(3) + mPts2(:,2) - 1;

        end

        function [mPts1,mPts2,mK1,mK2,mPose1,mPose2] = computeMatrices( ...
                objML,objMR,mPts1,mPts2)
            % Compute fundamental, intrinsic, and pose matrices

            % Choose random subset of points
            vIdx = randperm(size(mPts1,1), ...
                min([round(4/5*size(mPts1,1)) size(mPts1,1)]));
            mPts1 = mPts1(vIdx,:);
            mPts2 = mPts2(vIdx,:);

            % Remove outliers using ransac
            [~,lIn] = cv.findFundamentalMat(num2cell(mPts1,2)',num2cell(mPts2,2)', ...
                'Method','Ransac','RansacReprojThreshold',1,'Confidence',0.99);
            mPts1 = mPts1(logical(lIn),:);
            mPts2 = mPts2(logical(lIn),:);

            % Convert points to full scale
            mPts1 = mPts1*10;
            mPts2 = mPts2*10;

            % Compute fudamental matrix using 8 point algorithm
            [mF,~] = cv.findFundamentalMat(num2cell(mPts1,2)',num2cell(mPts2,2)', ...
                'Method','8Point');

            % Get data from mat files
            dF = objML.FocalLengthPixels;
            vPP1 = objML.PrincipalPointPixels;
            vPP2 = objMR.PrincipalPointPixels;

            % Define left camera intrinsic matrix
            mK1 = [dF  0 vPP1(1)
                0 dF vPP1(2)
                0  0  1];

            % Define right camera intrinsic matrix
            mK2 = [dF  0 vPP2(1)
                0 dF vPP2(2)
                0  0  1];

            % Compute relative camera pose matrices directly from the essential matrix
            % via SVD decomposition
            [mPose1,mPose2] = EXT_FUNC.decomposeEssentialMatrix(mK1,mK2,mF,mPts1,mPts2);

        end

        function [mPose1,mPose2] = decomposeEssentialMatrix(mK1,mK2,mF,mPts1,mPts2)
            % Get relative camera pose matrices from essential matrix using SVD. Four
            % solutions are possible, but in only one solution will the triangulated
            % points have positive depth.

            % Choose 500 random points
            mPts1 = [mPts1 ones(size(mPts1,1),1)]';
            vIdx = randperm(length(mPts1),min([500 length(mPts1)]));
            mPts1 = mPts1(:,vIdx);
            mPts2 = [mPts2 ones(size(mPts2,1),1)]';
            mPts2 = mPts2(:,vIdx);

            % Compute essential matrix
            mE = mK2' * mF * mK1;

            % Perform SVD on essential matrix
            [mU,~,mV] = svd(mE);
            mE = mU * diag([1 1 0]) * mV';
            [mU,~,mV] = svd(mE);

            % Define W matrix
            mW = [0 -1 0; 1 0 0; 0 0 1];

            % Compute 4 possible pose solutions
            cPoseR{1} = [mU * mW * mV' mU(:,3)];
            cPoseR{2} = [mU * mW * mV' -mU(:,3)];
            cPoseR{3} = [mU * mW' * mV' mU(:,3)];
            cPoseR{4} = [mU * mW' * mV' -mU(:,3)];

            % Loop for each possible solution
            mPose1 = [eye(3) zeros(3,1)];
            mPose2 = zeros(3,4);
            for iS = 1:4

                % Triangulate points with current solution
                mPts3D = triangulate(mPts1,mPts2,mPose1,cPoseR{iS},mK1,mK2,false);

                % Transform triangulated points to camera coordinate systems
                mPts1h = mPose1 * mPts3D;
                mPts2h = cPoseR{iS} * mPts3D;

                % Compute percentage of positive depth points in both cameras
                dPctPos1 = sum(mPts1h(3,:) > 0) / size(mPts1h,2);
                dPctPos2 = sum(mPts2h(3,:) > 0) / size(mPts2h,2);

                % Save correct pose solution with positive depth in both cameras and a
                % determinate of +1
                if dPctPos1 > 0.95 && dPctPos2 > 0.95 && ...
                        round(det(cPoseR{iS}(1:3,1:3))) == 1

                    % Save correct solution
                    mPose2(:,:,1) = cPoseR{iS};

                end
            end
        end

        function hF = plotCameraPoses(mP1,mP2)
            % Plot camera axes representing poses

            % Plot the left camera pose
            hF = figure;
            dScale = sum((mP1(1:3,4)-mP2(1:3,4)).^2)/2;
            vX1 = repmat(mP1(1,4),1,3);
            vY1 = repmat(mP1(2,4),1,3);
            vZ1 = repmat(mP1(3,4),1,3);
            vU1 = mP1(1,1:3);
            vV1 = mP1(2,1:3);
            vW1 = mP1(3,1:3);
            quiver3(vX1,vY1,vZ1,vU1,vV1,vW1,dScale),hold on
            vU1 = vU1 * dScale/2;
            vV1 = vV1 * dScale/2;
            vW1 = vW1 * dScale/2;
            text(vX1(1),vY1(1),vZ1(1),'Left Camera')
            text(vX1(1)+vU1(1),vY1(1)+vV1(1),vZ1(1)+vW1(1),'x');
            text(vX1(2)+vU1(2),vY1(2)+vV1(2),vZ1(2)+vW1(2),'y');
            text(vX1(3)+vU1(3),vY1(3)+vV1(3),vZ1(3)+vW1(3),'z');

            % Plot the right camera pose
            vX2 = repmat(mP2(1,4),1,3);
            vY2 = repmat(mP2(2,4),1,3);
            vZ2 = repmat(mP2(3,4),1,3);
            vU2 = mP2(1,1:3);
            vV2 = mP2(2,1:3);
            vW2 = mP2(3,1:3);
            quiver3(vX2,vY2,vZ2,vU2,vV2,vW2,dScale),hold on
            vU2 = vU2 * dScale/2;
            vV2 = vV2 * dScale/2;
            vW2 = vW2 * dScale/2;
            text(vX2(1),vY2(1),vZ2(1),'Right Camera')
            text(vX2(1)+vU2(1),vY2(1)+vV2(1),vZ2(1)+vW2(1),'x');
            text(vX2(2)+vU2(2),vY2(2)+vV2(2),vZ2(2)+vW2(2),'y');
            text(vX2(3)+vU2(3),vY2(3)+vV2(3),vZ2(3)+vW2(3),'z');

            % Define plot properties
            xlabel('x'), ylabel('y'), zlabel('z')
            xlim([min([vX1(1) vX2(1)])-dScale max([vX1(1) vX2(1)])+dScale])
            ylim([min([vY1(1) vY2(1)])-dScale max([vY1(1) vY2(1)])+dScale])
            zlim([min([vZ1(1) vZ2(1)])-dScale max([vZ1(1) vZ2(1)])+dScale])
            axis equal
            title('Hexagon Camera Poses')

        end

        function [mH1,mH2] = homographies(mK1,mK2,mPts1,mPts2)
            % Use nonlinear optimization to compute initial homography transformations

            % Center the points
            vMin = min(mPts1);
            vMax = max(mPts1);
            vOrigin1 = vMin + (vMax - vMin)/2;
            mPts1 = mPts1 - repmat(vOrigin1,size(mPts1,1),1);
            vMin = min(mPts2);
            vMax = max(mPts2);
            vOrigin2 = vMin + (vMax - vMin)/2;
            mPts2 = mPts2 - repmat(vOrigin2,size(mPts2,1),1);

            % Make points homogeneous
            mPts1 = [mPts1 ones(size(mPts1,1),1)]';
            mPts2 = [mPts2 ones(size(mPts2,1),1)]';

            % Set initial paramaters and solver options
            vVar = zeros(1,5);
            sOpt = optimset('TolFun',1E-6,'TolX',1E-6,'display','iter');

            % Define scaling
            dRxL = -90; dRxU = 90;
            dRyL = -1; dRyU = 1;
            dRzL = -90; dRzU = 90;
            mBnd = repmat(vVar,2,1) + [dRyL dRzL dRxL dRyL dRzL
                dRyU dRzU dRxU dRyU dRzU];

            % Scale the parameters
            vVar = zeros(size(vVar));
            for iS = 1:length(vVar)
                vVar(iS) = (vVar(iS) - mBnd(1,iS)) / (mBnd(2,iS) - mBnd(1,iS));
            end

            % Minimize cost function using nonlinear least-squares
            vVar = lsqnonlin(@(x) EXT_FUNC.optError(x,mK1,mK2,mPts1,mPts2,mBnd), ...
                vVar,zeros(1,5),ones(1,5),sOpt);

            % Make optimized homography matrices
            [mH1,mH2] = EXT_FUNC.makeHomographies(vVar,mK1,mK2,mPts1,mPts2,mBnd);

            % Make translation matrices
            mT1 = eye(3);
            mT2 = eye(3);
            mT1(1:2,3) = -vOrigin1;
            mT2(1:2,3) = -vOrigin2;

            % Apply translation matrices
            mH1 = mH1 * mT1;
            mH2 = mH2 * mT2;

        end

        function vError = optError(vVar,mK1,mK2,mPts1,mPts2,mBnd)

            % Make homography matrices
            [mH1o,mH2o] = EXT_FUNC.makeHomographies(vVar,mK1,mK2,mPts1,mPts2,mBnd);

            % Transform the points
            mPts1 = mH1o * mPts1;
            mPts2 = mH2o * mPts2;
            mPts1 = mPts1 ./ repmat(mPts1(3,:),3,1);
            mPts2 = mPts2 ./ repmat(mPts2(3,:),3,1);

            % Compute error
            vError = mPts1(2,:) - mPts2(2,:);

        end

        function [mH1,mH2] = makeHomographies(vVar,mK1,mK2,mPts1,mPts2,mBnd)

            % Descale the parameters
            for iF = 1:length(vVar)
                vVar(iF) = mBnd(1,iF) + vVar(iF) * (mBnd(2,iF) - mBnd(1,iF));
            end

            % Make rotation matrices
            mR1 = makeRotMat(0,vVar(1),vVar(2),'deg');
            mR2 = makeRotMat(vVar(3),vVar(4),vVar(5),'deg');

            % Make homography matrices
            dX = mean(mPts1(1,:)) - mean(mPts2(1,:));
            mH1 = mK1 * (mR1 / mK1); mH1(1,3) = 0;
            mH2 = mK2 * (mR2 / mK2); mH2(1,3) = dX;

        end
        
        function [mI1,mI2] = FilterImages(mI1,mI2,iE,sOpt)
            % Match image histograms, apply locally adaptive contrast and noise
            % filters

            % Create masks for empty pixels
            lM1 = mI1 ~= iE;
            lM2 = mI2 ~= iE;

            if sOpt.histmatch

                % Match image histograms
                mI2(lM2) = imhistmatch(mI2(lM2),mI1(lM1),256);

            end

            if sOpt.adapthisteq

                % Apply locally adaptive histogram equalization
                mI1 = adapthisteq(mI1,'ClipLimit',0.03,'NumTiles',[20 20]);
                mI2 = adapthisteq(mI2,'ClipLimit',0.03,'NumTiles',[20 20]);

            end

            if sOpt.wiener2

                % Apply noise removal filter
                mI1 = wiener2(mI1,[3 3]);
                mI2 = wiener2(mI2,[3 3]);

            end

            % Reassign empty pixels
            mI1(~lM1) = iE;
            mI2(~lM2) = iE;
        end

        function [cWindow] = DefineRegions(cM, rgi_dat)

            cWindow = {cM, rgi_dat};


        end

        function [] = DisparityLoop(cM,strExPath,cWindow,strRes,iBlkSz)
            % Rectify stereo images and compute disparity maps for each window

%             % Warn user if files exist from a previous run
%             cFileL = getFiles(strExPath,'Left.mat');
%             cFileR = getFiles(strExPath,'Right.mat');
%             if ~isempty(cFileL) || ~isempty(cFileR)
%                 strQ = questdlg([ ...
%                     'Existing ''Left.mat'' and ''Right.mat'' files ' ...
%                     'in folder will cause errors. ' ...
%                     'Please move them, then click continue.'], ...
%                     '','Continue','Cancel','Cancel');
%                 if ~strcmpi(strQ,'Continue')
%                     error('Process cancelled by user.')
%                 end
%             end

            % Initialize
            iNumWin = numel(cWindow);

            % Loop through each window
            for iW = 1:iNumWin

                try

                    % Update command window
                    disp(['processing window ' num2str(iW) '...'])

                    % Delete old directory and make new one for saving data
                    strSavePath = fullfile(strExPath, ...
                        strcat("Win-", num2str(iW)));
                    if exist(strSavePath,'dir')
                        warning('off','MATLAB:RMDIR:RemovedFromPath');
                        rmdir(strSavePath,'s');
                        warning('on','MATLAB:RMDIR:RemovedFromPath');
                    end
                    mkdir(strSavePath)

                    % Get full Hexagon image mat files
                    objML = cM{1};
                    objMR = cM{2};

                    % Make mat file for left Hexagon window
                    strFile = fullfile(strSavePath, 'Left.mat');
                    objL = matfile(strFile,'Writable',true);

                    % Make mat file for right Hexagon window
                    strFile = fullfile(strSavePath, 'Right.mat');
                    objR = matfile(strFile,'Writable',true);

                    % Save camera matrices in new mat files
                    objL.IntrinsicMatrix = objML.IntrinsicMatrix;
                    objR.IntrinsicMatrix = objML.IntrinsicMatrix;
                    objL.PoseMatrix = objML.LeftPoseMatrix;
                    objR.PoseMatrix = objML.RightPoseMatrix;

                    % Save source image info in mat files
                    objL.SourceImageInfo = struct('Path',objML.Properties.Source, ...
                        'SpatialTrans',objML.SpatialTrans);
                    objL.SourceImage = objML.Image10;
                    objR.SourceImageInfo = struct('Path',objMR.Properties.Source, ...
                        'SpatialTrans',objMR.SpatialTrans);
                    objR.SourceImage = objMR.Image10;

                    % Save window and ROI info in mat files
                    objL.Window = cWindow{iW}.left;
                    objR.Window = cWindow{iW}.right;
                    objL.WindowID = iW;
                    objR.WindowID = iW;
                    objL.RegionID = cWindow{iW}.region;
                    objR.RegionID = cWindow{iW}.region;

                    % Initialize
%                     cWin = {num2str(iW) num2str(iNumWin)};
                    objL.Accuracy = struct();

                    % Read Hexagon images
                    mWin = objL.Window;
                    objL.Image = objML.Image(mWin(3):mWin(4),...
                        mWin(1):mWin(2));
                    mWin = objR.Window;
                    objR.Image = objMR.Image(mWin(3):mWin(4),...
                        mWin(1):mWin(2));

                    % Rectify the stereo images
                    extStereoRect(objL,objR,strSavePath);

                    % Compute disparity map
                    extDisparity(objL,objR,strRes,iBlkSz);

                catch objExc

                    warning(objExc.message)
                    warning(['An error occurred during stereo rectification or ' ...
                        'disparity map computation. Skipping window ' num2str(iW) '...'])
                    objL.Error = objExc.message;
                    objR.Error = objExc.message;

                end
            end


        end

        function out = BundleAdjustLoop()


        end

        function out = TriangulateLoop()


        end


    end
end