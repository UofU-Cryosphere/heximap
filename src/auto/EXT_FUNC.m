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

        function out = SortImages()
        
        
        end

        function out = InitTrans()


        end

        function out = DefineRegion()


        end

        function out = DisparityLoop()


        end

        function out = BundleAdjustLoop()


        end

        function out = TriangulateLoop()


        end


    end
end