% this code was developed by Marius Pachitariu, Carsen Stringer and Georg Jaindl

classdef MariusMotionCorrector2 < scanimage.interfaces.IMotionCorrector    
    properties (SetObservable)
        runningAverageLength_s = 3;     % running average length in seconds for motion correction value
        correctionInterval_s = 10;      % interval for motion correction events in seconds
        correctionThreshold = [1 1 10]; % units are in [degrees degrees micron]
        thresholdExceedTime_s = 0;      % only correct if dR has been above the threshold for N seconds
        pC = 0.975;                     % time smoothing factor for the z-correlation
        zUpSampling = 10;               % upsampling factor for z-correlation
        separateZs = false;
        xyPixelCorrelationWindow = 100;
    end
    
    properties (SetObservable, Dependent)
        showPlots;  % [logical] shows / hides plots
    end
    
    properties (Access = private)
        motionHistory;  % stores the motion history
        lastUpdate;     % stores the time the last motion update happened
        drRefHistory;
        drRefHistoryTimestamps;
        drRefNeedsUpdate = true;
        lastDrRef = [];
        hConfidencePlot;
        hCorrelationPlot;
        hDeltaZPlot;
        lastTimeStamp = -Inf;
        zs = [];
        dx = [];
        dy = [];
        dz = [];
        xCMean = {};
        yCMean = {};
        zCMean = {};     % running average of z correlation
        xCnm = [];       % normalization factor for z correlation
        yCnm = [];       % normalization factor for z correlation
        zCnm = [];       % normalization factor for z correlation
        dZlastCorrection = [];
    end
    
    %% Lifecycle
    methods
        function obj = MariusMotionCorrector2()
            obj.initPlots();
        end
        
        function delete(obj)
            obj.deinitPlots();
        end
    end
    
    %% Abstract Method Realization
    methods        
        function start(obj)
            obj.lastUpdate = tic();
            obj.drRefNeedsUpdate = true;
            obj.resetRefHistory();
            obj.hConfidencePlot.reset();
            obj.hCorrelationPlot.reset();
            obj.hDeltaZPlot.reset();
            obj.resetZCmean();
            obj.lastTimeStamp = -Inf;
            obj.dZlastCorrection = [];
        end
        
        function abort(obj)
            % No-op
        end
        
        function updateMotionHistory(obj,motionHistory)
            % Note: in future releases the motionHistory might not be
            % sorted anymore. Use [motionHistory.historyIdx] to get the
            % history indices
            obj.motionHistory = motionHistory;
            obj.drRefNeedsUpdate = true;
            
            drRef = obj.getCorrection(); % update correction            
            
            timeToUpdate = toc(obj.lastUpdate) > obj.correctionInterval_s;
            if ~timeToUpdate
                return
            end
            
            historyMask = obj.drRefHistoryTimestamps >= obj.drRefHistoryTimestamps(end)-obj.thresholdExceedTime_s;
            
            drRefs = obj.drRefHistory(historyMask,:);

            threshExceed = false(1,3);
            for axIdx = 1:3
                da = drRefs(:,axIdx);
                da = da(~isnan(da));
                threshExceed(axIdx) = all(abs(da) >= obj.correctionThreshold(axIdx));
            end
            aboveCorrectionThreshold = any(threshExceed);
            
            if aboveCorrectionThreshold && timeToUpdate
                obj.lastUpdate = tic();
                obj.notify('correctNow');
            end
        end
        
        function drRef = getCorrection(obj)
            if ~obj.drRefNeedsUpdate
                drRef = obj.lastDrRef;
                return
            end
            
            if isempty(obj.motionHistory)
                drRef = [NaN NaN NaN];
                return
            end
            
            histLength = numel(obj.motionHistory);
            
            timestamps = [obj.motionHistory.frameTimestamp];
            drRefs = vertcat(obj.motionHistory.drRef);
            
            timestamp = timestamps(end);
            newmask = timestamps > obj.lastTimeStamp;
            newtimestamps = unique(timestamps(newmask));
            
            for idx = 1:length(newtimestamps)
                t = newtimestamps(idx);
                mask = timestamps == t;
                [dx_,dy_,dz_] = obj.updateTimestep(t,mask);
            end
            
            % filter for runningAverageLength
            historymask = timestamps>=timestamp-obj.runningAverageLength_s;
            
            drRef = drRefs(historymask,:);
            drRef = mean(drRef,1,'omitnan');
            
            %drRef(3) = dz_;
            drRef = [dx_,dy_,dz_];
            
            obj.drRefHistory = append(obj.drRefHistory,drRef,histLength);
            obj.drRefHistoryTimestamps = append(obj.drRefHistoryTimestamps,timestamp,histLength);
            
            drRef(isnan(drRef)) = 0;
            
            obj.lastDrRef = drRef;            
            obj.drRefNeedsUpdate = false;
            obj.lastTimeStamp = timestamp;
        end
        
        function correctedMotion(obj,dr,motionCorrectionVector)
            dz_ = dr(3);
            if ~isnan(dz_) && dz_~=0
                obj.resetZCmean();
                obj.dZlastCorrection = dz_;
            end
        end
    end
    
    methods (Hidden) 
        function s = saveUserData(obj)
            % override method in IMotionCorrector to store figure positions
            s = struct();
            s.hConfidencePlot_Position = obj.hConfidencePlot.hFig.Position;
            s.hCorrelationPlot_Position = obj.hCorrelationPlot.hFig.Position;
            s.hDeltaZPlot_Position = obj.hDeltaZPlot.hFig.Position;
        end
        
        function loadUserData(obj,s)
            % override method in IMotionCorrector to restore figure positions
            obj.hConfidencePlot.hFig.Position = s.hConfidencePlot_Position;
            obj.hCorrelationPlot.hFig.Position = s.hCorrelationPlot_Position;
            obj.hDeltaZPlot.hFig.Position = s.hDeltaZPlot_Position;
        end
    end
    
    %% Internal methods
    methods (Access = private)        
        function resetRefHistory(obj)
            obj.drRefHistory = [];
            obj.drRefHistoryTimestamps = [];
        end
        
        function resetZCmean(obj)
            obj.zs = [];
            obj.dz = [];
            obj.xCMean = {};
            obj.yCMean = {};
            obj.zCMean = {};
            obj.xCnm = [];
            obj.yCnm = [];
            obj.zCnm = [];
            obj.hCorrelationPlot.reset();
        end
        
        function initPlots(obj)
            obj.deinitPlots();
            
            obj.hCorrelationPlot = most.util.TimePlot('Correlation Plot',false);
            obj.hCorrelationPlot.xLabel = 'Z';
            obj.hCorrelationPlot.yLabel = 'Correlation';
            obj.hCorrelationPlot.historyLength = 1;
            obj.hCorrelationPlot.changedCallback = @eigenSetPlots;
            obj.hCorrelationPlot.YLim = [0,NaN];
            obj.hCorrelationPlot.legend = {{'z-corr','z-corr-mean'},'Location','northwest'};
            obj.hCorrelationPlot.yLimAnimationEnabled = false;
            
            obj.hConfidencePlot = most.util.TimePlot('Confidence Time Plot',false);
            obj.hConfidencePlot.xLabel = 'Acquisition Time [s]';
            obj.hConfidencePlot.yLabel = 'Confidence Value';
            obj.hConfidencePlot.historyLength = 100;
            obj.hConfidencePlot.changedCallback = @eigenSetPlots;
            obj.hConfidencePlot.legend = {{'confidence','confidence-mean'},'Location','northwest'};
            obj.hConfidencePlot.yLimAnimationEnabled = false;
            
            obj.hDeltaZPlot = most.util.TimePlot('Delta Z Plot',false);
            obj.hDeltaZPlot.xLabel = 'Acquisition Time [s]';
            obj.hDeltaZPlot.yLabel = 'Delta Z [um]';
            obj.hDeltaZPlot.historyLength = 100;
            obj.hDeltaZPlot.changedCallback = @eigenSetPlots;
            obj.hDeltaZPlot.lineSpecs = {...
                {},{},... % z, z_mean
                {'Marker','x','LineStyle','none','LineWidth',2,'Color','red'},... % correction marker
                {'LineStyle','--','Color','black'},{'LineStyle','--','Color','black'}}; % threshold lines 
            obj.hDeltaZPlot.legend = {{'delta z','delta z mean','delta z correction'},'Location','northwest'};
            obj.hDeltaZPlot.yLimAnimationEnabled = false;
            
            most.gui.tetherGUIs(obj.hCorrelationPlot.hFig,obj.hConfidencePlot.hFig,'righttop');
            most.gui.tetherGUIs(obj.hConfidencePlot.hFig,obj.hDeltaZPlot.hFig,'righttop');
            
            function eigenSetPlots(src)
                % this is to update the user interface
                if obj.showPlots ~= src.visible
                    obj.showPlots = src.visible;
                end
            end
        end
        
        function deinitPlots(obj)
            most.idioms.safeDeleteObj(obj.hConfidencePlot);
            most.idioms.safeDeleteObj(obj.hCorrelationPlot);
            most.idioms.safeDeleteObj(obj.hDeltaZPlot);
        end
        
        function [dx,dy,dz] = updateTimestep(obj,timestamp,mask)     
            % calculate Z correlation based on sliding average of Z correlation
            im_zs = [obj.motionHistory(mask).z];            % z of imaged plane
            im_corrs = {obj.motionHistory(mask).correlation};  % {x,y,z} correlations
            ref_zs = {obj.motionHistory(mask).zs};           % estimator reference zs
            
            zCorrAvailableMask = cellfun(@(c)numel(c)>=3&&~isempty(c{3}),im_corrs); % only use estimates that actually return correlations for z
            if ~all(zCorrAvailableMask)
                most.idioms.warn('Motion Corrector: Not all estimators return a z correlation');
            end
            
            % filter by zCorrAvailableMask
            im_zs = im_zs(zCorrAvailableMask);
            im_corrs = im_corrs(zCorrAvailableMask);
            ref_zs = ref_zs(zCorrAvailableMask);
            
            if isempty(im_zs)
                dx = 0;
                dy = 0;
                dz = 0;
                return % nothing to do here
            end
            
            % Sanity check; in ScanImage, one timestep is one frame is one z
            assert(all(im_zs == im_zs(1)),'Motion Corrector: zs needs to be unique');
            z = im_zs(1);
            
            % get reference zs
            assert(all(cellfun(@(cz)isequal(cz(:),ref_zs{1}(:)),ref_zs)),'Motion Corrector: Cannot process z correlations for estimators with different zs');
            ref_zs = ref_zs{1}(:);

            % get correlations
            xCorrs = cellfun(@(c)c{1}(:),im_corrs,'UniformOutput',false);
            xCorrs = cellfun(@(c)cropCorrelation(c,obj.xyPixelCorrelationWindow),xCorrs,'UniformOutput',false);
            xCorrs = horzcat(xCorrs{:});
            xCorrs = mean(xCorrs,2,'omitnan'); % average x correlations of current timestep
            
            yCorrs = cellfun(@(c)c{2}(:),im_corrs,'UniformOutput',false);
            yCorrs = cellfun(@(c)cropCorrelation(c,obj.xyPixelCorrelationWindow),yCorrs,'UniformOutput',false);
            yCorrs = horzcat(yCorrs{:});
            yCorrs = mean(yCorrs,2,'omitnan'); % average y correlations of current timestep
            
            zCorrs = cellfun(@(c)c{3}(:),im_corrs,'UniformOutput',false);
            zCorrLengths = cellfun(@(zc)numel(zc),zCorrs);
            assert(all(zCorrLengths == zCorrLengths(1)),'Motion Corrector: Cannot process z correlations of varying lengths.');
            zCorrs = horzcat(zCorrs{:});
            zCorrs = mean(zCorrs,2,'omitnan'); % average z correlations of current timestep
            
            if obj.separateZs
                idx = find(obj.zs == z);
                if isempty(idx)
                    idx = numel(obj.zs)+1;
                end
            else
                idx = 1;
            end
            
            if numel(obj.zs) < idx
                obj.xCMean{end+1} = zeros(size(xCorrs),'like',xCorrs);
                obj.xCnm(end+1) = 0;
                
                obj.yCMean{end+1} = zeros(size(yCorrs),'like',yCorrs);
                obj.yCnm(end+1) = 0;
                
                obj.zs(end+1) = z;
                obj.dz(end+1) = NaN;
                obj.zCMean{end+1} = zeros(size(zCorrs),'like',zCorrs);
                obj.zCnm(end+1) = 0;
            end
            
            obj.xCMean{idx} = obj.pC .* obj.xCMean{idx} + (1-obj.pC) .* xCorrs(:); % add correlation to running average X profile
            obj.xCnm(idx) = obj.pC .* obj.xCnm(idx) + (1-obj.pC); % normalizer for pC
            
            obj.yCMean{idx} = obj.pC .* obj.yCMean{idx} + (1-obj.pC) .* yCorrs(:); % add correlation to running average Y profile
            obj.yCnm(idx) = obj.pC .* obj.yCnm(idx) + (1-obj.pC); % normalizer for pC
            
            obj.zCMean{idx} = obj.pC .* obj.zCMean{idx} + (1-obj.pC) .* zCorrs(:); % add correlation to running average Z profile
            obj.zCnm(idx) = obj.pC .* obj.zCnm(idx) + (1-obj.pC); % normalizer for pC
            
            [xx_up,x_cc_up,x_cc_mean_up] = upsample(obj.xCMean{idx},obj.xCnm(idx),(1:numel(xCorrs)) - (numel(xCorrs)+1)/2,xCorrs);
            [yy_up,y_cc_up,y_cc_mean_up] = upsample(obj.yCMean{idx},obj.yCnm(idx),(1:numel(yCorrs)) - (numel(yCorrs)+1)/2,yCorrs);
            [zz_up,z_cc_up,z_cc_mean_up] = upsample(obj.zCMean{idx},obj.zCnm(idx),ref_zs,zCorrs);
            
            sf = obj.motionHistory(end).hRoi.scanfields(1);
            pixelSizeXY = sf.sizeXY ./ sf.pixelResolutionXY;
            
            % upsampled instantaneous correlation
            [cx,cx_idx] = max(x_cc_up);
            dx_instantenous = xx_up(cx_idx);
            dx_instantenous = dx_instantenous * pixelSizeXY(1);
            
            [cy,cy_idx] = max(y_cc_up);
            dy_instantenous = yy_up(cy_idx);
            dy_instantenous = dy_instantenous * pixelSizeXY(2);
            
            [c,c_idx] = max(z_cc_up);
            z_new = zz_up(c_idx);
            dz_instantenous = z - z_new;
            
            % upsampled time averaged correlation
            [cx_mean_up,cx_mean_up_idx] = max(x_cc_mean_up);
            dx_mean = xx_up(cx_mean_up_idx);
            dx_mean = dx_mean * pixelSizeXY(1);
            
            [cy_mean_up,cy_mean_up_idx] = max(y_cc_mean_up);
            dy_mean = yy_up(cy_mean_up_idx);
            dy_mean = dy_mean * pixelSizeXY(2);
            
            [c_mean_up,c_mean_up_idx] = max(z_cc_mean_up);
            z_mean_new = zz_up(c_mean_up_idx);
            dz_mean = z - z_mean_new;
            
            % return values
            obj.dx(idx) = dx_mean;
            obj.dy(idx) = dy_mean;
            obj.dz(idx) = dz_mean;

            dx = dx_mean;
            dy = dy_mean;
            dz = dz_mean;
            
            if obj.showPlots
                obj.hCorrelationPlot.addTimePoint([z_cc_up,z_cc_mean_up]);
                obj.hConfidencePlot.addTimePoint([c,c_mean_up],timestamp);
                correctionPoint = NaN;
                if ~isempty(obj.dZlastCorrection)
                    correctionPoint = obj.dZlastCorrection;
                    obj.dZlastCorrection = [];
                end
                obj.hDeltaZPlot.addTimePoint([dz_instantenous,dz_mean,correctionPoint,[-1,1]*obj.correctionThreshold(3)],timestamp);
            end
            
            function [pts_up,cc_up,cc_mean_up] = upsample(CMean,Cnm,pts,Corrs)
                % upsampling of correlation
                if numel(pts) == 1
                    pts_up = pts;
                    cc_up = Corrs(:);
                    cc_mean_up = CMean/Cnm;
                else
                    ii = (1:numel(pts))';
                    ii_up = linspace(1,ii(end),(numel(ii)-1)*obj.zUpSampling+1)';
                    pts = pts(:);
                    F_pts = griddedInterpolant(ii,pts,'linear');
                    pts_up = F_pts(ii_up);
                    F_cc = griddedInterpolant(pts,Corrs(:),'spline');
                    cc_up = F_cc(pts_up);
                    cc_mean = CMean/Cnm;
                    F_cc_mean = griddedInterpolant(pts,cc_mean,'spline');
                    cc_mean_up = F_cc_mean(pts_up);
                end
            end
            
            function corr = cropCorrelation(corr,windowSize)
                mask = (1:windowSize) - (windowSize+1)/2;
                mask = ceil(mask + (numel(corr)-1)/2);
                
                corr = corr(mask);
            end
        end
    end
    
    %% Property Getter/Setter
    methods
        function set.runningAverageLength_s(obj,val)
            validateattributes(val,{'numeric'},{'scalar','positive','nonnan','real','finite'});
            obj.runningAverageLength_s = val;
        end
        
        function set.correctionInterval_s(obj,val)
            validateattributes(val,{'numeric'},{'scalar','positive','nonnan','real','finite'});
            obj.correctionInterval_s = val;
        end
        
        function set.correctionThreshold(obj,val)
            validateattributes(val,{'numeric'},{'vector','row','numel',3,'nonnegative','nonnan','real'});
            obj.correctionThreshold = val;
        end
        
        function val = get.showPlots(obj)
            val = obj.hConfidencePlot.visible;
            val = val || obj.hCorrelationPlot.visible;
            val = val || obj.hDeltaZPlot.visible;
        end
        
        function set.showPlots(obj,val)
            validateattributes(val,{'numeric','logical'},{'scalar','binary'});
            obj.hConfidencePlot.visible = val;
            obj.hCorrelationPlot.visible = val;
            obj.hDeltaZPlot.visible = val;
        end
        
        function set.pC(obj,val)
            validateattributes(val,{'numeric'},{'scalar','nonnegative','real','finite','<=',1});
            obj.pC = val;
            obj.resetZCmean();
        end
        
        function set.zUpSampling(obj,val)
            validateattributes(val,{'numeric'},{'scalar','integer','positive','real'});
            obj.zUpSampling = val;
            obj.resetZCmean();
            obj.hCorrelationPlot.reset();
        end
        
        function set.separateZs(obj,val)
            validateattributes(val,{'numeric','logical'},{'scalar','binary'});
            obj.separateZs = val;
            obj.resetZCmean();
        end
    end
end

function vec = append(vec,v,veclength)
if isempty(vec)
    vec = v;
elseif size(vec,1) < veclength
    vec(end+1,:) = v;
else
    if size(vec,1) > veclength
        vec = vec(end-veclength+1:end,:);
    end
    vec = circshift(vec,-1,1);
    vec(end,:) = v;
end
end
