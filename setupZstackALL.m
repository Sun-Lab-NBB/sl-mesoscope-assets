function [zstack] = setupZstackALL(hSI, hSICtl, zum, zrange, zmirror, order, channel)

% Limited argument validation and default value assignment support. May not
% work on older MatLab versions, but good for R2022b+.
arguments
    hSI  % ScanImage handle object. Cannot be validated due to how MBF implemented the class.
    hSICtl  % ScanImage Controller. Cannot be validated due to how MBF implemented the class.
    zum (1,1) double {mustBePositive, mustBeFinite} = 20  % plane z-spacing in micrometers
    zrange (1,2) double {mustBeFinite} = [1050,1050]  % [min, max] z-range. For single-lane imaging set to the same
    zmirror double {mustBeFinite} = []  % optional mirroring positions
    order {mustBeMember(order, {'', 'smooth'})} = ''  % acquisition order
    channel (1,1) double {mustBePositive, mustBeInteger} = 1  % imaging channel
end

input('checklist: a) Scan phase ~0.8888. b) STIMULUS, c) POWER = 70%');

% Example call:
% [u, zstack] = setupZstackALL(hSI, hSICtl, 25, [300 700], [], '');


% % Parameter Definition
global Z

averageNVolumes = 20; % how many volumes to take

nzhalf = min(floor((zum-1)/2),12);

if isempty(zmirror)
    centerZs = zrange(1):zum:zrange(2);
else
    centerZs = [[zmirror(1)-nzhalf:-zum:zrange(1)] [zmirror(2)+nzhalf:zum:zrange(2)]];
end

centerZs = sort(centerZs);
refZs = centerZs(:) + [-nzhalf:nzhalf];
if strcmp(order, 'smooth')
    refZs = refZs';
end

if hSI.hFastZ.enableFieldCurveCorr==1
%     warning('field curvature correction is ON, turning off')
%     hSI.hFastZ.enableFieldCurveCorr = 0;
end

%% Set up Reference ROIs
hSI.abort();

hSI.hFastZ.hFastZs{1}.move(min(centerZs))

% Grab Reference Volumes
fprintf('Grabbing reference volume\n');

hSI.hStackManager.arbitraryZs = sort(refZs(:));
hSI.hStackManager.numVolumes = averageNVolumes;
hSI.hStackManager.enable = true;

hSI.hDisplay.displayRollingAverageFactor = averageNVolumes; % buffer N Volumes in frame averaging buffer

hSI.hMotionManager.enable = false;
hSI.hRoiManager.mroiEnable = true;
hSI.hStackManager.stackReturnHome = true;

hSI.startGrab();
while hSI.active
    pause(1); % wait for reference volume to be completed
end
hSI.abort();

fprintf('Grabbing reference volume completed.\n');

%% Set up Motion Estimators
fprintf('Setting up Motion Estimators.\n');

hSI.hMotionManager.clearAndDeleteEstimators();
hSI.hMotionManager.estimatorClassName = 'scanimage.components.motionEstimators.MariusMotionEstimator';
%%
hSI.hMotionManager.correctorClassName = 'scanimage.components.motionCorrectors.MariusMotionCorrector2';

roiDatas = hSI.hDisplay.getRoiDataArray();
zstack = copy(roiDatas);
arrayfun(@(rd)rd.onlyKeepChannels(channel),roiDatas);

% first dimension is roi index, second dimension is volume index
nRois = size(roiDatas,1);
nVolumes = size(roiDatas,2);

fprintf('Aligning %d stacks.\n',nVolumes);


Z = [];

for roiIdx = 1:nRois
    roi0 = copy(roiDatas(roiIdx,:));
    for j = 1:averageNVolumes
        roi0(j).imageData = [];
    end
    
    for iz = 1:numel(centerZs)
        id = find(ismember(refZs', centerZs(iz) + [-nzhalf:nzhalf]));
        
        roi1 = copy(roi0);
        for j = 1:averageNVolumes
            roi1(j).imageData{1}  = roiDatas(roiIdx, j).imageData{1}(id);
            roi1(j).zs  = roi1(j).zs(id);
        end
        alignedRoiData = hSI.hMotionManager.alignZStack(roi1);
        img = alignedRoiData.imageData{1};
        for j = 1:numel(img)
            Z{roiIdx, iz}(:,:,j) = img{j};
        end
        hSI.hMotionManager.addEstimator(alignedRoiData);
    end
    
end


hSI.hMotionManager.enable = true;
hSICtl.showGUI('MotionDisplay');
%%

% Set up acquisition, but don't start it yet
fprintf('Preparing system for acquisition.\n');

hSI.hStackManager.stackDefinition = 'arbitrary';
hSI.hStackManager.stackMode = 'fast';
hSI.hStackManager.stackFastWaveformType = 'step';
hSI.hStackManager.arbitraryZs = centerZs(:);
disp(hSI.hStackManager.arbitraryZs)

hSI.hStackManager.numVolumes = 200000;
hSI.hStackManager.enable = true;

hSI.hDisplay.displayRollingAverageFactor = 1; % buffer N Volumes in frame averaging buffer

hSI.hMotionManager.enable = true;

hSI.hStackManager.numVolumes = 100000;

tau = 100; % was 50
fs = hSI.hRoiManager.scanFrameRate;
hSI.hMotionManager.hMotionCorrector.pC  = exp(-1/(tau*fs));
hSI.hMotionManager.correctionEnableZ    = 1;
hSI.hMotionManager.correctionBoundsZ    = [-100 100];
hSI.hMotionManager.hMotionCorrector.thresholdExceedTime_s = 10;
hSI.hMotionManager.hMotionCorrector.correctionInterval_s = 10;
hSI.hMotionManager.hMotionCorrector.correctionThreshold = [.1 .1 .5];


