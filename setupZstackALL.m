function [zstack] = setupZstackALL(hSI, hSICtl, zum, zrange, zmirror, order, channel)

% Limited argument validation and default value assignment support. May not
% work on older MatLab versions, but good for R2022b+.
arguments
    hSI  % ScanImage handle object. Cannot be validated due to how MBF implemented the class.
    hSICtl  % ScanImage Controller. Cannot be validated due to how MBF implemented the class.
    zum (1,1) double {mustBePositive, mustBeFinite, mustBeInteger} = 20  % Plane z-spacing in micrometers
    zrange double {mustBePositive, mustBeInteger, validateZRange} = 1050  % Single value or [min, max]. The range of imaged z-planes.
    zmirror double {mustBePositive, mustBeInteger} = []  % Optional mirroring positions
    order {mustBeMember(order, {'', 'smooth'})} = ''  % PLane acquisition order.
    channel (1,1) double {mustBePositive, mustBeInteger} = 1  % The channel to use for motion detection.
    enableFieldCurveCorr (1,1) logical = true  % Determines whether to use curvature correction.
end

% Converts single zrange values to [min, max] format expected by the rest of the function.
if isscalar(zrange)
    fprintf('Single plane imaging at z = %d µm.\n', zrange);
    zrange = [zrange, zrange];
else
    fprintf('Z-stack imaging from %d to %d µm.\n', zrange(1), zrange(2));
end

% Validates zmirror has correct format if provided
if ~isempty(zmirror)
    if numel(zmirror) ~= 2
        error('zmirror must be empty or contain exactly 2 values.');
    end
    zmirror = zmirror(:)';  % Ensures row vector
end

% Instructs the user to verify important imaging parameters before generating the reference stack.
input('Check that a) Scan phase ~0.8888. b) Laser power = 70% before continuing.');

%% Parameter Definition
global Z  % Not sure why this is global, but keeping it this way for now

% Determines how many frames (volumes) to average at each z-position.
averageNVolumes = 20;

% Calculates reference half-width, maxing out at 12 reference planes on either end of the 
% imaged plane range. This determines how many zum-spaced planes are acquired above and below
% the each of the supported imaging plane(s) to support Z-drift correction.
nzhalf = min(floor((zum-1)/2),12);

% Generates z-plane imaging positions.
if isempty(zmirror)
    % If zmirror is not provided, creates a set of 'zum'-spaced planes from minimum to maximum imaging plane. 
    centerZs = zrange(1):zum:zrange(2);
else
    % If zmirror is provided, creates two plane sequences expanding outward from each of the imaging 
    % focal points
    centerZs = [[zmirror(1)-nzhalf:-zum:zrange(1)] [zmirror(2)+nzhalf:zum:zrange(2)]];
end

% Sorts the center-points of each target plane resolved above and generates a set of 
% planes above and below each imaging plane (reference Z-planes).
centerZs = sort(centerZs);
refZs = centerZs(:) + [-nzhalf:nzhalf];

% 'Smooth' acquisition order acquires all frames for the target plane before moving to the
% next one (Z1-Z1-Z1-Z2-Z2-Z2). Default acquisition order loops over planes (Z1-Z2-Z1-Z2-Z1-Z2)
if strcmp(order, 'smooth')
    refZs = refZs';
end

% Depending on the configuration, ensures that FieldCurvatureCorredction is either enabled or disabled.
% For Mesoscope, it should be enabled in most cases.
if hSI.hFastZ.enableFieldCurveCorr ~= enableFieldCurveCorr
    if enableFieldCurveCorr
        fprintf('Field curvature correction: Enabled.\n');
    else
        fprintf('Field curvature correction: Disabled.\n');
    end
    hSI.hFastZ.enableFieldCurveCorr = enableFieldCurveCorr;
end

%% Reference ROI setup
% If an acquisition is active, aborts it before changing system configuration.
hSI.abort();  

% Moves to the lowest plane to be imaged. Assumes that the fast-z is inverted, so smaller planes are 
% actually closest tot he surface of the brain.
hSI.hFastZ.hFastZs{1}.move(min(centerZs))

% Grabs the reference volumes
fprintf('Grabbing reference volume...\n');

% Configures the acquisition to operate on the set of reference Z-planes and acquire the requested number of 
% frames at each plane (20).
hSI.hStackManager.arbitraryZs = sort(refZs(:));
hSI.hStackManager.numVolumes = averageNVolumes;
hSI.hStackManager.enable = true;

% Buffers frames in frame averaging buffer.
hSI.hDisplay.displayRollingAverageFactor = averageNVolumes;

% Disables motion correction during z-stack acquisition and ensures MROI mode is active.
hSI.hMotionManager.enable = false;
hSI.hRoiManager.mroiEnable = true;
hSI.hStackManager.stackReturnHome = true;

% Ensures that the grabbed frames are saved as 'zstack_0000.tiff' file.
hSI.hChannels.loggingEnable = true;  % Enable data logging
hSI.hScan2D.logAverageFactor = 1;  % Save every frame (no averaging in saved data)
hSI.hScan2D.logFileStem = 'zstack';  % Base filename
hSI.hScan2D.logFileCounter = 0;  % Starting file counter

% Activates frame acquisition (starts grabbing)
hSI.startGrab();
while hSI.active
    pause(1); % waits for reference volume to be completed
end

fprintf('Reference volumes: Grabbed.\n');

%% Motion Estimators generation
fprintf('Setting up Motion Estimators...\n');

hSI.hMotionManager.clearAndDeleteEstimators();  % Removes existing estimators.

% Configures MotionEstimation plugin to use Marius code.
hSI.hMotionManager.estimatorClassName = 'scanimage.components.motionEstimators.MariusMotionEstimator';
hSI.hMotionManager.correctorClassName = 'scanimage.components.motionCorrectors.MariusMotionCorrector2';

% Loads the ROI stack data from the ROI manager
roiDatas = hSI.hDisplay.getRoiDataArray();
zstack = copy(roiDatas);

% Filters the ROI data to only contain the motion registration channel data.
arrayfun(@(rd)rd.onlyKeepChannels(channel),roiDatas);

% First dimension is roi index, second dimension is volume index
nRois = size(roiDatas,1);
nVolumes = size(roiDatas,2);

fprintf('Aligning %d stacks...\n',nVolumes);

%  Aligns z-stacks
Z = [];

% Loops over all ROIs
for roiIdx = 1:nRois

% Copies ROI data from the ROI manager
    roi0 = copy(roiDatas(roiIdx,:));
    
    % Precreates a template storage structure with correct dimensions, but no image data.
    for j = 1:averageNVolumes
        roi0(j).imageData = [];  % Clears image data, but keeps the metadata
    end
    
    % Finds reference planes for each target imaging position
    for iz = 1:numel(centerZs)
        id = find(ismember(refZs', centerZs(iz) + [-nzhalf:nzhalf]));
        
        % Extracts the data for relevant reference z-planes.
        roi1 = copy(roi0);
        for j = 1:averageNVolumes
            roi1(j).imageData{1}  = roiDatas(roiIdx, j).imageData{1}(id);
            roi1(j).zs  = roi1(j).zs(id);
        end
        
        % Aligns the requested number of frames (default is 20) for each plane to 
        % create a stable reference point and adds generated reference frame data to
        % the storage tensor.
        alignedRoiData = hSI.hMotionManager.alignZStack(roi1);
        img = alignedRoiData.imageData{1};
        for j = 1:numel(img)
            Z{roiIdx, iz}(:,:,j) = img{j};
        end
        
        % Generates and adds the motion estimator for the target ROI to the Motion Detection
        % manager
        hSI.hMotionManager.addEstimator(alignedRoiData);
    end
    
end

% Re-enables the Motion Detection plugin and shows it to user.
hSI.hMotionManager.enable = true;
hSICtl.showGUI('MotionDisplay');
%%

%% Prepares the system for acquisition
fprintf('Preparing system for acquisition...\n');

hSI.hStackManager.stackDefinition = 'arbitrary';  % Enables arbitrary stack traversal.
hSI.hStackManager.stackMode = 'fast';  % Enables fast-z (voice-coil)
% Uses step mode to optionally support acquiring multiple frames per ROI
hSI.hStackManager.stackFastWaveformType = 'step';  
hSI.hStackManager.arbitraryZs = centerZs(:);  % Configures z-stack manager to target the requested imaging planes.
disp(hSI.hStackManager.arbitraryZs)  % Ensures that the stack manager window is displayed

% Activates the stack manager and ensures the number of volumes is set to a very high value. Since each acquisition runs
% until it acquires the requested number of frames, this effectively makes the acquisition only stop in response to manual
% or external triggers.
hSI.hStackManager.enable = true;
hSI.hStackManager.numVolumes = 100000;

% Presets frame averaging to 5 to give better picture at runtime.
hSI.hDisplay.displayRollingAverageFactor = 5;

% Ensures that external trigger mode is enabled.
hSI.hScan2D.trigAcqTypeExternal = true;  % Enable external triggering

% Ensures that the motion estimator is enabled
hSI.hMotionManager.enable = true;

% Configures the motion estimation parameters
tau = 100; % Increased from the default value of 50
fs = hSI.hRoiManager.scanFrameRate;
hSI.hMotionManager.hMotionCorrector.pC  = exp(-1/(tau*fs));
hSI.hMotionManager.correctionEnableZ    = 1;  % Enables z-correction
hSI.hMotionManager.correctionBoundsZ    = [-100 100];  % Z-correction is performed within +- 100 um of target plane
hSI.hMotionManager.hMotionCorrector.thresholdExceedTime_s = 5;  % Decreased from the default value of 10
hSI.hMotionManager.hMotionCorrector.correctionInterval_s = 5;  % Decreased from the default value of 10
hSI.hMotionManager.hMotionCorrector.correctionThreshold = [.1 .1 .5];  # x, y, z


