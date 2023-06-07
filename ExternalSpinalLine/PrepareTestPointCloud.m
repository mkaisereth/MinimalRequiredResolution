function pc = PrepareTestPointCloud(eslParams, filePath, printLevel, showPcRoi)
if isempty(filePath)
    filePath = 'Photoneo_StaticBack_ClosedWindowShutters_20210623_100cm_fg.pcd';
end
if isstring(filePath) || ischar(filePath)
    pc = pcread(filePath);
end
if isa(filePath, 'pointCloud')
   pc = filePath; 
end

% convert to meter
if pc.YLimits(2)-pc.YLimits(1) > 10
   pcLoc = pc.Location;
   pcLoc = pcLoc/1000;
   pc = pointCloud(pcLoc);
   eslParams.DidConvertFromMToMM = 1;
end

if isempty(eslParams.ManualPcRoi) || eslParams.ManualPcRoi == 1
    addpath('..\Utils');
    if(~exist('pcRoi'))
        msgbox('You need to clone utils: https://gitlab.ti.bfh.ch/4d-spine/utils on relative path ..\Utils');
    end
    if exist('showPcRoi') && showPcRoi == 1
        pc = pcRoi(pc);
    end
end

% get ROI, if specified
if ~isempty(eslParams.ImgForegroundMin)
    ImgForegroundMin = eslParams.ImgForegroundMin;
    pcLoc = pc.Location;
    if size(ImgForegroundMin,1)>2 || size(ImgForegroundMin,2)>2
       pcLoc(pcLoc(:,1)<ImgForegroundMin(1),:) = [];
       pcLoc(pcLoc(:,2)<ImgForegroundMin(2),:) = [];
       pcLoc(pcLoc(:,3)<ImgForegroundMin(3),:) = [];
    else
        pcLoc(pcLoc(:,3)<ImgForegroundMin,:) = [];
    end
    pc = pointCloud(pcLoc);
end
if ~isempty(eslParams.ImgForegroundMax)
    ImgForegroundMax = eslParams.ImgForegroundMax;
    pcLoc = pc.Location;
    if size(ImgForegroundMax,1)>2 || size(ImgForegroundMax,2)>2
       pcLoc(pcLoc(:,1)>ImgForegroundMax(1),:) = [];
       pcLoc(pcLoc(:,2)>ImgForegroundMax(2),:) = [];
       pcLoc(pcLoc(:,3)>ImgForegroundMax(3),:) = [];
    else
        pcLoc(pcLoc(:,3)>ImgForegroundMax,:) = [];
    end
    pc = pointCloud(pcLoc);
end

global printFilter, global subjectNr;
if printLevel > 0 && (isempty(printFilter) || printFilter=="" || printFilter==subjectNr)
    f1 = figure('Name', 'OrigPc');
    pcshow(pc);
    title('original pointcloud')
    [pcFileFolderPath, pcFileName, ~] = fileparts(eslParams.FilePath);
    if ~exist(append(pcFileFolderPath, "/Output/"))
        mkdir(append(pcFileFolderPath, "/Output/"));
    end
    saveas(f1, append(pcFileFolderPath, "/Output/", pcFileName, "_origpc.fig"))
end

% prepare pointcloud
ds = eslParams.DownSample;
if eslParams.AbsoluteSampling > 0
    if pc.Count > 0
       tempDs = 1/(pc.Count/eslParams.AbsoluteSampling); % sample down to 10'000 pts
       if tempDs < ds
           ds = tempDs;
       end
       if ds > 1
           ds = 1;
       end
    end
end
if ds < 1
    pc = pcdownsample(pc,'random', ds);
end

% limit z resolution (if applicable)
if eslParams.LimitZResolution > 0
    pcLoc = pc.Location;
    pcLoc(:,3) = eslParams.LimitZResolution*round(pcLoc(:,3)/eslParams.LimitZResolution);
    pc = pointCloud(pcLoc);
end

% limit z frequency (if applicaple)
if eslParams.FrequencyFilter > 0
    pcLoc = pc.Location;
    
    % get horizontal slices
    [slices, xDistAvg, yDistAvg] = SliceHorizontally(pc, pcLoc, printLevel);
    
    filteredPc = [];
    for sliceInd = 1:size(slices,2)
        
        yVals = slices{sliceInd};
        yzVals = yVals(:,3);
        
        Ts = xDistAvg; % sampling period
        Fs = 1/Ts; % sampling frequency
        fc = eslParams.FrequencyFilter;

        %normalized cutoff frequency
        wn = fc/(Fs/2);
        [b,a] = butter(6,wn);

        if printLevel > 2
            figure;
            freqz(b,a)
        end
        if size(yzVals,1) <= 18
            disp("not enough data for filtering: " + sliceInd);
        else
            yzVals = filtfilt(b,a, double(yzVals));
        end
        filteredPc = [filteredPc;[yVals(:,1),yVals(:,2),yzVals(:,1)]];            
    end
    
    if printLevel > 1
        figure;
        pcshow(filteredPc);
    end
    pc = pointCloud(filteredPc);
end

% center the point cloud
centerX = (pc.XLimits(2)+pc.XLimits(1))/2;
centerY = (pc.YLimits(2)+pc.YLimits(1))/2;
centerZ = (pc.ZLimits(2)+pc.ZLimits(1))/2;

% make sure point cloud is in positive z, otherwise algorithm does not work
if pc.ZLimits(1) < 0.1
    pcLoc = pc.Location;
    pcLoc(:,3) = pcLoc(:,3)-pc.ZLimits(1)+0.1;
    pc = pointCloud(pcLoc);
    eslParams.DidShift0_1 = 1;
end

% rotate point cloud
theta = 0*pi/180;
rot = [cos(theta) 0 sin(theta); ...
               0          1  0;...
      -sin(theta) 0 cos(theta)];
trans = [0, 0, 0];
tform = rigid3d(rot,trans);
pc = pctransform(pc,tform);

% TODO
pc = pointCloud(pc.Location);

if printLevel > 1
    figure;
    pcshow(pc);
    title("original point cloud (centered, downsampled)")
end
end