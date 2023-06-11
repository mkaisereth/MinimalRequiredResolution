function [X,Y,Z,pc] = PrepareSmoothedTestPointCloud(eslParams, fileName, printLevel, showPcRoi)

if isempty(eslParams.DownSample)
    eslParams.DownSample = 0.4;
end
pc = PrepareTestPointCloud(eslParams, fileName, printLevel, showPcRoi);

% smooth the pointcloud
if isempty(eslParams.SmoothRadius)
    eslParams.SmoothRadius = 0.01;
end
disp('Going to do pcmedian ...');
% do some sanity checks
xDim = pc.XLimits(2)-pc.XLimits(1);
yDim = pc.YLimits(2)-pc.YLimits(1);
ptsPerXDim = pc.Count/xDim;
ptsPerXDimYDim = ptsPerXDim/yDim;
ptsPerRadius = ptsPerXDimYDim*eslParams.SmoothRadius;

limitPts = 10000; % it seems with 14000 points it takes a long time ...
if ptsPerRadius > limitPts
    pc = pcdownsample(pc, 'random', limitPts/ptsPerRadius);
    
    % check
    ptsPerXDim = pc.Count/xDim;
    ptsPerXDimYDim = ptsPerXDim/yDim;
    ptsPerRadius = ptsPerXDimYDim*eslParams.SmoothRadius;
    if abs(ptsPerRadius-limitPts) > 100
        disp('this did not work, something is wrong ...')
    end
end

pc = pcmedian(pc, 'Radius', eslParams.SmoothRadius);%0.01);
disp('Pcmedian done.');

if printLevel > 1
    figure;
    pcshow(pc);
    title("original point cloud (smoothed)")
end

if isempty(eslParams.DiscreteSteps)
    eslParams.DiscreteSteps = 120;
end
% approximate ordered point cloud
pcXMin = pc.XLimits(1);
pcXMax = pc.XLimits(2);
pcYMin = pc.YLimits(1);
pcYMax = pc.YLimits(2);
deltaX = (pcXMax-pcXMin)/eslParams.DiscreteSteps;
deltaY = (pcYMax-pcYMin)/eslParams.DiscreteSteps;
delta = min([deltaX,deltaY]);
pcXGrid = pcXMin:delta:pcXMax;
pcYGrid = pcYMin:delta:pcYMax;

[X,Y] = meshgrid(pcXGrid,pcYGrid);

Z = [];
delta2 = 5*delta;
for xi = 1:size(X,1)
    for yi = 1:size(X,2)
        inds = findPointsInROI(pc, [X(xi,yi)-delta2 X(xi,yi)+delta2 Y(xi,yi)-delta2 Y(xi,yi)+delta2 0 Inf]);
        roiPts = pc.Location(inds,:);
        zMedian = median(roiPts(:,3));
        Z(xi,yi) = zMedian;
    end
end

% smooth again TODO use more robust gradient
if isempty(eslParams.DiscreteKernel)
   eslParams.DiscreteKernel = 9; 
end
kernelSize=eslParams.DiscreteKernel;
yKernel = 33;
K = (1/(kernelSize*yKernel))*ones(yKernel, kernelSize);
Z = nanconv(Z,K,'edge', 'nanout');

if printLevel > 1
    figure;
    tempZ = Z;
    %tempZ(tempZ<1.0) = nan;
    surf(X,Y,reshape(tempZ,size(Z)));
    xlabel('x')
    ylabel('y')
    title('Smoothed point cloud rasterized')
end
end