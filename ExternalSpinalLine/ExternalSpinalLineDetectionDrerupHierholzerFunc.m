function [spLine, pLine, lStability] = ExternalSpinalLineDetectionDrerupHierholzerFunc(eslParams, fileName, printLevel, showPcRoi, evalStability)
%%
% According to paper: Back shape measurement using video rasterstereography
% and three-dimensional reconstruction of spinal shape (Drerup/Hierholzer
% 1994)
%
% 1. vertebral rotation is correlated with surface rotation
%
% 2. surface rotation is measured by the angle of the surface normal
% (horizontal component at the line of the spinous processes), i.e. by local data only
%
% According to Figure 2 the construction of the spinal midline requires three types of input data: 
% 1. the line of the spinous processes; 
% 2. surface rotation at locus of the line of the spinous processes; 
% 3. anatomical landmarks (needed for reference to the underlying skeletal structures). 
%
% Line of the spinous processes:
% - estimated by symmetry line of the back (horizontal profiles)
% - symmetry point is defined which divides the profile into two halves
% with minimum lateral asymmetry (w.r.t surface curvature): details:
% Analysis Of Left-Right Asymmetry Of The Back Shape Of Scoliotic Patients
% (Hierholzer 1986)
% - both principal curvatures K1, K2
% (https://faculty.sites.iastate.edu/jia/files/inline-files/surface-curvature.pdf)
%
%% Asymmetry according to Analysis Of Left-Right Asymmetry Of The Back Shape Of Scoliotic Patients (Hierholzer, 1986)
% a =  (H_L - H_R)^2 + (G_L^2 - 2*G_L*G_R*cos(2*e) + G_R^2)/2
% where
% H = (K1 + K2)/2 (mean curvature)
% G = (K1 - K2)/2 (difference of principal curvatures)
% K1 and K2 principal curvatures at Point P_L resp. P_R
% e = phi_L - phi_R (difference angle between corresponding principal
% directions; if the surface normals are not assumed to be parallel, the
% mean value for both principal direction should be taken)

[X,Y,Z,pc] = PrepareSmoothedTestPointCloud(eslParams, fileName, printLevel, showPcRoi);

% TODO
Z(isnan(Z(:))) = 0;

[gMap,dV,dL] = GaussianMap(X,Y,Z);

% assume x = 0 as mirror plane
sizeX = size(X,2);
sizeY = size(X,1);

asyms_mean = nan(sizeY,sizeX);

borderMarginX = 5;
borderMarginX2 = borderMarginX+5;
for mirror_xi = borderMarginX2+1:sizeX-borderMarginX2

    for yi = 5:sizeY-4
        a_temps = nan(sizeX,1);
       for xi = borderMarginX+1:sizeX-borderMarginX
           % only do left size up to mirror plane
           if xi >= mirror_xi
              break; 
           end
           % if nan, continue
           if Z(yi,xi) == 0
              continue; 
           end
           % find right mirror point
           xi_R = 2*mirror_xi-xi;
           if xi_R > sizeX-borderMarginX || Z(yi,xi_R) == 0
               continue;
           end
%            for xii = xi+1:sizeX-borderMarginX
%                if X(yi,xii)-X(yi,mirror_xi)>=X(yi,mirror_xi)-X(yi,xi)
%                    xi_R = xii;
%                    break;
%                end
%            end
%            if isnan(xi_R) % find counterpart, if none, stop
%                continue;
%            end
           if abs((xi_R-mirror_xi)-(mirror_xi-xi))>1 % sanity check
               disp('hö?');
           end

           V1_L = dV{yi,xi}(:,1);
           V2_L = dV{yi,xi}(:,2);
           L1_L = dL{yi,xi}(1);
           L2_L = dL{yi,xi}(2);
           if abs(L1_L) > abs(L2_L)
               K1_L = L1_L;
               K2_L = L2_L;
           else
               K1_L = L2_L;
               K2_L = L1_L;
           end
           if abs(K1_L) > 1000
               continue;
           end
           H_L = (K1_L+K2_L)/2;
           G_L = (K1_L-K2_L)/2;
           u = V1_L;
           v = [1;0]; % x-axis
           CosTheta = dot(u,v)/(norm(u)*norm(v));
           phi_L = acos(CosTheta);

           V1_R = dV{yi,xi_R}(:,1);
           V2_R = dV{yi,xi_R}(:,2);
           L1_R = dL{yi,xi_R}(1);
           L2_R = dL{yi,xi_R}(2);
           if abs(L1_R) > abs(L2_R)
               K1_R = L1_R;
               K2_R = L2_R;
           else
               K1_R = L2_R;
               K2_R = L1_R;
           end
           if abs(K1_R) > 1000
               continue;
           end
           
           H_R = (K1_R+K2_R)/2;
           G_R = (K1_R-K2_R)/2;
           u = V1_R;
           v = [1;0]; % x-axis
           CosTheta = dot(u,v)/(norm(u)*norm(v));
           phi_R = acos(CosTheta);

           e = phi_L-phi_R;

           a =  (H_L - H_R)^2 + (G_L^2 - 2*G_L*G_R*cos(2*e) + G_R^2)/2;
           if a > 1000
               %disp('large');
               continue;
           end
           a_temps(xi) = a;
           if ~isreal(a)
               disp('complex');
           end
       end
       asyms_mean(yi,mirror_xi) = mean(a_temps, 'omitnan');
       if ~isreal(asyms_mean(yi,mirror_xi))
           disp('complex');
           asyms_mean(yi,mirror_xi) = nan; % TOOD think about this
       end
    end

end

% make plot a bit nicer % TODO fix proper later
q_max = quantile(asyms_mean(:), 0.95); % cut the highest 5%
temp_asyms = asyms_mean(:);
temp_asyms(temp_asyms>q_max) = q_max;

q_min = quantile(asyms_mean(:), 0.05); % cut the lowest 5%
temp_asyms(temp_asyms<q_min) = q_min;

temp_asyms = reshape(temp_asyms, size(Z));
tempZ = Z;
% remove 0s (no data)
tempZ(tempZ<0.1) = nan;
tempZ = reshape(tempZ, size(Z));

global printFilter;
global subjectNr;
if printLevel > 0 && (isempty(printFilter) || printFilter=="" || printFilter==subjectNr)

    f2 = figure('Name', 'AsymmetryMap');
    surf(X,Y,tempZ, temp_asyms);
    xlabel('x')
    ylabel('y')
    title('asymmetrie map (minima are possible symmetry line candidates)(Drerup/Hierholzer, 1986)');
    colorbar;
    set(gca,'ColorScale','log')
    axis equal;
    
    % TODO
    [pcFileFolderPath, pcFileName, ~] = fileparts(eslParams.FilePath);
    saveas(f2, append(pcFileFolderPath, "/Output/", pcFileName, "_asymmMap.fig"))

end

% get the external spinal line (minimum asymmetry)
[amin, minInd] = min(asyms_mean,[],2);
qmax = quantile(asyms_mean,0.3,2);
minInds=sub2ind(size(asyms_mean),(1:size(asyms_mean,1))',minInd);
qmaxZ = tempZ;
qmaxZ(asyms_mean>=qmax) = nan;

% Not required anymore, use from above
%pc = PrepareTestPointCloud(downSample, fileName, printLevel);

%symmPts = [X(minInds),Y(minInds),tempZ(minInds)];
symmPts = [X(minInds),Y(minInds),qmaxZ];
if printLevel > 1
    figure;
    
    pcshow(pc);
    hold on;
    %pcshow(symmPts, [1 0 0]);
    %scatter3(X(minInds),Y(minInds),tempZ(minInds), 100, 'red', '.');
    scatter3(X(:),Y(:),qmaxZ(:), 100, 'red', '.');
    title("point cloud with quantile minima (Drerup/Hierholzer, 1986)");
end
% find the most centered local minima (TODO solve this with anatomical
% landmark conditions)
xMin = pc.XLimits(1);
xMax = pc.XLimits(2);
xMiddle = (xMin+xMax)/2;
lminInds = nan(size(asyms_mean,1),1);
for y=1:size(asyms_mean,1)
    lminVal = inf;
    lminInd = nan;
    for x=1:size(asyms_mean,2)
        if(~isnan(qmaxZ(y,x)) && abs(X(y,x)-xMiddle) < lminVal)
            lminVal = abs(X(y,x)-xMiddle);
            lminInd = x;
        end
    end
    lminInds(y) = lminInd;
end

% this is the closest min to the middle of the point cloud
qmaxZ2 = qmaxZ;
tempLminInds = lminInds;
tempLminInds(isnan(tempLminInds)) = 1;
minInds=sub2ind(size(asyms_mean),(1:size(asyms_mean,1))',tempLminInds);
if printLevel > 1
    figure;
    pcshow(pc);
    hold on;
    
    scatter3(X(minInds),Y(minInds),qmaxZ2(minInds)+0.001, 100, 'red', '.');
    title("point cloud with closest min to middle");
end

% now find the real min in that cluster
aminInds = nan(size(asyms_mean,1),1);
for y=1:size(asyms_mean,1)
    aminVal = inf;
    aminInd = nan;
    if ~isnan(lminInds(y))
        for x=lminInds(y):1:size(asyms_mean,2)
            if isnan(qmaxZ(y,x))
                break;
            end
            if asyms_mean(y,x) < aminVal
                aminVal = asyms_mean(y,x);
                aminInd = x;
            end
        end
        for x=lminInds(y):-1:1
            if isnan(qmaxZ(y,x))
                break;
            end
            if asyms_mean(y,x) < aminVal
                aminVal = asyms_mean(y,x);
                aminInd = x;
            end
        end
    end
    aminInds(y) = aminInd;
end

if evalStability == 1
    % evaluate stability of external spinal line
    % get 2nd derivate around external spinal line in x-direction
    minStabilityDeltaValue = 0;
    minStabilityMaxValue = -Inf;
    for deltaValue = 3:2:size(temp_asyms,2)/6 % allow up to 1/3 of overall width
        minStability = 0;
        for yi=1:size(aminInds,1)
            if isnan(aminInds(yi,:)) || aminInds(yi,:)-deltaValue <1 || aminInds(yi,:)+deltaValue > size(temp_asyms,2)
                % if no minima is found, this coun't as not stable
                continue;
            end
            diff1 = diff(temp_asyms(yi,aminInds(yi,:)-deltaValue:aminInds(yi,:)+deltaValue));
            diff2 = diff(diff1);
            % - because curvature is negative for good minima (z direction is to the front)
            currMinStability = -(mean(diff2(ceil(size(diff2,2)/2+1):size(diff2,2))) + mean(diff2(1:floor(size(diff2,2)/2))));
            if isnan(currMinStability)
                continue;
            end
            minStability = minStability + currMinStability;
        end
        minStability = minStability/size(aminInds,1);
        if minStability > minStabilityMaxValue
            minStabilityMaxValue = minStability;
            minStabilityDeltaValue = deltaValue;
        end
    end
    % check distance to next neighbor in y-direction
    minDistance = 0;
    minMissing = 0;
    lastMinInd = nan;
    for yi=2:size(aminInds,1)
        if isnan(aminInds(yi,:))
            % if no minima is found, this coun't as not stable
            yDiff = 2; % add some small penalty
            minMissing = minMissing + 1;
        else
            yDiff = abs(aminInds(yi,:)-aminInds(yi-1,:));
            if isnan(yDiff)
                % if last min is missing, try to use last one
                yDiff = abs(aminInds(yi,:)-lastMinInd);
                if isnan(yDiff)
                    yDiff = 2;
                    % treat as missing
                    minMissing = minMissing + 1;
                end
            end
            if ~isnan(aminInds(yi,:))
                lastMinInd = aminInds(yi,:);
            end
        end
        minDistance = minDistance + yDiff*yDiff; % penalize with square on distance
    end
    minDistance = minDistance/size(aminInds,1);

    disp(append("Found minStability with delta: ", num2str(minStabilityDeltaValue),...
        " and value: ", num2str(minStabilityMaxValue), " and distance: ", num2str(minDistance),...
        " missing: ", num2str(minMissing)));
    % minimum stability delta value: +/- num of pixels where max stability
    % occurs
    lStability.MinStabilityDeltaValue = minStabilityDeltaValue;
    % minimum stability max value: the value of maximal stability (at delta
    % value)
    lStability.MinStabilityMaxValue = minStabilityMaxValue;
    % the distance in x-direction to next minimum in y-direction
    lStability.MinDistance = minDistance;
    % how many minima were not found (e.g. because not enough points in
    % x-dir)
    lStability.MinMissing = minMissing;
else
    lStability = [];
end
tempAminInds = aminInds;
tempAminInds(isnan(tempAminInds)) = 1;
minInds=sub2ind(size(asyms_mean),(1:size(asyms_mean,1))',tempAminInds);

if printLevel > 0 && (isempty(printFilter) || printFilter=="" || printFilter==subjectNr)
    f3 = figure('Name', 'PcWithLine');
    pcshow(pc);
    hold on;
    
    scatter3(X(minInds),Y(minInds),qmaxZ(minInds), 100, 'red', '.');
    title("point cloud with external spinal line (Drerup/Hierholzer, 1986)");
    
    % TODO
    [pcFileFolderPath, pcFileName, ~] = fileparts(eslParams.FilePath);
    saveas(f3, append(pcFileFolderPath, "/Output/", pcFileName, "_eslMinima.fig"))
end

% now plot the spline resp. polynomial line
valleyPts = [X(minInds),Y(minInds),qmaxZ(minInds)];
vPtsAreNotNan = ~isnan(valleyPts);
vPtsOnly = valleyPts(vPtsAreNotNan(:,3), :);
spLine = ExternalSpinalSpline(pc, valleyPts, vPtsOnly, [], [], "(Drerup/Hierholzer, 1986)", printLevel, eslParams.FilePath);
pLine = ExternalSpinalPolynomial(pc, valleyPts, vPtsOnly, [], [], "(Drerup/Hierholzer, 1986)", printLevel, eslParams.FilePath);

[normals,curvature] = FindPointNormals(vPtsOnly,25,[0,0,0],true);
if printLevel > 1
    figure;
    surf(X,Y,tempZ);
    xlabel('x')
    ylabel('y')
    hold on;
    quiver3(vPtsOnly(:,1),vPtsOnly(:,2),vPtsOnly(:,3),...
        normals(:,1),normals(:,2),normals(:,3),'r');
    axis equal;
end
end