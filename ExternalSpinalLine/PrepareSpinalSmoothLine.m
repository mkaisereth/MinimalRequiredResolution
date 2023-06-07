function [x,y,z,t] = PrepareSpinalSmoothLine(vPtsOnly, printLevel)

% another approach - outlier removal
[vPtsOnlyInliers2, Ind] = rmoutliers(vPtsOnly(:,1), 'movmedian', 31, 'ThresholdFactor', 1);
vPtsOnlyInliers2 = [vPtsOnlyInliers2, vPtsOnly(~Ind,2), vPtsOnly(~Ind,3)];

% spline
if printLevel > 1
    fnplt(cscvn(vPtsOnlyInliers2'),'r',2)
    hold off;
end

% another approach - spline
x = vPtsOnlyInliers2(:,1);
y = vPtsOnlyInliers2(:,2);
z = vPtsOnlyInliers2(:,3);
t = 1:length(vPtsOnlyInliers2); % Assumed time stamp

end