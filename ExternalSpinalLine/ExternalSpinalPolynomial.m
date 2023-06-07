function pLine = ExternalSpinalPolynomial(pc, valleyPts, vPtsOnly, horizLine, horizInd, subtitle, printLevel, fileName)

[x,y,z,t] = PrepareSpinalSmoothLine(vPtsOnly, printLevel);

% TODO hack filter according to y coordinate
% yLimit1 = -0.24;
% yLimit2 = 0.32;
% lowerInds = y(:,1)<yLimit1;
% x(lowerInds) = [];
% y(lowerInds) = [];
% z(lowerInds) = [];
% t(lowerInds) = [];
% upperInds = y(:,1)>yLimit2;
% x(upperInds) = [];
% y(upperInds) = [];
% z(upperInds) = [];
% t(upperInds) = [];
% 
% vlowerInds = valleyPts(:,2)<yLimit1;
% valleyPts(vlowerInds, :) = [];
% vupperInds = valleyPts(:,2)>yLimit2;
% valleyPts(vupperInds, :) = [];

% another approach - polyfit
% make y evenly spaced
tmpVar = y;
deltaY = (max(tmpVar)-min(tmpVar))/length(tmpVar);
yy = min(tmpVar):deltaY:max(tmpVar)+0.1*deltaY;

px = polyfit(y,x,3);
xx = polyval(px,yy);
pz = polyfit(y,z,6);
zz = polyval(pz,yy);

pLine = [xx',yy',zz'];
if printLevel > 1
    figure
    scatter3(x,y,z)
    hold on
    plot3(xx,yy,zz-0.001)
end

global printFilter;
global subjectNr;
if printLevel > 0 && (isempty(printFilter) ||printFilter=="" || printFilter==subjectNr)
    f5 = figure('Name', 'PcWithPolyLine');
    pcshow(pc);
    hold on;
    if ~isempty(horizLine)
        pcshow(horizLine, [1 0 0]);
    end
    pcshow(valleyPts, [1 1 1]);
    plot3(xx,yy,zz-0.001, 'Color', 'r', 'LineWidth', 2)
    titleStr = "point cloud with external spinal line (polynomial)";
    if ~isempty(horizInd)
        titleStr = titleStr + "and horizontal line for h=" + horizInd;
    end
    if ~isempty(subtitle)
        titleStr = titleStr + subtitle;
    end
    title(titleStr);
    
        % TODO
    [pcFileFolderPath, pcFileName, ~] = fileparts(fileName);
    saveas(f5, append(pcFileFolderPath, "/Output/", pcFileName, "_eslPoly.fig"))
end
end