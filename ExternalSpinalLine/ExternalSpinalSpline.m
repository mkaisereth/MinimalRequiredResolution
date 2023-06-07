function spLine = ExternalSpinalSpline(pc, valleyPts, vPtsOnly, horizLine, horizInd, subtitle, printLevel, fileName)

[x,y,z,t] = PrepareSpinalSmoothLine(vPtsOnly, printLevel);

% Apply interpolation for each x,y and z 
tt = linspace(t(1),t(end));
%xx = interp1(t,x,tt,'spline');
%yy = interp1(t,y,tt,'spline');
%zz = interp1(t,z,tt,'spline');
xx = spline(t,x,tt);
yy = spline(t,y,tt);
zz = spline(t,z,tt);

spLine = [xx',yy',zz'];
% Visualize the result
if printLevel > 1
    figure
    scatter3(x,y,z)
    hold on
    plot3(xx,yy,zz-0.001)
end

global printFilter;
global subjectNr;
if printLevel > 0 && (isempty(printFilter) ||printFilter=="" || printFilter==subjectNr)
    f4 = figure('Name', 'PcWithSpline');
    pcshow(pc);
    hold on;
    if ~isempty(horizLine)
        pcshow(horizLine, [1 0 0]);
    end
    pcshow(valleyPts, [1 1 1]);
    plot3(xx,yy,zz-0.001, 'Color', 'r', 'LineWidth', 2)
    titleStr = "point cloud with external spinal line (spline)";
    if ~isempty(horizInd)
        titleStr = titleStr + "and horizontal line for h=" + horizInd;
    end
    if ~isempty(subtitle)
        titleStr = titleStr + subtitle;
    end
    title(titleStr);
    
        % TODO
    [pcFileFolderPath, pcFileName, ~] = fileparts(fileName);
    saveas(f4, append(pcFileFolderPath, "/Output/", pcFileName, "_eslSpline.fig"))
end

end