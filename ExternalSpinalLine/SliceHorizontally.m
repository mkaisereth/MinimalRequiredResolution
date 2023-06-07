% function to slice a pointcloud in horizontal slices
function [slices, xDistAvg, yDistAvg] = SliceHorizontally(pc, pcLoc, printLevel)

% get average distance horizontally and vertically
    xDistAvg = [];
    yDistAvg = [];
    for i=1:length(pcLoc)
        % first find vertical (y)
        numNeigbors = 2;
        foundy=0;
        while 1
            [inds, ~] = findNearestNeighbors(pc, pcLoc(i,:), numNeigbors, 'Sort', true);
            for nn=2:numNeigbors
                xDist = abs(pcLoc(i,1)-pc.Location(inds(nn),1));
                yDist = abs(pcLoc(i,2)-pc.Location(inds(nn),2));
                zDist = abs(pcLoc(i,3)-pc.Location(inds(nn),3));
                if yDist >= xDist %&& yDist >= zDist
                    yDistAvg(end+1) = yDist;
                    foundy=1;
                    break; % found one
                end
            end
            if foundy==1
                break;
            end
            numNeigbors = numNeigbors+2;
        end

        % first find vertical (y)
        numNeigbors = 2;
        foundx=0;
        while 1
            [inds, ~] = findNearestNeighbors(pc, pcLoc(i,:), numNeigbors, 'Sort', true);
            for nn=2:numNeigbors
                xDist = abs(pcLoc(i,1)-pc.Location(inds(nn),1));
                yDist = abs(pcLoc(i,2)-pc.Location(inds(nn),2));
                zDist = abs(pcLoc(i,3)-pc.Location(inds(nn),3));
                if xDist >= yDist %&& xDist >= zDist
                    xDistAvg(end+1) = xDist;
                    foundx=1;
                    break;
                end
            end
            if foundx==1
                break;
            end
            numNeigbors = numNeigbors+2;
        end
        
    end
    xDistAvg = mean(xDistAvg);
    yDistAvg = mean(yDistAvg);
    
    % slice in horizontal lines
    yMin = pc.YLimits(1);
    yMax = pc.YLimits(2);
    
    yVals = pc.Location(abs(pc.Location(:,2)-yMin)<0.5*yDistAvg,:);
    yValsMin = mean(yVals(:,2));
    
    counter2=1;
    totalCount = size(yValsMin:yDistAvg:yMax, 2);
    slices = {};
    for yVal = yValsMin:yDistAvg:yMax
        yVals = pc.Location(abs(pc.Location(:,2)-yVal)<0.5*yDistAvg,:);
        
        if length(yVals)==0
            yVals = pc.Location(abs(pc.Location(:,2)-yVal-0.2*yDistAvg)<0.5*yDistAvg,:);
        end
        
        [~,i] = sort(yVals(:,1)); % sort according to x
        yVals = yVals(i,:);
        
        % check again according to mean
        yMean = mean(yVals(:,2));
        yVals = pc.Location(abs(pc.Location(:,2)-yMean)<0.5*yDistAvg,:);
        [~,i] = sort(yVals(:,1)); % sort according to x
        yVals = yVals(i,:);
        
        if printLevel > 3 || (printLevel > 2 && counter2 == floor(totalCount/3))
            figure;
            pcshow(yVals,'MarkerSize', 34);
        end
        
        % sanity check
        if(size(yVals,1)<2)
           disp('here'); 
        end
        for i=2:size(yVals,1)
            if yVals(i,1)-yVals(i-1,1)>2*xDistAvg
                disp('missed a point!')
            end
        end
        
        yzVals = yVals(:,3);
        t = yVals(:,1);
        if printLevel > 3 || (printLevel > 2 && counter2 == floor(totalCount/3))
            figure;
            plot(t,yzVals);
        end
        
        slices{end+1} = [yVals(:,1), yVals(:,2), yzVals(:,1)];
        
        counter2 =counter2+1;
    end

end