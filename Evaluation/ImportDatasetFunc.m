% evaluates all pointclouds in a given folder w.r.t. ESL parameters (evaluation)

% get a list of all relevant folders
dirs = dir(baseFolderPath);
dirFlags = [dirs.isdir];
dirs = dirs(dirFlags);
%dirs = dirs(1:3);

% prepare some parameters for the external spinal line extraction
eslParams = ESLParams;
eslParams.DownSample = 1;%0.5;
eslParams.SmoothRadius = 0.01;
eslParams.DiscreteSteps = 120;
eslParams.DiscreteKernel = 9;
eslParams.ManualPcRoi = 0; % 0: don't use manual roi
eslParams.LimitZResolution = 0; % 0: don't use limit z resolution
eslParams.FrequencyFilter = 0; % 0: don't use frequency filter
asnParams = AsnParams;
asnParams.NumOfPeaks = 0; % 0: no noise
asnParams.Amplitude = 0.1;
asnParams.RandomNoise = 0.001;

FolderName = append(baseFolderPath,"\","Output");   % Your destination folder
if ~exist(FolderName, 'dir')
    mkdir(FolderName);
end

BaseFolder = FolderName + "/" + EvalParamStr;
FolderName = BaseFolder;   % Your destination folder
if ~exist(FolderName, 'dir')
    mkdir(FolderName);
end

inc=0;
startTime = now;

% loop over downsampling
sNumInd=0;
for sNum = sNums
    sNumInd = sNumInd+1;
    switch EvalParamStr
        case "DownSample"
            eslParams.DownSample = sNum;
        case "SmoothRadius"
            eslParams.SmoothRadius = sNum;
        case "DiscreteSteps"
            eslParams.DiscreteSteps = sNum;
        case "DiscreteKernel"
            eslParams.DiscreteKernel = sNum;
        case "SineNoiseNumOfPeaks"
            asnParams.NumOfPeaks = sNum;
        case "SineNoiseAmplitudes"
            asnParams.Amplitude = sNum;
            % needs peak as well
            asnParams.NumOfPeaks = 5;
        case "SineNoiseRandomNoise"
            asnParams.RandomNoise = sNum;
            % needs peak as well
            asnParams.NumOfPeaks = 5;
        case "LimitZResolution"
            eslParams.LimitZResolution = sNum;
        case "FrequencyFilter"
            eslParams.FrequencyFilter = sNum;
        otherwise
            msgbox(EvalParamStr + " is not a known parameter!");
            return;
    end
    
    % switch printlevel on for first and last
    if sNumInd == 1 || sNumInd == length(sNums)
        printLevel = 1;        
    else
        %printLevel = 0;
    end
    
    if saveSummary == 1
        dstr = datestr(now,'yyyymmdd_HHMMSSFFF');

        FolderName = BaseFolder + "\" + dstr;   % Your destination folder
        if ~exist(FolderName, 'dir')
            mkdir(FolderName);
        end
        
        tmpNameFilter = nameFilter;
        if ~isempty(nameFilter) && nameFilter ~= ""
            if nameTrans.isKey(nameFilter)
                tmpNameFilter = nameTrans(nameFilter);
            else
                tmpNameFilter = nameFilter;
            end
            tmpNameFilter = "_"+tmpNameFilter;
        end
        fid = fopen(append(FolderName,'\dists'+tmpNameFilter+'.csv'),'w');
        fprintf(fid, "DownSample,"+eslParams.DownSample+",SmoothRadius," + eslParams.SmoothRadius...
            +",DiscreteSteps,"+eslParams.DiscreteSteps+",DiscreteKernel,"+eslParams.DiscreteKernel...
            +",SineNoiseNumOfPeaks,"+asnParams.NumOfPeaks+",SineNoiseAmplitude,"+asnParams.Amplitude...
            +",SineNoiseRandomNoise,"+asnParams.RandomNoise+",LimitZResolution,"+eslParams.LimitZResolution...
            +",FrequencyFilter,"+eslParams.FrequencyFilter+"\r\n");
        % save some used external spinal line parameters
        fprintf(fid, "Subject,MeanDistance [m],StdDistance [m],MinDistance [m],MaxDistance [m],"...
            +"MinStabilityDeltaValue,MinStabilityMaxValue,MinDistance,MinMissing\r\n");       
    end
    dirNums = 3:size(dirs,1);
    for dirNum = dirNums
        clearvars -except dirs dirNum printLevel saveOutput fid saveSummary eslParams FolderName sNum BaseFolder inc...
            sNums dirNums EvalParamStr startTime asnParams pcExtension nameFilter hasEOS hasReferenceSymmetryLine symmPcE_orig...
            nameTrans baseFolderPath parNameCols EvalParamStrs nameFilters sNumInd printFilter doEval
        
        maxInc = size(dirNums,2)*size(sNums,2);
        inc = inc+1;
        if inc > 1
            eta = ((now-startTime)/(inc-1)*(maxInc-inc))*24*60;
        else
            eta = 0;
        end
        disp(inc + "/" + maxInc + " eta: " + round(eta,1) + " min");
        
        if dirs(dirNum).name == "." || dirs(dirNum).name == ".."
            continue;
        end

        subjectPath = append(dirs(dirNum).folder,'\',dirs(dirNum).name);
        subjNameSplit = split(dirs(dirNum).name,"_");
        global subjectNr;
        subjectNr = convertCharsToStrings(subjNameSplit{end});
        addpath(subjectPath);

        % read the formetric stl file
        stlFilesList = [];
        for pcExInd = 1:length(pcExtension)
            fPath = subjectPath + "\*"+nameFilter+"*."+pcExtension(pcExInd);
            fPath = replace(fPath, "**", "*")
            dirRes = dir(fPath);
            for dri = 1:length(dirRes)
                if isempty(stlFilesList)
                    stlFilesList = dirRes(dri);
                else
                    stlFilesList(end+1) = dirRes(dri);
                end
            end
        end
        
        if size(stlFilesList,1) > 1
           disp('warning, multiple stl found ...'); 
        end
        for sfli = 1:length(stlFilesList)
            pause(1)
            close all;
            pause(1)
            if endsWith(stlFilesList(sfli).name, 'stl')
                pc = stlread(stlFilesList(sfli).name);
                pcLoc = pc.Points;
            else
                pc = pcread(stlFilesList(sfli).name);
                pcLoc = pc.Location;
            end
            eslParams.FilePath = append(stlFilesList(sfli).folder,"\",stlFilesList(sfli).name);

            % check if pointcloud is huge ...
            if length(pcLoc)>250000
                pc = pointCloud(pcLoc);
                ds = 250000/pc.Count;
                pcLoc = pcdownsample(pc,'random',ds).Location;
            end

            % check whether pc is in mm or m
            if max(pcLoc(:,3)) < 10
                % this must be m, otherwise we are measuring a back at 10mm
                % distance
                pcLoc = pcLoc*1000;
            end

            % normalize the z values
            pcLoc(:,3) = pcLoc(:,3)-mean(pcLoc(:,3));
            pc = pointCloud(pcLoc);

            % add sine noise
            if asnParams.NumOfPeaks > 0
                pc = AddSineNoise(pc, asnParams, printLevel);
            end

            if printLevel > 1
                figure;
                pcshow(pc);
                axis equal;
                title("subject "+subjectNr);
            end

            % my external spinal line
            pc2 = pointCloud(pc.Location);
            if printLevel > 1
                figure;
                pcshow(pc2);
            end

            % z must be positive
            pc2Loc = pc2.Location;
            pc2Loc(:,3) = pc2Loc(:,3)-pc2.ZLimits(1);

            % convert to m
            if pc2.ZLimits(2)-pc2.ZLimits(1) > 10
                pc2Loc = pc2Loc/1000;
            end
            % add distance again
            pc2Loc(:,3) = pc2Loc(:,3)+0.1;
            pc3 = pointCloud(pc2Loc);

            % Do the external spinal line detection according to Drerup/Hierholzer
            [spLine, pLine, lStability] = ExternalSpinalLineDetectionDrerupHierholzerFunc(eslParams, pc3, printLevel, 0, 1);
            
            % use baseline as reference
            if(length(sNums)<=1 || abs(sNum - sNums(1)) < 0.5*abs(sNums(1)-sNums(2)))
                symmPcE_orig{dirNum,sfli} = pLine;
            end
            symmPcE = symmPcE_orig{dirNum,sfli};
            
            % save all open figures
            if printLevel <= 1 && saveOutput == 1
                FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
                for iFig = 1:length(FigList)
                  FigHandle = FigList(iFig);
                  FigName   = subjectNr + "_" + get(FigHandle, 'Name') +"_"+num2str(dirNum)+ "_" + num2str(sfli);
                  pause(0.1)
                  savefig(FigHandle, fullfile(FolderName, FigName + '.fig'));
                  pause(0.1)
                end
            end

            if printLevel > 0  && (printFilter=="" || printFilter==subjectNr)
                fh = findobj('Type', 'Figure', 'Name', 'PcWithPolyLine');
                figure(fh);
                hold on;
            end
            if printLevel > 0  && (printFilter=="" || printFilter==subjectNr)
                pcshow(symmPcE, 'blue', 'MarkerSize', 64);
            end
            % compare to Formetric external line
            dists = zeros(1,size(pLine,1));
            for i=1:size(pLine,1)
                [ind,dist] = findNearestNeighbors(pointCloud(symmPcE), pLine(i,:), 1);
                dists(i) = dist;
            end
            distAvg = mean(dists);
            distStd = std(dists);
            distMin = min(dists);
            distMax = max(dists);

            if saveSummary == 1
                fprintf(fid, ...
                    subjectNr + "," + distAvg + "," + distStd + "," + distMin + "," + distMax + ","...
                    + lStability.MinStabilityDeltaValue + "," + lStability.MinStabilityMaxValue +","...
                    + lStability.MinDistance + "," + lStability.MinMissing + "," + stlFilesList(sfli).name + "\r\n");
            end
        end
        rmpath(subjectPath);
    end
    if saveSummary == 1
        fclose(fid);
    end
end