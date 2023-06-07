% fourier transform analysis (frequency analysis) of horizontal slices of point clouds
close all; clear all; clc;

dt = datetime('now');
dt.Format = 'yyyyMMddHHmmss';

printLevel = 1;
doFilters = [0 1 1 1 1 1 1 1 1 1 1 1]; % 1: filter with butterworth
filterFreqs = [0 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005]; % 0.01 [1/mm]
% photoneo average filter frequency 1.3 => filterFreq/0.65 == [0,1]% => [0-0.65]
% => filterFreq=0.1 ==> 0.15%

baseFolderPath = "C:\Users\ETH\Downloads\Upright20220721_sample";

%pcExtension = "stl";
pcExtension = ["ply", "pcd"];

%nameFilter = "";
%nameFilter = "D415";
%nameFilter = "point_cloud";
nameFilter = "Photoneo";
%nameFilter = "AS_smooth";
%nameFilter = "ADIToF"
%nameFilter = "Astra"

figuresTitle = "Upright sample";

% get a list of all relevant folders
%dirs = [dir("Milano Datasets - ScolioSIM 15-10-2021"); dir("Milano Datasets - ScolioSIM 29-11-2021")];
dirs = dir(baseFolderPath);
dirFlags = [dirs.isdir];
dirs = dirs(dirFlags);

dirNums = 3:size(dirs,1);

delete(baseFolderPath + "\ErrorToOrig.txt");
delete(baseFolderPath + "\ErrorToOrig_MAE.txt");
needDeleteErrorFile = {};

num = 0;
numTotal = length(dirNums)*length(filterFreqs);
startTime = now;

for filterFreqsId=1:length(filterFreqs)
    close all;
    doFilter = doFilters(filterFreqsId);
    filterFreq = filterFreqs(filterFreqsId);
    
    filterString = "";
    if doFilter == 1
        filterString = "_" + string(filterFreq);
    end
    
    zMeans = nan(1,size(dirNums,2));
    zStds = nan(1,size(dirNums,2));
    zMins = nan(1,size(dirNums,2));
    zMaxs = nan(1,size(dirNums,2));
    relZMins = nan(1,size(dirNums,2));
    relZMaxs = nan(1,size(dirNums,2));
    
    yMags = {};
    p2s = {};
    errorsToOriginalSAE = {};
    errorsToOriginalMAE = {};
    errorsToOriginalNums = {};
    counter = 1;
    for dirNum = dirNums
        disp(num2str(dirNum));
        if dirs(dirNum).name == "." || dirs(dirNum).name == ".."
            continue;
        end

        % estimate eta
        num = num+1;
        if num > 1
            eta = ((now-startTime)/(num-1)*(numTotal-num))*24*60;
        else
            eta = 0;
        end
        disp(num + "/" + numTotal + " eta: " + round(eta,1) + " min");
    
        subjectPath = append(dirs(dirNum).folder,'\',dirs(dirNum).name);
        subjNameSplit = split(dirs(dirNum).name,"_");
        subjectNr = subjNameSplit(end);
        addpath(subjectPath);

        if length(needDeleteErrorFile)<dirNum || isempty(needDeleteErrorFile{dirNum})
            delete(subjectPath + "\ErrorToOrig.txt");
            delete(subjectPath + "\ErrorToOrig_MAE.txt");
            needDeleteErrorFile{dirNum} = 0;
        end
    
        % read the formetric stl file
        stlFilesList = [];
        for pcExInd = 1:length(pcExtension)
            fPath = subjectPath + "\*"+nameFilter+"*."+pcExtension(pcExInd)
            fPath = replace(fPath, "**", "*");
            dirRes = dir(fPath);
            for dri = 1:length(dirRes)
                if isempty(stlFilesList)
                    stlFilesList = dirRes(dri);
                else
                    stlFilesList(end+1) = dirRes(dri);
                end
            end
        end
        if length(stlFilesList) > 1
           disp('warning, multiple stl found ...'); 
        end
        for sfli = 1:length(stlFilesList)
            close all;
            if endsWith(stlFilesList(sfli).name, 'stl')
                pc = stlread(stlFilesList(sfli).name);
                pcLoc = pc.Points;
            else
                pc = pcread(stlFilesList(sfli).name);
                pcLoc = pc.Location;
            end
    
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
            
            filteredPc = [];
            if printLevel > 1
                figure;
                pcshow(pcLoc);
                axis equal;
                title("subject "+subjectNr);
            end
    
            pc = pointCloud(pcLoc);
    
            % get horizontal slices
            [slices, xDistAvg, yDistAvg] = SliceHorizontally(pc, pcLoc, printLevel);
    
            % save the frequency domain values for later
            yMags{counter} = {};
            p2s{counter} = {};
            errorsToOriginalSAE{counter} = {};
            errorsToOriginalMAE{counter} = {};
            errorsToOriginalNums{counter} = 0;
    
            counter2=1;
            for sliceInd = 1:size(slices,2)
    
                yVals = slices{sliceInd};
                yzVals = yVals(:,3);
    
                if size(yzVals,1) <= 30
                    disp("not enough data: " + sliceInd);
                    continue;
                end
                % sanity check
                if max(diff(yVals(:,1))) > 15*xDistAvg
                    disp("xDistAvg sanity check failed: " + sliceInd);
                    continue;
                end
    
                % make it regular
                queryPoints = min(yVals(:,1)):xDistAvg:max(yVals(:,1));
                [yValsXUnique,ia] = unique(yVals(:,1), 'stable');
                yValsYUnique = yVals(ia,3);
                vq2 = interp1(yValsXUnique,yValsYUnique,queryPoints);
    
                %sanity check
                if max(vq2) > 1000 || min(vq2) < -1000
                    disp('fail');
                end
    
                yMeanTmp = mean(yVals(:,2));
                yValsTemp = [queryPoints', yMeanTmp*ones(length(vq2), 1), vq2'];
                yVals = yValsTemp;
                yzVals = yVals(:,3);
    
                % save original for error calculation later
                yzValsOriginal = yzVals;
    
                % filter (if applicable)
                if doFilter == 1
                    Ts = xDistAvg; % sampling period [mm]
                    Fs = 1/Ts; % sampling frequency [1/mm]
                    fc = filterFreq;
                    %fs = 1000;
                    %normalized cutoff frequency
                    wn = fc/(Fs/2);
                    [b,a] = butter(6,wn);
    
                    if printLevel > 3 || (printLevel > 2 && sliceInd == floor(size(slices,2)/3))
                        figure;
                        freqz(b,a)
                    end
                    if size(yzVals,1) <= 30
                        disp("not enough data for filtering: " + sliceInd);
                    else
                        yzVals = filtfilt(b,a, double(yzVals));
                    end
                    filteredPc = [filteredPc;[yVals(:,1),yVals(:,2),yzVals(:,1)]];
                    t = yVals(:,1);
                    if printLevel > 3  || (printLevel > 2 && sliceInd == floor(size(slices,2)/3))
                        figure;
                        plot(t,yzVals);
                    end
                end
                % do a fourier transform
                if ~isempty(yzVals)
                    Ts = xDistAvg; % sampling period [mm]
                    Fs = 1/Ts; % sampling frequency [1/mm]
                    L = size(yzVals,1); % length
                    Y = fft(yzVals); % fourier transform
                    f = (0:L-1)*Fs/L;        
                    n = L;
                    if length(1:n/2+1)> length(f)
                        disp('here')
                    end
                    f = f(1:n/2+1);
                    absY = abs(Y(1:n/2+1));
                    fshift = (-n/2:n/2-1)*(Fs/n);
                    yshift = fftshift(Y);
    
                    yMags{counter}{counter2}{1} = f;
                    yMags{counter}{counter2}{2} = absY;
                    errorsToOriginalSAE{counter}{counter2}{1} = sum(abs(yzValsOriginal-yzVals));
                    errorsToOriginalMAE{counter}{counter2}{1} = sum(abs(yzValsOriginal-yzVals));
                    errorsToOriginalNums{counter} = errorsToOriginalNums{counter} + size(abs(yzValsOriginal-yzVals),1);
    
                    if printLevel > 3 || (printLevel > 2 && sliceInd == floor(size(slices,2)/3))
                        figure;
                        plot(f,absY)
                        title("Fouriertransform of a single slice");
                        xlabel('Frequency (1/mm)')
                        ylabel('Magnitude')
                    end
                    % another approach
                    n = 2^nextpow2(L);
                    Y = fft(yzVals,n);
                    f = Fs*(0:(n/2))/n;
                    P = abs(Y/n);
                    p2s{counter}{counter2}{1} = f;
                    p2s{counter}{counter2}{2} = P(1:n/2+1);
    
                    if printLevel > 3 || (printLevel > 2 && sliceInd == floor(size(slices,2)/3))
                        figure;
                        plot(f,P(1:n/2+1))
                        title("Fouriertransform of a single slice");
                        xlabel('Frequency (1/mm)')
                        ylabel('|P(f)|')
                    end
    
                    if counter2 == 1
                        disp('huhu')
                    end
                    counter2 =counter2+1;
                end
            end
    
            % calculate error for overall subject
            sumError = 0;
            % loop all slices
            for i=1:length(errorsToOriginalSAE{counter})
                sumError = sumError + errorsToOriginalSAE{counter}{i}{1};
            end
            % divide SAE by slices
            errorsToOriginalSAE{counter} = sumError/length(errorsToOriginalSAE{counter});
            % divide SAE (per slice) by number of points per slice => MAE per point
            errorsToOriginalMAE{counter} = errorsToOriginalSAE{counter}/max(1,errorsToOriginalNums{counter});
            % save it to csv
            fid = fopen(subjectPath + "\ErrorToOrig.txt", 'a+');
            fprintf(fid, string(doFilter) + "," + filterFreq + "," + string(errorsToOriginalSAE{counter}) + "\n");
            fclose(fid);
            fid2 = fopen(subjectPath + "\ErrorToOrig_MAE.txt", 'a+');
            fprintf(fid2, string(doFilter) + "," + filterFreq + "," + string(errorsToOriginalMAE{counter}) + "\n");
            fclose(fid2);
    
            if doFilter == 1 && printLevel > 1
                figure;
                pcshow(filteredPc);
            end
    
            % show summary for this subject
            allFis = [];
            for tempc2 = 1:length(yMags{counter})
                for i=1:length(yMags{counter}{tempc2}{1})
                    fi = yMags{counter}{tempc2}{1}(i);
                    % check whether already in there
                    alreadyContained = 0;
                    for fisInd=1:length(allFis)
                        if abs(allFis(fisInd)-fi)<0.0001
                            alreadyContained = 1;
                            break;
                        end
                    end
                    if alreadyContained == 0
                        allFis(end+1) = fi;
                    end
                end
            end
            allFis = sort(allFis);
    
            allyMags = zeros(length(allFis),1);
            allnumOfMags = zeros(length(allFis),1);
            % now go trough again and assign fises
            for tempc2 = 1:length(yMags{counter})
                for i=1:length(yMags{counter}{tempc2}{1})
                    fi = yMags{counter}{tempc2}{1}(i);
                    for fisInd=1:length(allFis)
                        if abs(allFis(fisInd)-fi)<0.0001
                            allnumOfMags(fisInd) = allnumOfMags(fisInd)+1;
                            allyMags(fisInd) = allyMags(fisInd)+yMags{counter}{tempc2}{2}(i);
                            break;
                        end
                    end
                end
            end
            allyMagsAvg = zeros(length(allyMags),1);
            for i=1:length(allyMags)
                allyMagsAvg(i) = allyMags(i)/allnumOfMags(i);
            end
            if dirNum == 6
                disp('here')
            end
            if printLevel > 1
                fig9 = figure;
                plot(allFis, allyMagsAvg);
                title("Subject "+subjectNr + " allyMagsAvg");
                xlabel("Frequency [1/mm]");
                ylabel("Magnitude");
                saveas(fig9, subjectPath + "\allyMagsAvg"+filterString+".fig");
            end
    
            counter=counter+1;
        end
        rmpath(subjectPath);
        if contains(subjectNr, "ksm")
            disp("ksm")
        end
    end
    
    %% show overall summary
    allFis = [];
    for tempc =1:length(yMags)
        for tempc2 = 1:length(yMags{tempc})
            for i=1:length(yMags{tempc}{tempc2}{1})
                fi = yMags{tempc}{tempc2}{1}(i);
                % check whether already in there
                alreadyContained = 0;
                for fisInd=1:length(allFis)
                    if abs(allFis(fisInd)-fi)<0.0001
                        alreadyContained = 1;
                        break;
                    end
                end
                if alreadyContained == 0
                    allFis(end+1) = fi;
                end
            end
        end
    end
    allFis = sort(allFis);
        
    allyMags = zeros(length(allFis),1);
    allnumOfMags = zeros(length(allFis),1);
    % now go trough again and assign fises
    for tempc = 1:length(yMags)
        for tempc2 = 1:length(yMags{tempc})
            for i=1:length(yMags{tempc}{tempc2}{1})
                fi = yMags{tempc}{tempc2}{1}(i);
                for fisInd=1:length(allFis)
                    if abs(allFis(fisInd)-fi)<0.0001
                        allnumOfMags(fisInd) = allnumOfMags(fisInd)+1;
                        allyMags(fisInd) = allyMags(fisInd) + yMags{tempc}{tempc2}{2}(i);
                        break;
                    end
                end
            end
        end
    end
    for i=1:length(allyMags)
        allyMagsAvg(i) = allyMags(i)/allnumOfMags(i);
    end
    
    if printLevel > 0
        fig9 = figure;
        plot(allFis, allyMagsAvg);
        title(figuresTitle + " - " + length(yMags) + " samples");
        xlabel("Frequency [1/mm]");
        ylabel("Magnitude");
        saveas(fig9, baseFolderPath + "\allyMagsAvg_samples"+filterString+".fig");
    end

    % save overall error to csv
    overallError = 0;
    for i=1:length(errorsToOriginalSAE)
        overallError = overallError + errorsToOriginalSAE{i};
    end
    overallError = overallError/length(errorsToOriginalSAE);

    fid = fopen(baseFolderPath + "\ErrorToOrig.txt", 'a+');
    fprintf(fid, string(doFilter) + "," + filterFreq + "," + string(overallError) + "\n");
    fclose(fid);
    % MAE
    overallError = 0;
    for i=1:length(errorsToOriginalMAE)
        overallError = overallError + errorsToOriginalMAE{i};
    end
    overallError = overallError/length(errorsToOriginalMAE);

    fid2 = fopen(baseFolderPath + "\ErrorToOrig_MAE.txt", 'a+');
    fprintf(fid2, string(doFilter) + "," + filterFreq + "," + string(overallError) + "\n");
    fclose(fid2);
end