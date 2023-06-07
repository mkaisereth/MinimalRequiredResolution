close all; clear all; clc;

baseFolderPath = "C:\Users\ETH\Downloads\Upright20220721_sample";

% get a list of all relevant folders
dirs = dir(baseFolderPath);
dirFlags = [dirs.isdir];
dirs = dirs(dirFlags);

dirNums = 3:size(dirs,1);

numOfSubjects = 0;
for dirNum = dirNums
    disp(num2str(dirNum));
    if dirs(dirNum).name == "." || dirs(dirNum).name == ".." || strcmp(dirs(dirNum).name, "Output")
        continue;
    end

    subjectPath = append(dirs(dirNum).folder,'\',dirs(dirNum).name);
    subjNameSplit = split(dirs(dirNum).name,"_");
    subjectNr = subjNameSplit(end);
    addpath(subjectPath);

    % read the .txt file
    fileContent = readmatrix(subjectPath + "\ErrorToOrig_MAE.txt");

    figure;
    plot(fileContent(:,3));
    xticks(1:length(fileContent(:,3)));
    xticklabels(fileContent(:,2));
    title("Error to original shape");
    ylabel('Error');
    xlabel('cutoff frequency [1/mm]');
    numOfSubjects = numOfSubjects+1;
end

% read the summary .txt file
fileContent = readmatrix(baseFolderPath + "\ErrorToOrig_MAE.txt");

% convert m to mm
fileContent(:,3) = fileContent(:,3)*1000;

fig9 = figure;
plot(fileContent(:,3));
xticks(1:length(fileContent(:,3)));
xticklabels(fileContent(:,2));
title("Error to original shape - " + string(numOfSubjects) + " subjects");
ylabel('Error [mm]');
xlabel('cutoff frequency [1/mm]');

saveas(fig9, baseFolderPath + "\ErrorToOrig_MAE.fig");