% evaluates all pointclouds in given folders w.r.t. ESL parameters (evaluation)
close all; clear all; clc;

doEval = true;

printLevel = 3;
saveOutput = 1;
saveSummary = 1;
global printFilter;
printFilter = "154531966";

baseFolderPath = "C:\Users\ETH\Downloads\Upright20220721_sample";

addpath('ExternalSpinalLine');
addpath('Evaluation');

%pcExtension = "stl";
pcExtension = ["ply", "pcd"];

%nameFilters = ["Photoneo", "point_cloud", "AS_smooth", "Astra", "D415", "ADIToF"];
nameFilters = ["Photoneo"];

hasEOS = 0;
hasReferenceSymmetryLine = 0;

nameTrans = containers.Map;
nameTrans("point_cloud") = "SLTida";
nameTrans("AS_smooth") = "vltAS";

parNameCols = containers.Map;
parNameCols("DownSample") = 2;
parNameCols("SmoothRadius") = 4;
parNameCols("DiscreteSteps") = 6;
parNameCols("DiscreteKernel") = 8;
parNameCols("SineNoiseNumOfPeaks") = 10;
parNameCols("SineNoiseAmplitudes") = 12;
parNameCols("SineNoiseRandomNoise") = 14;
parNameCols("LimitZResolution") = 16;
parNameCols("FrequencyFilter") = 18;

EvalParamStrs = containers.Map;
%EvalParamStrs("DownSample") = 1:-0.1:0.1;
%EvalParamStrs("SmoothRadius") = [0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2];
%EvalParamStrs("DiscreteSteps") = [480, 240,120,100,80,60,40];
%EvalParamStrs("DiscreteKernel") = [1,3,5,9,13,21,29];
%EvalParamStrs("SineNoiseNumOfPeaks") = 1:20;
%EvalParamStrs("SineNoiseAmplitudes") = [0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2];
%EvalParamStrs("SineNoiseRandomNoise") = [0.05,0.1,0.5,1,1.5,2];
%EvalParamStrs("LimitZResolution") = [0.001,0.002,0.005,0.01,0.02];
%EvalParamStrs("FrequencyFilter") = [0.076, 0.05,0.02,0.01,0.009, 0.0076]*1000; % PrepareTestPointCloud is in m!
EvalParamStrs("FrequencyFilter") = [0.1, 0.05,0.02,0.01,0.005]*1000; % PrepareTestPointCloud is in m!

for EvalParamStrInd = 1:length(EvalParamStrs)
    kys = EvalParamStrs.keys;
    EvalParamStr = kys{EvalParamStrInd};
    sNums = EvalParamStrs(EvalParamStr);
    for nameFilter = nameFilters
        
        if doEval
            ImportDatasetFunc;
        end

        % copy the folder
        FolderName = append(baseFolderPath,"\","Output");   % Your destination folder
        BaseFolder = FolderName + "/" + EvalParamStr;
        fNum = 1;
        dstr = datestr(now,'yyyymmdd');
        newFolder = append(BaseFolder, "_", dstr);
        newFolderwNum = append(newFolder, "_", num2str(fNum));
        if doEval
            while exist(newFolderwNum)
                fNum = fNum+1;
                newFolderwNum = append(newFolder, "_", num2str(fNum));
            end
            movefile(BaseFolder,newFolderwNum);
        end
        % do the eval of boxplots
        parName = EvalParamStr;
        parCol = parNameCols(parName);
        id = append(dstr, "_", num2str(fNum));
        BoxPlotDistsFunc;
    end
end