% Create boxplots from distance evaluation between ground truth ESL and
% calculated ESL

if nameTrans.isKey(nameFilter)
    tmpNameFilter = nameTrans(nameFilter);
else
    tmpNameFilter = nameFilter;
end

tmpNameFilter = "*" + tmpNameFilter + "*";
if ~isempty(nameFilter) && nameFilter ~= ""
    tmpNameFilter = "*";
end

outputFolderPath = baseFolderPath + "\Output\";
csvList = dir(outputFolderPath+parName + "_" + id+"\**\"+tmpNameFilter+".csv");

if size(csvList,1)>0
    DownSamples = zeros(1,size(csvList,1));
    MeanDists = [];
    MinStabilityDeltaValues = [];
    MinStabilityMaxValues = [];
    MinDistances = [];
    MinMissings = [];

    for i=1:1:size(csvList,1)
        fid = fopen(csvList(i).folder +"\"+csvList(i).name);
        txt = textscan(fid,'%s','delimiter','\n'); 
        line1 = txt{1}{1};
        line1Split = split(line1, ",");
        DownSamples(i) = str2num(line1Split{parCol});

        mDists = zeros(size(txt{1},1)-2, 1);
        mStabDeltaVals = zeros(size(txt{1},1)-2, 1);
        mStabMaxVals = zeros(size(txt{1},1)-2, 1);
        minDists = zeros(size(txt{1},1)-2, 1);
        mMissings = zeros(size(txt{1},1)-2, 1);
        for k=1:1:size(mDists,1)
            linek = txt{1}{k+2};
            linekSplit = split(linek, ",");
            mDists(k) = str2num(linekSplit{2});
            mStabDeltaVals(k) = str2num(linekSplit{6});
            mStabMaxVals(k) = str2num(linekSplit{7});
            minDists(k) = str2num(linekSplit{8});
            mMissings(k) = str2num(linekSplit{9});
        end
        MeanDists = [MeanDists, mDists];
        MinStabilityDeltaValues = [MinStabilityDeltaValues, mStabDeltaVals];
        MinStabilityMaxValues = [MinStabilityMaxValues, mStabMaxVals];
        MinDistances = [MinDistances, minDists];
        MinMissings = [MinMissings, mMissings];
    end

    % m to mm if required
    MeanDists = MeanDists*1000;

    %% mean distances
    DownSamples = DownSamples/1000;
    if size(MeanDists,1) > 0
        f1 = figure;
        boxplot(MeanDists);
        xticklabels(DownSamples);
        xlabel("cutoff frequency [1/mm]")
        title(append(parName, ": mean distances - " + nameFilter + " - " + size(MeanDists,1) + " samples"));
        saveas(f1,outputFolderPath+parName + "_" + id+"\MeanDists1.png");
        saveas(f1,outputFolderPath+parName + "_" + id+"\MeanDists1.fig");

        f2 = figure;
        plot(mean(MeanDists,1), '--o');
        xticks(1:size(DownSamples,2));
        xticklabels(DownSamples);
        xlabel("cutoff frequency [1/mm]")
        title(append(parName, ": mean distances - " + nameFilter + " - " + size(MeanDists,1) + " samples"));
        saveas(f2,outputFolderPath+parName + "_" + id+"\MeanDists2.png");
        saveas(f2,outputFolderPath+parName + "_" + id+"\MeanDists2.fig");

        % min stability delta values
        f3 = figure;
        boxplot(MinStabilityDeltaValues);
        xticklabels(DownSamples);
        xlabel("cutoff frequency [1/mm]")
        title(append(parName, ": minimum stability delta values - " + nameFilter + " - " + size(MeanDists,1) + " samples"));
        saveas(f3,outputFolderPath+parName + "_" + id+"\MinStabilityDeltaValues1.png");
        saveas(f3,outputFolderPath+parName + "_" + id+"\MinStabilityDeltaValues1.fig");

        f4 = figure;
        plot(mean(MinStabilityDeltaValues,1), '--o');
        xticks(1:size(DownSamples,2));
        xticklabels(DownSamples);
        xlabel("cutoff frequency [1/mm]")
        title(append(parName, ": minimum stability delta values - " + nameFilter + " - " + size(MeanDists,1) + " samples"));
        saveas(f4,outputFolderPath+parName + "_" + id+"\MinStabilityDeltaValues2.png");
        saveas(f4,outputFolderPath+parName + "_" + id+"\MinStabilityDeltaValues2.fig");

        % min stability max values
        f5 = figure;
        boxplot(MinStabilityMaxValues);
        xticklabels(DownSamples);
        xlabel("cutoff frequency [1/mm]")
        title(append(parName, ": minimum stability max values - " + nameFilter + " - " + size(MeanDists,1) + " samples"));
        saveas(f5,outputFolderPath+parName + "_" + id+"\MinStabilityMaxValues1.png");
        saveas(f5,outputFolderPath+parName + "_" + id+"\MinStabilityMaxValues1.fig");

        f6 = figure;
        plot(mean(MinStabilityMaxValues,1), '--o');
        xticks(1:size(DownSamples,2));
        xticklabels(DownSamples);
        xlabel("cutoff frequency [1/mm]")
        title(append(parName, ": minimum stability max values - " + nameFilter + " - " + size(MeanDists,1) + " samples"));
        saveas(f6,outputFolderPath+parName + "_" + id+"\MinStabilityMaxValues2.png");
        saveas(f6,outputFolderPath+parName + "_" + id+"\MinStabilityMaxValues2.fig");

        % min stability distance values
        f7 = figure;
        boxplot(MinDistances);
        xticklabels(DownSamples);
        xlabel("cutoff frequency [1/mm]")
        title(append(parName, ": minimum stability distance values - " + nameFilter + " - " + size(MeanDists,1) + " samples"));
        saveas(f7,outputFolderPath+parName + "_" + id+"\MinStabilityDistanceValues1.png");
        saveas(f7,outputFolderPath+parName + "_" + id+"\MinStabilityDistanceValues1.fig");

        f8 = figure;
        plot(mean(MinDistances,1), '--o');
        xticks(1:size(DownSamples,2));
        xticklabels(DownSamples);
        xlabel("cutoff frequency [1/mm]")
        title(append(parName, ": minimum stability distance values - " + nameFilter + " - " + size(MeanDists,1) + " samples"));
        saveas(f8,outputFolderPath+parName + "_" + id+"\MinStabilityDistanceValues2.png");
        saveas(f8,outputFolderPath+parName + "_" + id+"\MinStabilityDistanceValues2.fig");

        % min stability missing values
        f9 = figure;
        boxplot(MinMissings);
        xticklabels(DownSamples);
        xlabel("cutoff frequency [1/mm]")
        title(append(parName, ": minimum stability missing values - " + nameFilter + " - " + size(MeanDists,1) + " samples"));
        saveas(f9,outputFolderPath+parName + "_" + id+"\MinStabilityMissingValues1.png");
        saveas(f9,outputFolderPath+parName + "_" + id+"\MinStabilityMissingValues1.fig");

        f10 = figure;
        plot(mean(MinMissings,1), '--o');
        xticks(1:size(DownSamples,2));
        xticklabels(DownSamples);
        xlabel("cutoff frequency [1/mm]")
        title(append(parName, ": minimum stability missing values - " + nameFilter + " - " + size(MeanDists,1) + " samples"));
        saveas(f10,outputFolderPath+parName + "_" + id+"\MinStabilityMissingValues2.png");
        saveas(f10,outputFolderPath+parName + "_" + id+"\MinStabilityMissingValues2.fig");
    end
end