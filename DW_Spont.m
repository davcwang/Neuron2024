function [spontMat] = DW_Spont(Expt,NPXSpikes)

% Script for analyzing points in neural state space 
% 1) Open NPXSpikes or NPXSpikesSpont for synced probes

% Initialize parameters
Fs = 30000; % sampling rate
bin = 0.010  * Fs; % bin time in seconds * Fs = bin time in samples

% Define shortcuts
ss = NPXSpikes.ss; % array of spike sample times for all spikes
clu = NPXSpikes.clu; % array of cids for each spike in ss
cids = NPXSpikes.cids'; % array of sorted good unit cids
freq = NPXSpikes.CluFreq;
depth = NPXSpikes.CluDepth;

% ExptDur = ss(end); % duration (in samples) of whole recording, defined by the time of the last recorded spike
ExptDur = 15 * 60 * Fs; % 15min spont act

% Get TRAPidx
[TRAPcids, TRAPidx, age, ExptNo] = DW_GetTRAPcidsidx(Expt,NPXSpikes); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get binned and smoothed spike rate during spont epoch

% Create array of spike counts vs bins
unitSpikeMat = zeros(length(NPXSpikes.cids), ExptDur/bin); % create matrix where rows are units, columns are each time bin

% Populate unitSpikeMat over ExptDur
for i = 1:bin:ExptDur-bin % loop through first index (sample) in each bin (here it's 100ms, so each i corresponds to 100ms, 200ms, etc) until the 2nd last bin (did this bc there was some random error if I went all the way to ExptDur)

    binNum = round(i / bin) + 1; % this is the bin number, for assigning to the correct unitSpikeMat column

    % Find all spikes within this interval
    firstSpike = find( ss >= i, 1); % ss index of first spike within this interval
    lastSpike = find(ss <= i+bin, 1, 'last'); % ss index of last spike within this interval
    intervalSpikesSS = ss(firstSpike:lastSpike); % create an array of spike times within this interval
    intervalSpikesCLU = clu(firstSpike:lastSpike); % create an array of unsorted cluster IDs corresponding to each spike within this interval

    for j = 1:length(intervalSpikesSS) % loop through spikes within interval

        idx = find(cids == intervalSpikesCLU(j)); % find idx of each spike, since this find() returns the index in cids of the current clu, which is the idx

        if ~isempty(idx) % idx will only = [] if the current clu ID is not in cids (ie. it's not a good curated unit)
            
            unitSpikeMat(idx,binNum) = unitSpikeMat(idx,binNum) + 1; % increment the number of spikes belonging to that unit within this bin

        end
    end
end

% Convolve each row of unitSpikeMat
convSpikeMat = zeros(length(NPXSpikes.cids), ExptDur/bin);

for i = 1:length(unitSpikeMat(:,1)) % loop over each cid

    convSpikeMat(i,:) = smooth(unitSpikeMat(i,:),10,'lowess'); % unclear if the conv is doing the right thing, since conv is over already binned data (so if convBin is std 300, it's not actually 0.01s but it's 300 bins...?)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate smoothed z-score spike rate

% Create a z-score matrix from the binned spike counts
ZConvSpikeMat = zeros(size(convSpikeMat));
ZUnitSpikeMat = zeros(size(unitSpikeMat));

for i = 1:length(convSpikeMat(:,1))
    
    % Get Z-score for smoothed mat
    mu = mean(convSpikeMat(i,:));
    sigma = std(convSpikeMat(i,:));
    ZConvSpikeMat(i,:) = (convSpikeMat(i,:) - mu) / sigma;

    % Get Z-score for non-smoothed mat
    mu = mean(unitSpikeMat(i,:));
    sigma = std(unitSpikeMat(i,:));
    ZUnitSpikeMat(i,:) = (unitSpikeMat(i,:) - mu) / sigma;

end

% Sort cells based on similarity in firing in ZConvSpikeMat
T = clusterdata(ZConvSpikeMat,'Criterion','distance','Cutoff',0.1);
T(:,2) = 1:length(T(:,1)); % make the second column the idx of the cells for easy reference later
ZConvTemp = zeros(length(ZConvSpikeMat(:,1)),length(ZConvSpikeMat(1,:))+1); % make temp matrix of unit# x bin# + 1 (+1 for column to add cluster IDs)
ZConvTemp(:,1) = T(:,1); % make a new matrix where the first column is cluster ID from clusterdata
ZConvTemp(:,2:end) = ZConvSpikeMat; % populate rest of columns with Z scores
ZConvSorted = sortrows(ZConvTemp); % sort z-score matrix by cluster ID
idxSorted = sortrows(T); % sort T the same way, now the 2nd column is the unit idx sorted in the same way the z-score matrix was sorted
idxSorted = idxSorted(:,2); % only keep the 2nd column, now can sort other stuff by this (eg. ZScoreMatP) to look at if these cells cluster together
ZConvSorted = ZConvSorted(:,2:end); % get rid of 1st row that's currently storing the cluster IDs

% Plot as heatmap
testrangesec = [1 800]; % enter seconds of spont recording to visualize
testrange = (testrangesec(1)*30000/bin):(testrangesec(2)*30000/bin);
figure
hold on
subplot(2,1,1)
imagesc(ZConvSorted(:,testrange))
colormap(colormap_BlueWhiteRed);
title(strcat("Spontaneous period Z-scores, sorted, ",int2str(bin/Fs*1000),"ms bins"))
caxis([-2 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Participation in population peaks

% Get rid of silent cells during spont act to prevent mean for popAct returning NaNs
ZConvSorted(isnan(ZConvSorted)) = 0;

% Initiate spont mat
% Col 1:3 = distToHighPk, distToLowPk, distToAllPk (these are total distances)
% Col 4:6 = avgToHighPk, avgToLowPk, avgToAllPk
% Col 7:9 = highCoPk, lowCoPk, allCoPk
spontMat = zeros(length(NPXSpikes.cids),9);
unitPksCell = {}; % for storing pks of individual units

% Get peaks of individual units
smoothBin = 20;
unitPksCell = DW_GetUnitPks(NPXSpikes, ZConvSorted, smoothBin, unitPksCell); % unitPksCell: each row = sorted row (from ZConvSorted), col 1 = pks, 2 = locs, 3 = FWHM

% Generate average population activity and find peaks, locations, FWHM from population activity
ZpopAct = mean(ZConvSorted,1); % average population activity
smoothZpopAct = smooth(ZpopAct,smoothBin,"lowess"); % 
stdPopAct = std(smoothZpopAct); % for calculating peaks that are a certain std above baseline

% Find population high peaks (>30% cells) and remove those under a certain ensemble participation
threshold = stdPopAct + median(smoothZpopAct); % find indices of high pks > std above median
[popHighPks, popHighLocs, popHighWidths] = DW_GetPopPks(smoothZpopAct,threshold); % find pks 
unitsInHighPk = zeros(length(NPXSpikes.cids),length(popHighPks)); % for storing all the units that fire within each pk (then divide by pks later)
[spontMat, unitsInHighPk] = DW_GetUnitCoPks(NPXSpikes,spontMat,unitPksCell,popHighPks,popHighLocs,popHighWidths,unitsInHighPk, 1, 7); % find unitsInHighPk to threshold based on ensemble participation; increment dist count in col 1 and coPk count in col 7
avgUnitsInHighPk = mean(unitsInHighPk,1); % get avgUnitsInPk as fraction of total units active in 
ensemblesize = 0.30; % set threshold ensemble participation for high pk as 30% of units
popHighPks = popHighPks(find(avgUnitsInHighPk > ensemblesize)); % only keep Pks above ensemlesize
popHighLocs = popHighLocs(find(avgUnitsInHighPk > ensemblesize));
popHighWidths = popHighWidths(find(avgUnitsInHighPk > ensemblesize));

% Find population any peaks and remove those already found in high peaks such that only "low peaks" remain
threshold = * stdPopAct + median(smoothZpopAct); % find indices of high pks > std above median
[popLowPks, popLowLocs, popLowWidths] = DW_GetPopPks(smoothZpopAct,threshold); % find pks 
unitsInLowPk = zeros(length(NPXSpikes.cids),length(popLowPks)); % for storing all the units that fire within each pk (then divide by pks later)
popLowPks( arrayfun(@(y)find(popLowPks == y), popHighPks) ) = []; % remove popLowPks already in popHighPks
popLowLocs( arrayfun(@(y)find(popLowLocs == y), popHighLocs) ) = [];
popLowWidths( arrayfun(@(y)find(popLowWidths == y), popHighWidths) ) = [];

% Re-initialize spontMat and unintsInPks (since they're filled in with trash rn) and fill with high and low pk data
spontMat = zeros(length(NPXSpikes.cids),9);
unitsInHighPk = zeros(length(NPXSpikes.cids),length(popHighPks)); % for storing all the units that fire within each pk (then divide by pks later)
unitsInLowPk = zeros(length(NPXSpikes.cids),length(popLowPks)); % for storing all the units that fire within each pk (then divide by pks later)
[spontMat, unitsInHighPk] = DW_GetUnitCoPks(NPXSpikes,spontMat,unitPksCell,popHighPks,popHighLocs,popHighWidths,unitsInHighPk, 1, 7); % find unitsInHighPk to threshold based on ensemble participation; increment dist count in col 1 and coPk count in col 7
[spontMat, unitsInLowPk] = DW_GetUnitCoPks(NPXSpikes,spontMat,unitPksCell,popLowPks,popLowLocs,popLowWidths,unitsInLowPk, 2, 8); % find unitsInLowPk to threshold based on ensemble participation; increment dist count in col 2 and coPk count in col 8

% Convert from #bin to sec
spontMat(:,1) = spontMat(:,1) * (bin/Fs); % convert from bin to s
spontMat(:,2) = spontMat(:,2) * (bin/Fs); % convert from bin to s

% Get combined metrics for distToPk and CoPk
spontMat(:,3) = spontMat(:,1) + spontMat(:,2); % col 3 = distToAllPks
spontMat(:,9) = spontMat(:,7) + spontMat(:,8); % col 9 = allCoPk

% Calculate avgDistToPeak
spontMat(:,4:6) = spontMat(:,1:3) ./ spontMat(:,7:9);

% Normalize CoPks to fraction of total high or low peaks
spontMat(:,7) = spontMat(:,7) / length(popHighPks);
spontMat(:,8) = spontMat(:,8) / length(popLowPks);
spontMat(:,9) = spontMat(:,9) / (length(popHighPks) + length(popLowPks));

% Normalize unitsInPk by total number of peaks
avgUnitsInHighPk = mean(unitsInHighPk,1); % normalize by number of units to get what percentage of total units particpated in each pop pk
avgUnitsInLowPk = mean(unitsInLowPk,1); % normalize by number of units to get what percentage of total units particpated in each pop pk

% Resort spontMat, unitsInPks by idxSorted, then get rid of that idx column
spontMat(:,10) = idxSorted;
spontMat = sortrows(spontMat,10);
spontMat = spontMat(:,1:9);
spontMat(isnan(spontMat)) = 0;
unitsInHighPk(:,end+1) = idxSorted;
unitsInHighPk = sortrows(unitsInHighPk,size(unitsInHighPk,2));
unitsInHighPk = unitsInHighPk(:,1:(end-1));
unitsInLowPk(:,end+1) = idxSorted;
unitsInLowPk = sortrows(unitsInLowPk,size(unitsInLowPk,2));
unitsInLowPk = unitsInLowPk(:,1:(end-1));