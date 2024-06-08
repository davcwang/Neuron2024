%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for plotting optotag summary (can be adaptable to any event trigger; useful for plotting responses per cell type)

% Get expt data and TRAPidx
Expt = 'DW-XI-12-1-PPC';
[TRAPcids, TRAPidx, ~, ~, ~, ~, ~, ~, ~] = DW_GetTRAPcidsidx(Expt,NPXSpikes); 

% Define time window parameters
binSize = 1; % 0.002 = 2 ms bin size
windowSize = 2; % 0.08 = 80 ms window centered around each event time

% Initiate matrix (each row is cell's average 
optotagMat = zeros(length(NPXSpikes.cids),windowSize/binSize);

% Get opto event times, OptoProt2NPX to bin indices
if contains(Expt,'XI-8')
    optoTimes = [preOptoEventsNPX' postOptoEventsNPX']'/30000;
else
    optoTimes = OptoProt1NPX / 30000;
end
eventBins = round(optoTimes / binSize);

% Loop through all units
for i = 1:length(NPXSpikes.cids)

    % Convert spikeTimes to bin indices
    spikeTimes = NPXSpikes.SpikeTimes.tsec{i};
    spikeBins = round(spikeTimes / binSize);
    
    % Initialize binary matrix
    numBins = round(windowSize / binSize);
    binaryMatrix = zeros(length(optoTimes), numBins);
    
    % Populate binary matrix of firing rates over each bin
    for j = 1:length(optoTimes) % loop over j opto events
        centerBin = eventBins(j);
        startTime = centerBin - round(numBins / 2);
        endTime = centerBin + round(numBins / 2) - 1;
        
        % Check spike times within the window
        spikeIndices = find(spikeBins >= startTime & spikeBins <= endTime);
        
        % Count the number of spikes in each bin
        relativeIndices = spikeBins(spikeIndices) - startTime + 1;
        for k = 1:length(relativeIndices)
            binaryMatrix(j,relativeIndices(k)) = binaryMatrix(j,relativeIndices(k)) + 1;
        end
    
    end
    
    % Get variables for z score
    meanRate = mean(binaryMatrix,'all');
    stdDevRate = std(binaryMatrix(:));
    if sum(sum(binaryMatrix)) == 0
        zScore = zeros(1,length(windowSize/binSize));
    else
        zScore = (mean(binaryMatrix,1) - meanRate) / stdDevRate;
    end
    
    % Smooth
    smoothZ = smooth(zScore,3);
    
    % Populate optotagMat
    optotagMat(i,:) = smoothZ;
end

% Plot
imagesc(optotagMat);
colormap(colormap_BlueWhiteRed);
caxis([-1 1]) %plot the caxis as 2*std of all elements in the matrix, allowing better visualization of the medium values
set(gca,'xtick',[1:40],'xticklabel',[-40:2:40]);
ylabel('Unit Number');
title('Z-score optotag')
colorbar;

% Get TRAPed and nonTRAPed
TRAPoptotagMat = optotagMat(TRAPidx,:);
nonTRAPoptotagMat = optotagMat;
nonTRAPoptotagMat(TRAPidx,:) = 0;