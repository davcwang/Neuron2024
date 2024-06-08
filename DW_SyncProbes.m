function DW_SyncProbes(FileName,MainDir,myKsDir,NPXSpikes1,NPXSpikes2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Block 1: Generate eventsNPX matrix

% First do this for bank 1?
Banks = 0;

%Set directories
RecordingDir = [MainDir,'\experiment1\recording1'];
EventsDir = [RecordingDir,'\events'];

%Extract folder names in events directory
cd(EventsDir)
fnames = dir;
fnames = {fnames.name};

eventsNPX = []; %binary events

% Loop through all folders
for ii = 1:length(fnames)
           
    if contains(fnames{ii},'Neuropix')
        
        % Navigate to Neuropix\TTL folder
        cd([EventsDir,'\',fnames{ii}])
        tempF = dir;
        tempF = {tempF.name};
        cd([EventsDir,'\',fnames{ii},'\',tempF{end}])
        
        % Read in NPX timestamps as double
        tempeve = double(readNPY('timestamps.npy'));
        tempbool = tempeve>1; % crop off any non-sample time elements
        eventsNPX = [eventsNPX,tempeve(tempbool)]; % append NPX real sample times to the empty eventsNPX array

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Block 2: Sync probes and output synced file

cd(myKsDir)

% Create copies of NPXSpikes structs
NPXSpikes1Sync = NPXSpikes1;
NPXSpikes2Sync = NPXSpikes2;

% Function for syncing spike  times between two probes (this isn't perfect, it will be off by +/- 1 sample... but that's fine
eventsNPXdesync = eventsNPX(:,1) - eventsNPX(:,2); % array of temporal drift between two probes, for the NPX TTLs

c = polyfit(eventsNPX(:,1), eventsNPXdesync, 1); % fit a 1st order polynomial for the temporal drift vs event time (found in eventsNPX(:,1))
% c returns [slope, y-int]

syncSpikes1 = NPXSpikes1.ss; % array for storing synchronized spikes from probe 1 (this will just be the regular probe 1 sample times, since probe 1 is the reference)
syncSpikes2 = uint64(zeros(length(NPXSpikes2.ss),1)); % array for storing synchronized spikes from probe 2; cast as uint64 bc that's the original class

for i = 1:length(NPXSpikes2.ss) % loop through all the samples in 2nd probe

    spike = double(NPXSpikes2.ss(i)); % the reference (ie. independent variable) sample time is the current event/spike in the 2nd probe; cast as double so we can do algebra on it
    drift = round(c(1)*spike + c(2)); % calculate the amount of drift based the mx+b (where m = c(1), x = sample, b = c(2))
    
    syncSpikes2(i) = uint64(spike + drift); % store correct sample in syncSpikes2; recast as uint64 to store as original class

end

% Set NPXSpikesSync to new times
NPXSpikes2Sync.ss = syncSpikes2;
NPXSpikes2Sync.st = double(syncSpikes2)/30000;

% Create new cids (cids in NPXSpikes1 become 1000+cid, in NPXSpikes2 become 2000+cid)
NPXSpikes1Sync.cids = NPXSpikes1.cids + 10000;
NPXSpikes2Sync.cids = NPXSpikes2.cids + 20000;
NPXSpikes1Sync.clu = NPXSpikes1.clu + 10000;
NPXSpikes2Sync.clu = NPXSpikes2.clu + 20000;

% Concat info about cids into new struct NPXSpikesSync
NPXSpikesSync = NPXSpikes1; % create a new struct to store everything in
NPXSpikesSync.cids = [NPXSpikes1Sync.cids NPXSpikes2Sync.cids]; % idx of cids2 will just be tacked on after idx of cids1; same for all below
NPXSpikesSync.cgs = [NPXSpikes1Sync.cgs NPXSpikes2Sync.cgs];
NPXSpikesSync.CluDepth = [NPXSpikes1Sync.CluDepth NPXSpikes2Sync.CluDepth];
NPXSpikesSync.CluFreq = [NPXSpikes1Sync.CluFreq NPXSpikes2Sync.CluFreq];
NPXSpikesSync.CluAmp = [NPXSpikes1Sync.CluAmp NPXSpikes2Sync.CluAmp];

% Combine and sort ss and clu
tempss = [NPXSpikes1Sync.ss' NPXSpikes2Sync.ss']; % concat ss
tempclu = [NPXSpikes1Sync.clu' NPXSpikes2Sync.clu']; % concat clu
tempssclu = [tempss ; tempclu]; % make matrix with ss and clu
tempssclut = [tempss ; tempclu]'; % transpose for sortrows
tempsorted = sortrows(tempssclut,1)'; % sort based on 1st column (ss), retranspose into rows

NPXSpikesSync.ss = tempsorted(1,:)'; % assign the sorted ss
NPXSpikesSync.clu = tempsorted(2,:)'; % assign the sorted clu

cd(MainDir)

save('NPXSpikesSync','NPXSpikesSync');

% Things that weren't combined: spikeTemplates, tempScalingAmps, xcoods,
% ycoods, temps, winv, pcFeat, pcFeatInd, CluContamPct, SpikeTimes,
% ValveTimes