function [unitPksCell] = DW_GetUnitPks(NPXSpikes, ZConvSorted, smoothBin, unitPksCell)

% Find pks for each unit to compare to pop peaks
%for i = 1:length(NPXSpikes.cids) % loop over all units to generate cell array with each unit's pks, locs, widths
for i = 1:length(ZConvSorted(:,1)) % loop over all units to generate cell array with each unit's pks, locs, widths

    smoothUnitSpikes = smooth(ZConvSorted(i,:),smoothBin,"lowess"); % smooth unit spikes in same way as smoothZpopAct
    stdUnitAct = std(smoothUnitSpikes); % for calculating peaks in this unit's activity

    [pks,loc,width] = findpeaks(smoothUnitSpikes); % find pks of i-th unit (i = idx)

    highPks = find(pks > 1 * std(stdUnitAct) + median(smoothUnitSpikes)); % find indices of pks > median + std
    unitPks = pks(highPks);
    unitLocs = loc(highPks); % find locations (in bin number) of real peaks of this unit 
    unitWidths = width(highPks); % find widths are real peaks of this unit

    % Store in cell array (row = unit idx, column = 1 is pks, 2 is locs, 3 is widths)
    unitPksCell{i,1} = unitPks;
    unitPksCell{i,2} = unitLocs;
    unitPksCell{i,3} = unitWidths;

end