function [spontMat, unitsInPk] = DW_GetUnitCoPks(NPXSpikes,spontMat,unitPksCell,popPks,popLocs,popWidths,unitsInPk, distCol, coPkCol)

% Loop through all pop peaks to see if unit pks co-occur, how much they differ in pk time, and how many units co-peak within each peak
for i = 1:length(popPks)
    
    popPkLoc = popLocs(i); % find loc of current pop pk
    popPkWidth = popWidths(i); % find width at current pop peak

    % for u = 1:length(NPXSpikes.cids) % loop through all units to look at their individual pks
    for u = 1:length(unitPksCell(:,1)) % loop through all units to look at their individual pks
        
        unitLocs = unitPksCell{u,2}; % get locs of pks of current unit

        % Find ind of unit pks that are within 25% width of pop pk
        prePks = find(unitLocs > (popPkLoc - 0.5*popPkWidth)); 
        postPks = find(unitLocs < (popPkLoc + 0.5*popPkWidth));
        coPks = intersect(prePks,postPks);

        if length(coPks) == 1 % if there is a co peak,             
            spontMat(u,distCol) = spontMat(u,distCol) + (unitLocs(coPks) - popPkLoc); % get distance from unit's pk to pop pk (positive means unit pks after pop pk)
            unitsInPk(u,i) = unitsInPk(u,i) + 1; % increment count in units in pop peak
            spontMat(u,coPkCol) = spontMat(u,coPkCol) + 1; % increment count of pop pks that this unit has a peak in
        
        elseif length(coPks) > 1 % if there is more than one co pk
            for p = 1:length(coPks) % loop through all co pks and get distance for each copk              
                spontMat(u,distCol) = spontMat(u,distCol) + (unitLocs(coPks(p)) - popPkLoc); % get distance from unit's pk to pop pk (positive means unit pks after pop pk)           
            end
            unitsInPk(u,i) = unitsInPk(u,i) + 1; % increment count in units in pop peak
            spontMat(u,coPkCol) = spontMat(u,coPkCol) + 1; % increment count of pop pks that this unit has a peak in
        end
    
    end
end

