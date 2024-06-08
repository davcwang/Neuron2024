function [popPks, popLocs, popWidths] = DW_GetPopPks(smoothZpopAct,threshold)

[pks,loc,width] = findpeaks(smoothZpopAct); % find pk values, locations, FHWM

thresPks = find(pks > threshold); % find indices of high pks > std above median

popPks = pks(thresPks);
popLocs = loc(thresPks); % find locations (in bin number) of real peaks
popWidths = width(thresPks); % find widths are real peaks
