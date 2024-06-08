function [odorParamsMat odorParamsMatP] = DW_OdorResponseParams(Expt, NPXSpikes, responseIndexMatP, OO, OID, Raster, N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting odor response parameters (peak rate, time to peak, peak width)
% 1) Open reseponIndexMatP, NPXSpikes, .mat file
% 2) Run first block below to get KDFStruct and PSTHStruct

% Odor IDs
OID = ["Et ac"; "Iso ac"; "Hex"; "Et bu"; "Maternal ur"; "Male ur"; "MO"; "VA ventrum"; "GO ventrum"; "Male Ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"; "Unmask 2-hex"; "Unmask VA vent"; "Unmask Male ur"; "Post MO"; "Post iso ac"; "Post et bu"; "Post Maternal ur"; "Post Male ur"];
odorsToPlot = OO; % inherit order of odors from input

% Get KDFstruct to look at cell-odor response parameters
KDFRaster = Raster;
KDFparams.PST = [-1,2];
KDFparams.kernelsize = (0.05); % 50ms kernel
KDFparams.trials = 1:10; % 10 trials per odor
[KDFStruct.KDF, KDFStruct.KDFtrials, KDFStruct.KDFt, KDFStruct.KDFe] = KDFmaker_Beast(KDFRaster,KDFparams.PST,KDFparams.kernelsize,KDFparams.trials);

%Get PSTHstruct to look at average cell-odor pair responses
plotparams.PSTH.Axes = 'on';
plotparams.PSTHparams.Axes = 'on';
plotparams.OnlyData = true;
plotparams.VOI = 1:size(Raster,1); % this is usually 1:length(Odors), but we aren't inputting Odors from the .mat file, so just get this from Raster length
plotparams.PSTHparams.PST = [-1,2]; % default is [-1,2] (-1s pre odor, and 2s post odor onset)
plotparams.PSTHparams.KernelSize = 0.05; % 50ms bins
PSTHstruct = NPX_RasterPSTHPlotter(Raster,NPXSpikes.SpikeTimes,plotparams);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After running first block to get KDFStruct and PSTHStruct, run this second block to get odorParamsMatP

odorParamsMat = zeros(14,3,length(NPXSpikes.cids)); % odors x params ([peak rate, time to peak, 50% width]) x units
odorParamsMatP = zeros(14,3,length(NPXSpikes.cids)); % odors x params ([peak rate, time to peak, 50% width]) x units
binsize = mean(diff(KDFStruct.KDFt)); % get binsize as average difference in timesteps in KDF

for u = 1:length(NPXSpikes.cids) % loop over all units u

    for v = 1:size(Raster,1) % loop over all odors (where v = valve)

        odorInd = find(odorsToPlot==v); % indexes odor (ie. when v = 1 [first sigma odor], odorInd = 2, so we're looking at the first sigma odor in ZScoreMatP and saving this to the 2nd column in odorParamsMatP

        [pk, loc, width] = findpeaks(KDFStruct.KDF{v,:,u}); % find peaks from chronux KDF

        pk(pk < 0) = []; % remove pks before odor onset

        if length(pk) > 1 % if multiple pks, find global pk and bin number
            loc = loc(find(pk == max(pk))); % find loc of global max
            width = width(find(pk == max(pk)));
        end

        % Calculate peak rate and time to peak from pk and loc
        peakRate = PSTHstruct.KDF{v,u}(loc); % get peak rate at loc from PSTH of v-th odor and u-th unit
        timeToPeak = KDFStruct.KDFt(loc); % get time of peak from KDFt
        FWHM = width * binsize;

        % Code for finding std of firing rate to see if peak rate is above some std above mean pre-odor firing rate
        preFiringt = find(PSTHstruct.KDFt < 0,1,'last');
        preFiringRate = mean(PSTHstruct.KDF{v,u}(1:preFiringt));
        preFiringSTD = std(PSTHstruct.KDF{v,u}(1:preFiringt));

        % Only save params if there's activation
        % if isempty(pk) || responseIndexMatP(u,odorInd) < 0.4 || peakRate < (preFiringRate + 3*preFiringSTD) % if there are no peaks, r <0.4, or peak rate is greater than 3*stdev from baseline firing rate
        if isempty(pk) || responseIndexMatP(u,v) < 0.4 || peakRate < (preFiringRate + 3*preFiringSTD) % if there are no peaks, r < 0.4, or peak rate is greater than 3*stdev from baseline firing rate
            % do nothing
        else % save params
    
            odorParamsMat(v,:,u) = [peakRate timeToPeak FWHM]; % populate odorParamsMat with params for visually checking with odor PSTH .pdfs
            % odorParamsMatP(odorInd,:,u) = [peakRate timeToPeak FWHM]; % populate odorParamsMatp with params
        
        end
    end
end

