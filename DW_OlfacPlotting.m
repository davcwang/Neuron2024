%Open eventsNPX and StNPX .mat file
Expt = 'DW-XI-12-1-PPC';

% Get expt spec info
[TRAPcids, TRAPidx, age, ExptNo, N, OID, OO, OOforA, GOVA] = DW_GetTRAPcidsidx(Expt,NPXSpikes); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sync APC and PPC probes
% 1) Open PPC NPXSpikes and rename NPXSpikes1
% 2) Open APC NPXSpikes and rename NPXSpikes2
% 3) Run the following
% 4) If using concat prepro .mat's, first load the eventsNPX from each .mat,
% then manually make a 2-column eventsNPX, then manually run Block 2 in DW_SyncProbes

DW_SyncProbes(FileName,MainDir,myKsDir,NPXSpikes1,NPXSpikes2);

% This outputs an NPXSpikesSync with concat cids (cid 1xxxx is the xxxxx-th cid from probe 1, 2xxxx for probe 2; idx's are just concat'ed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NPX_MasterScript for curated units

% %Generates struct from curated data; requires curating
NPXSpikes = loadKSdirGoodUnits(myKsDir); %Sorted Units, generates data struct
NPXSpikes.SpikeTimes = NPX_GetBeastCompatSpikeTimes(NPXSpikes);

%Some more parameters (this is a repeat of the above for uncurated units)
Probe = 1; 
CeventsMCCt = CeventsMCC(:,Probe); %TTLs, usually final valve opening and closings.
PREXmatFV = CeventsMCCt(StMCC == -1); 

%Generate Raster params and input for plotting odor-response PSTHs
[PREXOdorTimes,Odors] = NPX_PREX2Odor(PREXmatFV,OlfacMat,1);
NPXSpikes.ValveTimes = NPX_GetBeastCompatValveTimes(PREXOdorTimes);
Raster = NPX_RasterAlign(NPXSpikes.ValveTimes,NPXSpikes.SpikeTimes); %cell array of spikes aligned to the 1st param.

%Get PSTHstruct to look at average cell-odor pair responses
plotparams.PSTH.Axes = 'on';
plotparams.PSTHparams.Axes = 'on';
plotparams.OnlyData = true;
plotparams.VOI = 1:14; % default is 1:length(Odors)
plotparams.PSTHparams.PST = [-1,2]; % default is [-1,2] (-1s pre odor, and 2s post odor onset)
plotparams.PSTHparams.KernelSize = 0.02; % 20ms bins

[ltd_pre,td_pre] = NPX_GetTD(Raster,[-temprange,0],binsize);
[ltd_post,td_post] = NPX_GetTD(Raster,[0,temprange],binsize);
spikes_pre = td_pre(:,:);
spikes_post = td_post(:,:);

%Generate final PSTHs for saving to PDF
NPX_GetPSTHpdf(Raster,NPXSpikes.SpikeTimes,'Odor'); %Plot PSTHs
NPX_PlotAveragePSTH(PSTHstruct,valvesToPlot,idxToPlot,1,[],1); %averaged PSTH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create odor output data
DW_OlfacFigures(Expt,NPXSpikes,spikes_pre,spikes_post,OO,OOforA,Raster,TRAPidx,age,ExptNo)