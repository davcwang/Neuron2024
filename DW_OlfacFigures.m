function DW_OlfacFigures(Expt,NPXSpikes,spikes_pre,spikes_post,OO,OOforA,RasterforA,TRAPidx,age,ExptNo)

T = 10; % number of trials
N = 14;
OID = ["MO"; "Iso ac"; "Et ac"; "Hex"; "Et bu"; "Familiar ventrum"; "Novel ventrum"; "Male Ventrum"; "Familiar milk"; "Novel milk"; "Familiar AF"; "Novel AF"; "Maternal ur"; "Male ur"];
OdorTicks = OID; % OID is already in the right order

% Generate ZScoreMat
ZScoreMat = zeros(N,length(spikes_pre(1,:))); % Generate z-score matrix for odors
spikeMat = zeros(N,length(spikes_pre(1,:))); % Generate spike matrix for odors (for LFS)
for i = 1:N %loop over odors    
    for j = 1:length(spikes_pre(1,:)) %loop over units

        %Get total spike counts for this odor across all trials
        spikeSum = sum(spikes_post(T*(i-1)+1:T*i,j)');
        spikeMat(i,j) = spikeSum;
        
        %Calculate Z-score relative to one time bin of baseline pre-odor
        muVec = spikes_pre(T*(i-1)+1:T*i,j)';
        postVec = spikes_post(T*(i-1)+1:T*i,j)';
       
        post = mean(postVec);
        mu = mean(muVec);
        sigma = std(muVec);
        
        z = (post - mu) / sigma;

        % If sigma = 0 for this cell-odor pair (thus returning z = NaN or Inf), recalculate sigma using the baseline firing across all odors
        if sigma == 0
            sigma = std(spikes_pre(:,j)); % if sigma = 0 (ie. std of this unit's baseline firing for this odor is 0), set sigma = std of this unit for all odors
            z = (post - mu) / sigma; % recalculate z with new sigma
            if sigma == 0 % if sigma is still 0 (ie. this unit had no spikes or the exact same number of spikes at baseline)
                z = 0;
            end
        end
        
        ZScoreMat(i,j) = z;
    end
end

%Plot z-scored to baseline heatmap
figure(1)
ZScoreMatP = ZScoreMat';
imagesc(ZScoreMatP);
colormap(colormap_BlueWhiteRed);
caxis([-std(ZScoreMatP(:)) std(ZScoreMatP(:))]) %plot the caxis as 2*std of all elements in the matrix, allowing better visualization of the medium values
set(gca,'xtick',[1:N],'xticklabel',OdorTicks);
ylabel('Unit Number');
title('Z-score over baseline')
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rank sum for odor response vs baseline within cells, for sparseness matrix

if contains(Expt,'13') % if adult, do 1s time windows
    MCSR_PST = [0,1];
    MCSR_PSTpre = [-1,0];
else % if neonate, do 2s time windows
    MCSR_PST = [0,2];
    MCSR_PSTpre = [-2,0];
end
NPXMCSR = NPX_GetMultiCycleSpikeRate(RasterforA,1:10,MCSR_PST); % 1:10 for 10 trials
NPXMCSRPr = NPX_GetMultiCycleSpikeRate(RasterforA,1:10,MCSR_PSTpre); % for pre-odor baseline
OdorScores = NPX_SCOmakerPreInh(NPXMCSR,NPXMCSRPr);

% Calculate and plot response index (-1 is unequivocal inh, +1 is exc)
responseIndexMat = (OdorScores.auROC*2 - 1)'; % calculate response index as 2 * auROC - 1 

figure(2)
responseIndexMatP = responseIndexMat; % P stands for Perm
imagesc(responseIndexMatP);
colormap(colormap_BlueWhiteRed);
caxis([-1 1]) % always btwn -1 and 1 bc those are bounds of r
set(gca,'xtick',[1:N],'xticklabel',OdorTicks);
ylabel('Unit Number');
title('Response index r')
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get odor parameters (peak latency, peak rate, peak width)

[odorParamsMat,odorParamsMatP] = DW_OdorResponseParams(Expt, NPXSpikes, responseIndexMatP, OO, OID, RasterforA);

% Collapse odor response across odors to avg and std to generally characterize how a cell responds to odors
avgOdorParams = zeros(size(odorParamsMat,3),3); 
stdOdorParams = zeros(size(odorParamsMat,3),3);

% Reshape odorParamsMatP to be a units x odors*3 matrix for storage in matrix A, so each row is a unit and each triplet of columns is the 3 odor params for an odor
odorParamsReshape = zeros(size(odorParamsMat,3),14*3);

% Reshape odorParamsMatP and get avg/std of odor response
for u = 1:size(odorParamsMat,3)
    for o = 1:size(odorParamsMat,1)
        odorParamsReshape(u, (3*(o-1) + 1) : 3*o ) = odorParamsMat(o,:,u); % add the 3 params into concatted rows of odorParamsReshape
    end

    tempParamsMat = odorParamsMat(:,:,u); % temp matrix for odor params of u-th unit
    tempParamsMat = tempParamsMat(~all(tempParamsMat == 0, 2),:); % remove empty rows using "logical indexing"

    if isempty(tempParamsMat) % if there are no odor responses, set to avg and std 0
        avgOdorParams(u,:) = 0; % take avg
        stdOdorParams(u,:) = 0; % take std
    else
        avgOdorParams(u,:) = mean(tempParamsMat,1); % take avg
        stdOdorParams(u,:) = std(tempParamsMat,1); % take std
    end
end

figure(4)
imagesc(odorParamsReshape);
colormap(colormap_BlueWhiteRed);
caxis([0 2]) % always btwn -1 and 1 bc those are bounds of r
ylabel('Unit Number');
title('Odor Params Reshape')
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lifetime sparseness

SrMat = spikeMat';
emptyCols = find(OOforA == 0);
SrMat(:,emptyCols) = [];

N = length(SrMat(1,:)); %number of odors
C = length(SrMat(:,1)); %number of cells

Sr = zeros(C,1); %test response matrix

for i = 1:C
    
    term1 = 1 / (1 - 1/N);
    num = ( sum(SrMat(i,:)) / N )^2;
    denom = sum(SrMat(i,:).^2) / N;
    term2 = ( 1 - num  / denom );

    Sr(i) = term1*term2;

end

figure(3)
hist(Sr)
ylim([0 C])
ylabel('Number of cells')
xlim([0 1])
xlabel('Lifetime sparseness')
title('Lifetime sparseness')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Population sparseness to test if odor intensity is decently matched

% Get rid of empty columns when the specific odors weren't included in this expt
SpMat = spikeMat';
emptyCols = find(OOforA == 0);
SpMat(:,emptyCols) = [];
SpMat = SpMat'; % tranpose so we perform analysis wrt odors instead of units

N = length(SpMat(1,:)); %number of odors
C = length(SpMat(:,1)); %number of cells

Sp = zeros(C,1); %test response matrix

for i = 1:C
    
    term1 = 1 / (1 - 1/N);
    num = ( sum(SpMat(i,:)) / N )^2;
    denom = sum(SpMat(i,:).^2) / N;
    term2 = ( 1 - num  / denom );

    Sp(i) = term1*term2;

end

figure(4)
histogram(Sp)
%ylim([0 N])
ylabel('Number of cells')
xlim([0 1])
xlabel('Population sparseness')
title('Population sparseness')