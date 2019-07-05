CM_files = dir('**/coherenceResponseData.mat');
RF_files = dir('*/*RFmapping_results.mat');

ct=1;
for f = 1:length(RF_files) % the idea is that we can first put everything into big cell arrays, letting us do individual recording analyses..
    load([RF_files(f).folder '\' RF_files(f).name])
    load([CM_files(f).folder '\' CM_files(f).name])
 
%% aggregate the data
    azi_all_cell{f} = azi_pref;
    alt_all_cell{f} = alt_pref;

    for c = 1:size(azi_pref,2)
        temp_fit = azi_fit{c}(1:40);
        temp_fit(temp_fit<0.01) = 0;
        azitempfits(c,:) = temp_fit;
    end
    
    azi_fits_cell{f} = [azitempfits];
    
    for c = 1:size(alt_pref,2)
        temp_fit2 = alt_fit{c}(1:30);
        temp_fit2(temp_fit2<0.01) = 0;
        alttempfits(c,:) = temp_fit2;
    end
    
    alt_fits_cell{f} = [alttempfits];
    
    azitempfits = [];
    alttempfits = [];
    
    CC_mean_coherence_cell{f} = ([CC_mean_coherence]);
    isCohCC_cell{f} = [isCoherenceResponsive];
    
    isSpat_cell{f} = [isSpatiallyTuned];
    
    azi_allp_cell{f} = [azi_p];
    alt_allp_cell{f} = [alt_p];
    
    azi_corr_aggregate(f) = azi_corr;
    alt_corr_aggregate(f) = alt_corr;
    
    
    % calculate some stuff
    isANOVA_temp = azi_allp_cell{f}<0.01 | alt_allp_cell{f}<0.01;
    isSpat = isSpatiallyTuned & isANOVA_temp;
    temp = corrcoef(alt_all_cell{f}(isSpat),CC_mean_coherence_cell{f}(isSpat));
    altCC_session(f) = temp(1,2);
    temp = corrcoef(azi_all_cell{f}(isSpat),CC_mean_coherence_cell{f}(isSpat));
    aziCC_session(f) = temp(1,2);
    
    azi_mean(f) = mean(azi_all_cell{f}(isSpat));
    alt_mean(f) = mean(alt_all_cell{f}(isSpat));
    
    frac_included(f) = mean(isSpat);
    numCells(f) = sum(isSpat);
    
    % calculating DSI
    RDK_resp = mean(mean(RespVec,3),4);
    RDK_resp = RDK_resp - min(RDK_resp(:)); %rectify
    [~,pref_dir] = max(RDK_resp,[],2);
    anti_dir = wrapN(pref_dir+size(RDK_resp,2)/2,size(RDK_resp,2));
    
    for c = 1:size(RDK_resp,1)
      DSI_cell{f}(c) = (RDK_resp(c,pref_dir(c)) - RDK_resp(c,anti_dir(c))) ./ (RDK_resp(c,pref_dir(c)) + RDK_resp(c,anti_dir(c)));
    end
    
    % shuffle analysisssss
    x_pos_cell{f} = roi_centroids(:,1);
    y_pos_cell{f} = roi_centroids(:,2);
    
    x_locations = roi_centroids(isSpat,1);
    y_locations = roi_centroids(isSpat,2);
    
    for iter = 1:1000 % for now, maybe higher later?
    shuff_x = x_locations(datasample(1:length(x_locations),length(x_locations),'Replace',false));
    shuff_y = y_locations(datasample(1:length(y_locations),length(y_locations),'Replace',false));
    alt_shuff(f,iter) = corr(alt_pref(isSpat)',shuff_x);
    azi_shuff(f,iter) = corr(azi_pref(isSpat)',shuff_y);
    end
    
    
    %% new stuff as of 05Jul2019
    pref_dir_cell{f} = pref_dir;
    reliability_cell{f} = reliability;
end

% %% weird normalization thing
% gm = mean([CC_mean_coherence_cell{:}]);
% 
% temp = cellfun(@(x) (x-mean(x))+gm,CC_mean_coherence_cell,'UniformOutput',false);
% CC_mean_coherence_cell = temp;

%% 

% rejecting recordings whose RFs aren't on the field
alt_thresh = prctile(alt_shuff,99,2);
azi_thresh = prctile(azi_shuff,99,2);

isGoodAlt = (alt_corr_aggregate > alt_thresh');% & (alt_corr_aggregate>0.5);
isGoodAzi = (azi_corr_aggregate > azi_thresh');% & (azi_corr_aggregate>0.5);

isGoodRecording = isGoodAlt;% & isGoodAzi;

%% turning from session to altogether (all cells analysis...)
azi_pref = cat(2,azi_all_cell{isGoodRecording});
alt_pref = cat(2,alt_all_cell{isGoodRecording});

y_pos = cat(1,x_pos_cell{isGoodRecording})';
isSpatiallyTuned = cat(2,isSpat_cell{isGoodRecording});
CC_mean_coherence = cat(2,CC_mean_coherence_cell{isGoodRecording});
azi_p = cat(2,azi_allp_cell{isGoodRecording});
alt_p = cat(2,alt_allp_cell{isGoodRecording});

DSI = cat(2,DSI_cell{isGoodRecording});

reliability = cat(1,reliability_cell{isGoodRecording});

% also ones that pass anova...
isANOVA = azi_p<0.01 | alt_p<0.01;

isReliable = (reliability>0.1)';

% spearman's?
corr(azi_pref(isSpatiallyTuned&isANOVA&isReliable)',CC_mean_coherence(isSpatiallyTuned&isANOVA&isReliable)','Type','Pearson')
corr(alt_pref(isSpatiallyTuned&isANOVA&isReliable)',CC_mean_coherence(isSpatiallyTuned&isANOVA&isReliable)','Type','Pearson')

corr(alt_pref(isSpatiallyTuned&isReliable)',CC_mean_coherence(isSpatiallyTuned&isReliable)','Type','Pearson')


% visualize
scatter(alt_pref(isSpatiallyTuned&isReliable)',CC_mean_coherence(isSpatiallyTuned&isReliable)')
scatplot(double(alt_pref(isSpatiallyTuned&isANOVA&isReliable'))',double(CC_mean_coherence(isSpatiallyTuned&isANOVA&isReliable'))')



%% Weighted representation?
azi_locs = 1:40;
alt_locs = 1:30;

alt_fits_all = cat(1,alt_fits_cell{isGoodRecording});
azi_fits_all = cat(1,azi_fits_cell{isGoodRecording});

alt_fits = alt_fits_all(isSpatiallyTuned&isANOVA&isReliable,:);
azi_fits = azi_fits_all(isSpatiallyTuned&isANOVA&isReliable,:);
CC_sig   = CC_mean_coherence(isSpatiallyTuned&isANOVA&isReliable);

for c = 1:length(alt_fits)
    alt_fitxCC(c,:) = rescale(alt_fits(c,:)) * CC_sig(c);
    weights_alt(c,:) = rescale(alt_fits(c,:));
end

for c = 1:length(azi_fits)
    azi_fitxCC(c,:) = rescale(azi_fits(c,:)) * CC_sig(c);
    weights_azi(c,:) = rescale(azi_fits(c,:));
end

% 
% for c = 1:length(alt_fits_all)
%     alt_fitxCC(c,:) = rescale(alt_fits_all(c,:)) * CC_mean_coherence(c);
%     weights_alt(c,:) = rescale(alt_fits_all(c,:));
% end
% 
% for c = 1:length(azi_fits_all)
%     azi_fitxCC(c,:) = rescale(azi_fits_all(c,:)) * CC_mean_coherence(c);
%     weights_azi(c,:) = rescale(azi_fits_all(c,:));
% end

alt_weightedCC = sum(alt_fitxCC)./sum(weights_alt);
azi_weightedCC = sum(azi_fitxCC)./sum(weights_azi);

alt_mat = repmat(alt_weightedCC,40,1)';
azi_mat = repmat(azi_weightedCC,30,1);

screen_mat = alt_mat+azi_mat;

% also let's try binned stuff, both fraction and the value w/in binz

bin_sz=3;
azi_groups = discretize(azi_pref,[0:bin_sz:length(azi_locs),length(azi_locs)]); % in case our bin size doesn't reach the length....
alt_groups = discretize(alt_pref,[0:bin_sz:length(alt_locs),length(alt_locs)]);

isIncluded = isANOVA & isSpatiallyTuned & isReliable;
for ii = 1:length(unique(azi_groups))
    isCurrBin = azi_groups == ii;
    azi_binned{ii} = CC_mean_coherence(isCurrBin&isIncluded);
end

for ii = 1:length(unique(alt_groups))
    isCurrBin = alt_groups == ii;
    alt_binned{ii} = CC_mean_coherence 

 errorbar(1:bin_sz:30,cellfun(@mean,alt_binned),cellfun(@(x) std(x)/sqrt(size(x,2)),alt_binned),'o')
 hold on
errorbar(1:bin_sz:40,cellfun(@mean,azi_binned),cellfun(@(x) std(x)/sqrt(size(x,2)),azi_binned),'o')

  

isLeftHalf = alt_pref<15;

isRightHalf = ~isLeftHalf;
    
mean(CC_mean_coherence(isRightHalf&isIncluded))



dat = alt_binned;
temp = cat(2,dat{:});
groups = [];

for ii = 1:length(dat)
    g_temp = ii * ones(1,length(dat{ii}));
    groups = [groups g_temp];
end




    

