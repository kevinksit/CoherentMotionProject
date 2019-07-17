%retinotopicDependency_aggregateAnalysis;

% cleaned up version of RF_CM_aggregateAnalysis.m

% Written 10Jul2019 KS

%% Analysis parameters
anova_threshold = 0.05; % Threshold for determining of the cells are significantly tuned, according to ANOVA;

shuffle_iterations = 1000; % Iterations for shuffling the correlation between preference and location on imaging field
shuffle_threshold = 95; % percentile of the distribution to cut off to be considered well tuned

reliability_threshold = 0.2; % threshold for reliability to be included, just to get rid of crap cells
%% Get your files
CM_files = dir('**/coherenceResponseData.mat');
RF_files = dir('*/*RFmapping_results.mat');

if ~(length(CM_files) == length(RF_files))
    fprintf('You have a serious issue.. your recordings are mismatched... \n')
    return
end


%% initialize cell arrays
azi_pref_cell              = {};
alt_pref_cell              = {};
isSpatiallyTuned_cell      = {};
isANOVA_cell               = {};

CC_mean_coherence_cell     = {};
isCoherenceResponsive_cell = {};
pref_dir_cell              = {};
reliability_cell           = {};

isGoodAlt                  = logical.empty();
isGoodAzi                  = logical.empty();

%% Aggregate your data
for f = 1:length(RF_files) % going through the RF files...
    
    load([RF_files(f).folder '\' RF_files(f).name]) % load the appropriate data
    load([CM_files(f).folder '\' CM_files(f).name])
    
    % Retiontopic preferences
    azi_pref_cell              = cellArrayAggregator(azi_pref_cell,'azi_pref');
    alt_pref_cell              = cellArrayAggregator(alt_pref_cell,'alt_pref');
    isSpatiallyTuned_cell      = cellArrayAggregator(isSpatiallyTuned_cell,'isSpatiallyTuned');
    
    % coherence correlation
    CC_mean_coherence_cell     = cellArrayAggregator(CC_mean_coherence_cell,'CC_mean_coherence');
    isCoherenceResponsive_cell = cellArrayAggregator(isCoherenceResponsive_cell,'isCoherenceResponsive');
    pref_dir_cell              = cellArrayAggregator(pref_dir_cell,'pref_dir');
    reliability_cell           = cellArrayAggregator(reliability_cell,'reliability');
    
    % Calculations that need to be done within each recording...
    isANOVA_cell{f}            = alt_p < anova_threshold;
    
    % Shuffling for field significance
    x_locations                = roi_centroids(isSpatiallyTuned,1);
    y_locations                = roi_centroids(isSpatiallyTuned,2);

    for iter = 1:shuffle_iterations
        shuff_x = x_locations(datasample(1:length(x_locations),length(x_locations),'Replace',false));
        shuff_y = y_locations(datasample(1:length(y_locations),length(y_locations),'Replace',false));
        alt_shuff(iter) = corr(alt_pref(isSpatiallyTuned)',shuff_x);
        azi_shuff(iter) = corr(azi_pref(isSpatiallyTuned)',shuff_y);
    end
    
    isGoodAzi(f)               = azi_corr > prctile(azi_shuff,shuffle_threshold,2);
    isGoodAlt(f)               = alt_corr > prctile(alt_shuff,shuffle_threshold,2);
end

% Some further processing and turning from cell arrays to the corresponding vectors

isGoodRecording = true(1,27); % determines which recordings are good to include

% putting everything together

azi_pref = cat(2,azi_pref_cell{isGoodRecording});
alt_pref = cat(2,alt_pref_cell{isGoodRecording});
isSpatiallyTuned = cat(2,isSpatiallyTuned_cell{isGoodRecording});
isANOVA = cat(2,isANOVA_cell{isGoodRecording});

CC_mean_coherence = cat(1,CC_mean_coherence_cell{isGoodRecording});
isCoherenceResponsive = cat(1,isCoherenceResponsive_cell{isGoodRecording})';
reliability = cat(1,reliability_cell{isGoodRecording});

isReliable = (max(reliability,[],2) > reliability_threshold)';


%% Analyse
isIncluded = isSpatiallyTuned;

% Calculate the correlation
corr(azi_pref(isIncluded)',max(CC_mean_coherence(isIncluded,:),[],2),'Type','Pearson')
corr(alt_pref(isIncluded)',max(CC_mean_coherence(isIncluded,:),[],2),'Type','Pearson')

% Visualize? 

scatplot(alt_pref(isIncluded)',max(CC_mean_coherence(isIncluded,:),[],2)')


%% Binning stuff again... round 2

alt = alt_pref(isIncluded);
cc = max(CC_mean_coherence(isIncluded,:),[],2);

alt_bins = discretize(alt,0:3:30);

for ii = 1:length(unique(alt_bins))
    
    binned_resp(ii) = mean(cc(alt_bins == ii));
end

plot(binned_resp)

%% Direciton analysis
pref_dir = cat(1,pref_dir_cell{isGoodRecording});
pref_dir = pref_dir(isIncluded);

temp = tabulate(pref_dir);
theta = cat(1,0,temp(:,1))*(pi/4);
rho = cat(1,temp(:,2),temp(1,2));
p = polarplot(theta,rho);
pax = gca;
set(p,'Marker','o','MarkerSize',5,'LineStyle','--','LineWidth',2);
set(pax,'ThetaZeroLocation','top','ThetaDir','clockwise');
title('Directional histogram')

%% Quartile analysis...
alt = alt_pref(isIncluded);
cc = max(CC_mean_coherence(isIncluded,:),[],2);

alt_bins = discretize(alt,0:3:30);

for ii = 1:length(unique(alt_bins))
    
    quartile_resp(ii,:) = prctile(cc(alt_bins == ii),[10:10:100]);
end
    
% im imagining some kind of pseudocolored thing where the dots are labeled based on their percintile
for ii = 1:length(unique(alt_bins)) % num altitudes...
   temp = cc(alt_bins == ii);
   temp2 = alt(alt_bins == ii);
   for j = 1:size(quartile_resp,2)
       if j == 1
           color_resp{ii,j} = temp(temp < quartile_resp(ii,j));
           alt_resps{ii,j} = temp2(temp < quartile_resp(ii,j));
       else
   color_resp{ii,j} = temp(temp < quartile_resp(ii,j) & temp>quartile_resp(ii,j-1));
              alt_resps{ii,j} = temp2(temp < quartile_resp(ii,j) & temp>quartile_resp(ii,j-1));

   
       end
   end
end


colors = parula(10);


    figure;
    hold on
for ii = 1:size(color_resp,1)
    for jj = 1:size(color_resp,2)
        scatter(alt_resps{ii,jj},color_resp{ii,jj},'filled','MarkerFaceColor',colors(jj,:));
    end
end
