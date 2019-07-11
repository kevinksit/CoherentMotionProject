%retinotopicDependency_aggregateAnalysis;

% cleaned up version of RF_CM_aggregateAnalysis.m

% Written 10Jul2019 KS

%% Analysis parameters
anova_threshold = 0.01; % Threshold for determining of the cells are significantly tuned, according to ANOVA;

shuffle_iterations = 1000; % Iterations for shuffling the correlation between preference and location on imaging field
shuffle_threshold = 99; % percentile of the distribution to cut off to be considered well tuned

%% Get your files
CM_files = dir('**/coherenceResponseData.mat');
RF_files = dir('*/*RFmapping_results.mat');

if ~(length(CM_files) == length(RF_files))
    fprintf('You have a serious issue.. your recordings are mismatched...\n')
    return
end


%% initialize cell arrays
azi_pref_cell              = {};
alt_pref_cell              = {};
isSpatiallyTuned_cell      = {};

CC_mean_coherence_cell     = {};
isCoherenceResponsive_cell = {};
pref_dir_cell              = {};
reliability_cell           = {};

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
    isANOVA{f}                 = azi_p < anova_threshold & alt_p < anova_threshold;
    
    % shuffle analyses
    x_locations                = roi_centroids(isSpatiallyTuned,1);
    y_locations                = roi_centroids(isSpatiallyTuned,2);

    for iter = 1:shuffle_iterations % for now, maybe higher later?
        shuff_x = x_locations(datasample(1:length(x_locations),length(x_locations),'Replace',false));
        shuff_y = y_locations(datasample(1:length(y_locations),length(y_locations),'Replace',false));
        alt_shuff(iter) = corr(alt_pref(isSpatiallyTuned)',shuff_x);
        azi_shuff(iter) = corr(azi_pref(isSpatiallyTuned)',shuff_y);
    end
    
    isGoodAzi(f)               = azi_corr > prctile(azi_shuff,shuffle_threshold,2);
    isGoodAlt(f)               = alt_corr > prctile(alt_shuff,shuffle_threshold,2);
end

% Some further processing and turning from cell arrays to the corresponding vectors

%start here