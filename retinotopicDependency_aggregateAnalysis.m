%retinotopicDependency_aggregateAnalysis;

% cleaned up version of RF_CM_aggregateAnalysis.m
% Written 10Jul2019 KS

%% Analysis parameters
anova_threshold = 0.05; % Threshold for determining of the cells are significantly tuned, according to ANOVA;

shuffle_iterations = 1; % Iterations for shuffling the correlation between preference and location on imaging field
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
pref_mag_cell              = {};
reliability_cell           = {};

isGoodAlt                  = logical.empty();
isGoodAzi                  = logical.empty();

azi_fits_cell              = {};
alt_fits_cell              = {};
azi_gof_cell               = {};
alt_gof_cell               = {};
%% Aggregate your data
for f = 1:length(RF_files) % going through the RF files...
    
    load([RF_files(f).folder '\' RF_files(f).name]) % load the appropriate data
    load([CM_files(f).folder '\' CM_files(f).name])
    
    % Retiontopic preferences
    azi_pref_cell              = cellArrayAggregator(azi_pref_cell,'azi_pref');
    alt_pref_cell              = cellArrayAggregator(alt_pref_cell,'alt_pref');
    isSpatiallyTuned_cell      = cellArrayAggregator(isSpatiallyTuned_cell,'isSpatiallyTuned');
    azi_fits_cell              = cellArrayAggregator(azi_fits_cell,'azi_fit');
    alt_fits_cell              = cellArrayAggregator(alt_fits_cell,'alt_fit');
    alt_gof_cell               = cellArrayAggregator(alt_gof_cell,'altgof');
    azi_gof_cell               = cellArrayAggregator(azi_gof_cell,'azigof');
    
    % coherence correlation
    CC_mean_coherence_cell     = cellArrayAggregator(CC_mean_coherence_cell,'CC_mean_coherence');
    isCoherenceResponsive_cell = cellArrayAggregator(isCoherenceResponsive_cell,'isCoherenceResponsive');
    pref_dir_cell              = cellArrayAggregator(pref_dir_cell,'pref_dirVS');
    pref_mag_cell              = cellArrayAggregator(pref_mag_cell,'pref_mag');
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
% % z score?
% narrow_thresh = 2; fit_thresh = 0.3;
% for ii = 1:length(isGoodRecording)
%     temp = CC_mean_coherence_cell{ii};
%     
%     for jj = 1:length(azi_fits_cell{ii})
%         azi_width(jj) = azi_fits_cell{ii}{jj}.c1;
%         alt_width(jj) = alt_fits_cell{ii}{jj}.c1;
%     end
%     
%     alt_rsquare = [alt_gof_cell{ii}.rsquare];
%     azi_rsquare = [azi_gof_cell{ii}.rsquare];
%     
%     isTooNarrow = (azi_width<narrow_thresh) | (alt_width<narrow_thresh);
%     isFit = (alt_rsquare>fit_thresh) | (azi_rsquare>fit_thresh);
%     isEdge = round(alt_pref_cell{ii}) == 30 | round(alt_pref_cell{ii}) == 1;
%     isSpatiallyTuned = isFit & ~isTooNarrow & ~isEdge;
%     
%     
%     temp = max(temp(isSpatiallyTuned,:),[],2);
%     CC_temp{ii} = zscore(temp);
%     clear azi_width alt_width
% end

% putting everything together


azi_pref = cat(2,azi_pref_cell{isGoodRecording});
alt_pref = cat(2,alt_pref_cell{isGoodRecording});
isSpatiallyTuned = cat(2,isSpatiallyTuned_cell{isGoodRecording});
isANOVA = cat(2,isANOVA_cell{isGoodRecording});
azi_fit = cat(2,azi_fits_cell{isGoodRecording});
alt_fit = cat(2,alt_fits_cell{isGoodRecording});
azigof  = cat(2,azi_gof_cell{isGoodRecording});
altgof  = cat(2,alt_gof_cell{isGoodRecording});
prefdir = cat(1,pref_dir_cell{isGoodRecording});


CC_mean_coherence = cat(1,CC_mean_coherence_cell{isGoodRecording});
isCoherenceResponsive = cat(1,isCoherenceResponsive_cell{isGoodRecording})';
reliability = cat(1,reliability_cell{isGoodRecording});


isEdge = round(alt_pref) == 30 | round(alt_pref) == 1;


isReliable = (max(reliability,[],2) > reliability_threshold)';

%% Recalculator isSpatiallyTuned...

narrow_thresh = 2;
fit_thresh = 0.3;

for ii = 1:length(azi_fit)
azi_width(ii) = azi_fit{ii}.c1;
alt_width(ii) = alt_fit{ii}.c1;
end
alt_rsquare = [altgof.rsquare];
azi_rsquare = [azigof.rsquare];

isTooNarrow = (azi_width<narrow_thresh) | (alt_width<narrow_thresh);
isFit = (alt_rsquare>fit_thresh) | (azi_rsquare>fit_thresh);
isEdge = round(azi_pref) == 40 | round(azi_pref) == 1;
isSpatiallyTuned = isFit & ~isTooNarrow & ~isEdge;


%% Analyse
isIncluded = isSpatiallyTuned; 

% Calculate the correlation
corr(azi_pref(isIncluded)',max(CC_mean_coherence(isIncluded,:),[],2),'Type','Pearson')
corr(alt_pref(isIncluded)',max(CC_mean_coherence(isIncluded,:),[],2),'Type','Pearson')


% Visualize? 

scatter(azi_pref(isIncluded)',max(CC_mean_coherence(isIncluded,:),[],2)')


%% Binning stuff again... round 2

alt = alt_pref(isIncluded);
cc = max(CC_mean_coherence(isIncluded,:),[],2);
isCoh = isCoherenceResponsive(isIncluded);

alt_bins = discretize(alt,0:3:30);

for ii = 1:length(unique(alt_bins))
    
    binned_resp(ii) = mean(cc(alt_bins == ii));
    std_binned_resp(ii) = std(cc(alt_bins == ii))/sqrt(sum(alt_bins == ii));
end

errorbar(1:length(binned_resp),binned_resp,std_binned_resp)


%% Direciton analysis
pref_dir = cat(1,pref_dir_cell{isGoodRecording});

%% Direction selecitivty
data = CC_mean_coherence;
[pref,pref_idx] = max(CC_mean_coherence,[],2);
wrapN = @(x, N) (1 + mod(x-1, N));
N = 8;

null_idx = wrapN(pref_idx+4,8);

for ii = 1:size(CC_mean_coherence,1)
null(ii) = CC_mean_coherence(ii,null_idx(ii));
end
DSI = (pref - null') ./ (pref + null');

isDirSelective = DSI>0.25;

%% Discarding cells with no consecutive positive binz
for t = 1:length(pref_dir)
a = CC_mean_coherence(t,:) > 0;
a0 = double(a);

ii= strfind(a0,[1 0]);% Find the end of any consecutive 1's in a0
a1 = cumsum(a);
i1 = a1(ii); % Cumulative sum at end of any consecutive 1's in a0
% Places the amount to subtract during cumulative-sum 1-element past the
% consecutive 1's in a to produce only the cumulative sum of consecutive
% 1's in a0.  If this is confusing, output a0 after this step.
try
a0(ii+1) = -[i1(1),diff(i1)];
% output vector
out = cumsum(a0); 

max_consec(t) = max(out);
catch
    max_consec(t) = 0;
end
end
pref_dir = pref_dir(max_consec>2 & isDirSelective');


bin_sz = pi/6;
% minorbinning...
dir_bins = discretize(pref_dir,[0:bin_sz:2*pi]);
temp = tabulate(dir_bins);
theta = [0:bin_sz:2*pi];
rho = cat(1,temp(:,2),temp(1,2));


% shuff
iter = 5000;
cts = zeros(iter,length([0:bin_sz:2*pi])-1);
for ii = 1:iter
    pref_dir_samp = pref_dir(randi(length(pref_dir),1,length(pref_dir)));
    dir_bins = discretize(pref_dir_samp,[0:bin_sz:2*pi]);
    temp = tabulate(dir_bins);
    cts(ii,:) = temp(:,2);
end

low_bound = prctile(cts,5);
high_bound = prctile(cts,95);

rho2 = [low_bound low_bound(1)];
rho3 = [high_bound high_bound(1)];

p = polarplot(theta,rho);
hold on
polarplot(theta,rho2,'k');%,'LineStyle','--');
polarplot(theta,rho3,'k');%,'LineStyle','--');
pax = gca;
set(p,'Marker','o','MarkerSize',5,'LineWidth',2);
set(pax,'ThetaZeroLocation','top','ThetaDir','clockwise');
%title('Directional histogram')


%{
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

%}
