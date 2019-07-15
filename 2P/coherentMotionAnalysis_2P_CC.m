
function [] = coherentMotionAnalysis_2P_CC(data,Stimdat)

%% For analysis of coherent motion data from 2P
% Changelog Updated 09Jul2019 KS, cleaned up pretty significantly
%           Updated 11Jul2019 KS, added code to use a single lag correction
%           Updated 15Jul2019 KS, pref_dir calculated via CC instead of magnitude, and plotting supported, added input args

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define necessary parameters
visual_threshold    = 2; % zscore to be considered visually responsive
coherence_threshold = 0.05; % p value to be considered coherence responsive
maximum_offset      = 20; % maximum offset allowed for 
plot_flag           = 1; % visualize the data or no?

%% Load data
if nargin == 0
[fn_stim,pn_stim] = uigetfile('.mat');
Stimdat = importdata([pn_stim fn_stim]);
[fn_data,pn_data] = uigetfile('.mat');
data = importdata([pn_data,fn_data]);
end

% Preprocessing, turning into spikes
if ~isfield(data,'spikes')
    spikeInference([pn_data fn_data],'Yes');
    data = importdata([pn_data,fn_data]);
end

% Extract the necessary components
repeats = Stimdat.repeats;
directions = length(Stimdat.dot_parameters.coherence_direction);
on_time = Stimdat.on_time; 
pre_time = 1;
off_time = Stimdat.off_time - pre_time;

fs=10;
on_frames = on_time*fs;
pre_frames = pre_time*fs;
off_frames = off_time*fs;

dir_frames = on_frames+off_frames+pre_frames;
rep_frames = dir_frames * directions;

numCells = size(data.DFF,1);



%% Sorting data into response matrices

% Raw responses (neural traces)
RespVec_raw = zeros(numCells,on_frames,directions,repeats,'single');
off_resp = zeros(numCells,directions,repeats,'single');
on_resp = zeros(numCells,on_frames,directions,repeats,'single');

for rep = 1:repeats
    for dir = 1:directions
        curr_frame = (rep-1)*rep_frames + (dir-1)*dir_frames+ pre_frames;
        off_resp(:,dir,rep) = mean(data.spikes(:,curr_frame+1 : curr_frame+off_frames),2);
        on_resp(:,:,dir,rep) = data.DFF(:,curr_frame+off_frames+1 : curr_frame+off_frames+on_frames);
        RespVec_raw(:,:,dir,rep) = on_resp(:,:,dir,rep) - off_resp(:,dir,rep);
    end
end

% Blocked response vector, with single numbers per response bin

coherence = Stimdat.dot_parameters.coherence; 
temp = tabulate(coherence); %only taking real coherences, to account for the drift time
isRealCoherence = temp(:,3)> Stimdat.transition_time * fs; 
unique_coherences = temp(isRealCoherence,1);

RespVec = zeros(numCells,directions,length(unique_coherences),repeats,'single');

for coh = 1:length(unique_coherences)
    for rep = 1:size(RespVec_raw,4)
        for dir = 1:size(RespVec_raw,3)
            for c = 1:size(RespVec_raw,1)
                isCoherence = coherence == unique_coherences(coh);
                RespVec(c,dir,coh,rep) = mean(RespVec_raw(c,isCoherence,dir,rep)); % sort into the right bin
            end
        end
    end
end

%% Finding preferred direction of cells...
big_coherences = find(unique_coherences > 0.4, 1,'first');
RespVec_dir = squeeze(mean(RespVec(:,:,big_coherences:end,:),3)); % Only taking the top coherences, because others won't have signal...

% unscrambling directions
[~,dir_sort] = sort(Stimdat.dot_parameters.coherence_direction);

% Sorting the other two response vectors
RespVec = RespVec(:,dir_sort,:,:);
RespVec_raw = RespVec_raw(:,:,dir_sort,:);

%% Calculating visual responsiveness
meanOnResp = mean(mean(mean(on_resp,4),2),3);
meanOffDev = std(squeeze(mean(off_resp,2)),[],2);
isVisuallyResponsive = (meanOnResp ./ meanOffDev) > visual_threshold;

%% Correlation analysis

% Calculate cross correlation for each cell
CC_mean_coherence = zeros(numCells,1);
p_values          = zeros(numCells,1);
reliability       = zeros(numCells,1); 

RespVec_summed = zeros(1,size(RespVec_raw,2));

% Determine the offset for the dataset
for c = 1:size(RespVec_dir,1)
    [~,temp_pref] = max(mean(RespVec_dir(c,:,:),3));
    RespVec_summed = mean(RespVec_raw(c,:,temp_pref,:),4) + RespVec_summed;
end

[~,idx] = max(xcorr(RespVec_summed,coherence,maximum_offset));

% This is the amount of lag... should always be positive
lag = (idx) - maximum_offset;

if lag < 0
    fprintf('Warning... something''s wrong, your lag is negative! Setting lag to 0... \n');
    lag = 0;
end

% Correct for lag in the coherence trace
coherence_lagCorrected = padarray(coherence,[0 lag],'pre','replicate');
coherence_lagCorrected = coherence_lagCorrected(1:end-lag);

for c = 1:size(RespVec_raw,1)
    for dir = 1:size(RespVec_raw,3)
    % calculate the cross correlation
    cell_trace = detrend(mean(RespVec_raw(c,:,dir,:),4));
    [cc , p] = corrcoef(coherence_lagCorrected,cell_trace); % calculate correlation
    
    CC_mean_coherence(c,dir) = cc(1,2);
    p_values(c,dir)          = p(1,2);
    
    % calculate reliability
    reliability_all_reps = zeros(repeats,1);
    for rep = 1:size(RespVec_raw,4)
        temp = corrcoef(RespVec_raw(c,:,dir,rep),mean(RespVec_raw(c,:,dir,1:end~=rep),4));
        reliability_all_reps(rep)  = temp(1,2);
    end
    
    reliability(c,dir) = mean(reliability_all_reps);
    end
end

% Calculating final stuff
[~,pref_dir] = max(CC_mean_coherence,[],2);
isCoherenceResponsive = (p_values < coherence_threshold) & isVisuallyResponsive;


%save the data
save coherenceResponseData.mat isCoherenceResponsive isVisuallyResponsive RespVec pref_dir CC_mean_coherence reliability 


fprintf('Percent coherence responsive: %.1f%% \n',round(mean(isCoherenceResponsive(isVisuallyResponsive))*100,2)); % Here, the %.1f says, floating point, with 1 decimal point precision

%% Display the data...
if plot_flag 
    figure
    
    subplot(1,2,1)
    histogram(max(CC_mean_coherence,[],2));
    hold on
    line(repmat(mean(max(CC_mean_coherence,[],2)),1,2),ylim,'LineWidth',2,'Color','r'); 
    title('Coherence CC')
    ylabel('Number Cells')
    xlabel('Coherence CC')
    
    subplot(1,2,2)
    temp = tabulate(pref_dir);
    theta = cat(1,0,temp(:,1))*(pi/4);
    rho = cat(1,temp(:,2),temp(1,2));
    p = polarplot(theta,rho);
    pax = gca;
    set(p,'Marker','o','MarkerSize',5,'LineStyle','--','LineWidth',2);
    set(pax,'ThetaZeroLocation','top','ThetaDir','clockwise');
    title('Directional histogram')
end
    
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD JUNK CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%

%This was for checking for coherence via slope, it didn't work too hot... 09Jul2019

% checking for reliability, a la Juavinett & Callaway 2015
for c = 1:size(RespVec_coh,1)
    RespVec_prefDir(c,:,:) = squeeze(RespVec(c,pref_dir(c),:,:));
    %  dprime(c) = (mean(RespVec_coh(c,:,:),3) - mean(off_resp(c,:),2))/(std(RespVec_coh(c,:,:),[],3) + std(off_resp(c,:),[],2));
    
    unique_coherences_expanded = repmat(unique_coherences,[repeats,1]);
    RespVec_expanded = reshape(RespVec_prefDir(c,:,:),1,[]);
    temp_lm = fitlm(unique_coherences_expanded,RespVec_expanded);
    p(c) = temp_lm.Coefficients.pValue(2);%anova1(squeeze(RespVec_prefDir(c,:,:))',[],'off');
    slope(c) = temp_lm.Coefficients.Estimate(2);
end

isVisuallyResponsive = zScoreVec > 2;% & dprime > 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Old version of finding cross correlation... 09Jul2019


       for rep = 1:size(RespVec_directions,3)
            temp = corrcoef(RespVec_directions(c,:,rep),mean(RespVec_directions(c,:,1:end~=rep),3));
            CC(c,rep) = temp(1,2);
        end
        [temp2,lag_temp] = xcorr(mean(RespVec_directions(c,:,:),3),coherence,10); %calculate unscaled xcorr
        [~,max_CC_idx] = max(temp2);  % get the best one
        lag(c) = lag_temp(max_CC_idx);
        if lag(c) < 0 % essentialy, we should not haev ANY negative lag, so this gets rid of that...
            max_CC_idx = find(lag_temp == 0);
            lag(c) = lag_temp(max_CC_idx);
        end
    
        %%NOT TO MEAN!!!
        temp3 = corrcoef(mean(RespVec_directions(c,:,:),3),coherence); % find teh true CC at the 0 lag
        true_CC = temp3(1,2);
        CC_mean_raw(c) = temp2(max_CC_idx);
        CC_mean_coherence(c) = CC_mean_raw(c) .* true_CC/temp2(lag_temp == 0); % scale the xcorrs to the 0 corr
    
        for rep = 1:size(RespVec_directions,3)
            temp2 = xcorr(RespVec_directions(c,:,rep),coherence,10); %calculate unscaled xcorr
            %%NOT TO MEAN!!!
            temp3 = corrcoef(RespVec_directions(c,:,rep),coherence); % find teh true CC at the 0 lag
            true_CC = temp3(1,2);
            CC_ind_raw(c,rep) = temp2(max_CC_idx);
            CC_ind_coherence(c,rep) = CC_ind_raw(c,rep) .* true_CC/temp2(lag_temp == 0); % scale the xcorrs to the 0 corr
        end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%


%Was for plotting: 

for c = 1:size(RespVec_directions)
    imagesc(squeeze(RespVec_directions(c,:,:))');
    pause
end


figure; hold on
for c = 1:size(RespVec_prefDir,1)
    if isVisuallyResponsive(c)
        switch isCoherenceResponsive(c)
            case 1
                plot(mean(RespVec_prefDir(c,:,:),3),'r','LineWidth',2)
            case 0
                plot(mean(RespVec_prefDir(c,:,:),3),'Color',[0.7 0.7 0.7])
        end
    end
end

%}


