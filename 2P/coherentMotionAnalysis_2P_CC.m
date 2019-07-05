
function [] = coherentMotionAnalysis_2P_CC()

% in either case, let's take the thing and ya...

%% is something wrong with the sorting? idk why we have V1 resps that are sucky..
testing_flag=0;
% load the data
[fn_stim,pn_stim] = uigetfile('.mat');
Stimdat = importdata([pn_stim fn_stim]);
[fn_data,pn_data] = uigetfile('.mat');
data = importdata([pn_data,fn_data]);


if ~isfield(data,'spikes')
    spikeInference([pn_data fn_data],'Yes');
    data = importdata([pn_data,fn_data]);
end


repeats = Stimdat.repeats;
directions = length(Stimdat.dot_parameters.coherence_direction);
on_time = Stimdat.on_time; %replace 60 with presentations
pre_time = 1;
off_time = Stimdat.off_time - pre_time;

fs=10;
on_frames = on_time*fs;
pre_frames = pre_time*fs;
off_frames = off_time*fs;

dir_frames = on_frames+off_frames+pre_frames;
rep_frames = dir_frames * directions;


% sorting
RespVec_raw = zeros(size(data.DFF,1),on_frames,directions,repeats,'single');

for rep = 1:repeats
    for dir = 1:directions
        curr_frame = (rep-1)*rep_frames + (dir-1)*dir_frames+ pre_frames;
        off_resp(:,dir,rep) = mean(data.spikes(:,curr_frame+1 : curr_frame+off_frames),2);
        on_resp(:,:,dir,rep) = data.DFF(:,curr_frame+off_frames+1 : curr_frame+off_frames+on_frames);
        RespVec_raw(:,:,dir,rep) = on_resp(:,:,dir,rep) - off_resp(:,dir,rep); % let's keep this off for now...
    end
end
coherence = Stimdat.dot_parameters.coherence;


meanOnResp = mean(mean(mean(on_resp,4),2),3)';
meanOffDev = std(squeeze(mean(off_resp,2)),[],2)';

off_resp = squeeze(mean(off_resp,2));

%same after dis
% Magnitude of response...
temp = tabulate(coherence);
isRealCoherence = temp(:,3)>5;
unique_coherences = temp(isRealCoherence,1);

for coh = 1:length(unique_coherences)
    for rep = 1:size(RespVec_raw,4)
        for dir = 1:size(RespVec_raw,3)
            for c = 1:size(RespVec_raw,1)
                
                %    shifted_coherence = padarray(coherence,[0 lag(c)],0,'pre');
                %shift the coherence to account for lag...
                
                isCoherence = coherence == unique_coherences(coh);
                %isCoherence = coherence == unique_coherences(ii);
                % trimming
                %isCoherence = isCoherence(1:length(coherence));
                RespVec(c,dir,coh,rep) = mean(RespVec_raw(c,isCoherence,dir,rep));
            end
        end
    end
end

% from the main RespVec, we can create other respvecs....

big_coherences = find(unique_coherences > 0.4,1,'first');

RespVec_coh = squeeze(mean(RespVec,2)); % meaning across directions
RespVec_dir = squeeze(mean(RespVec(:,:,big_coherences:end,:),3)); % meaning across coherences

% sorting the directions correctly...
 [~,dir_sort] = sort(Stimdat.dot_parameters.coherence_direction);
 RespVec_dir = RespVec_dir(:,dir_sort,:);
 
% finding the preferred direction of the cell
[~,pref_dir] = max(mean(RespVec_dir,3),[],2);

%% sorting everything...

RespVec = RespVec(:,dir_sort,:,:);
RespVec_raw = RespVec_raw(:,:,dir_sort,:);

zScoreVec = meanOnResp ./ meanOffDev;


%% checking for reliability, a la Juavinett & Callaway 2015
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


isCoherenceResponsive = (p<0.05) & isVisuallyResponsive;

%% Correlation Analysis
for ii = 1:size(RespVec_raw,1)
    RespVec_directions(ii,:,:) = RespVec_raw(ii,:,pref_dir(ii),:);
    %RespVec_directions(ii,:,:) = squeeze(mean(RespVec_raw(ii,:,:,:),3));
end%%

% crosscorrelationnn
for c = 1:size(RespVec_directions,1)
    %     for rep = 1:size(RespVec_directions,3)
    %         temp = corrcoef(RespVec_directions(c,:,rep),mean(RespVec_directions(c,:,1:end~=rep),3));
    %         CC(c,rep) = temp(1,2);
    %     end
    %     [temp2,lag_temp] = xcorr(mean(RespVec_directions(c,:,:),3),coherence,10); %calculate unscaled xcorr
    %     [~,max_CC_idx] = max(temp2);  % get the best one
    %     lag(c) = lag_temp(max_CC_idx);
    %     if lag(c) < 0 % essentialy, we should not haev ANY negative lag, so this gets rid of that...
    %         max_CC_idx = find(lag_temp == 0);
    %         lag(c) = lag_temp(max_CC_idx);
    %     end
    %
    %     %%NOT TO MEAN!!!
    %     temp3 = corrcoef(mean(RespVec_directions(c,:,:),3),coherence); % find teh true CC at the 0 lag
    %     true_CC = temp3(1,2);
    %     CC_mean_raw(c) = temp2(max_CC_idx);
    %     CC_mean_coherence(c) = CC_mean_raw(c) .* true_CC/temp2(lag_temp == 0); % scale the xcorrs to the 0 corr
    %
    %     for rep = 1:size(RespVec_directions,3)
    %         temp2 = xcorr(RespVec_directions(c,:,rep),coherence,10); %calculate unscaled xcorr
    %         %%NOT TO MEAN!!!
    %         temp3 = corrcoef(RespVec_directions(c,:,rep),coherence); % find teh true CC at the 0 lag
    %         true_CC = temp3(1,2);
    %         CC_ind_raw(c,rep) = temp2(max_CC_idx);
    %         CC_ind_coherence(c,rep) = CC_ind_raw(c,rep) .* true_CC/temp2(lag_temp == 0); % scale the xcorrs to the 0 corr
    %     end
    %
    %% Coherence CC
    cell_trace = detrend(mean(RespVec_directions(c,:,:),3));
    temp = xcorr(cell_trace,coherence,10,'coeff');
    CC_mean_coherence(c) = max(temp);
    
    %% reliability
    for rep = 1:size(RespVec_directions,3)
        temp = corrcoef(RespVec_directions(c,:,rep),mean(RespVec_directions(c,:,1:end~=rep),3));
        reliability_all_reps(c,rep)  = temp(1,2);
    end
    
end

CC_mean_coherence(CC_mean_coherence>1) = 1;
CC_mean_coherence(CC_mean_coherence<-1) = -1;
reliability = mean(reliability_all_reps,2);

% for c = 1:size(RespVec_directions)
%     imagesc(squeeze(RespVec_directions(c,:,:))');
%     pause
% end


% figure; hold on
% for c = 1:size(RespVec_prefDir,1)
%     if isVisuallyResponsive(c)
%         switch isCoherenceResponsive(c)
%             case 1
%                 plot(mean(RespVec_prefDir(c,:,:),3),'r','LineWidth',2)
%             case 0
%                 plot(mean(RespVec_prefDir(c,:,:),3),'Color',[0.7 0.7 0.7])
%         end
%     end
% end

save coherenceResponseData.mat isCoherenceResponsive isVisuallyResponsive RespVec pref_dir CC_mean_coherence reliability


sum(isCoherenceResponsive)/sum(isVisuallyResponsive)
