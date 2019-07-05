% coherentMotionAnalysis
clear;

%% is something wrong with the sorting? idk why we have V1 resps that are sucky..
testing_flag=0;
% load the data
[fn_stim,pn_stim] = uigetfile('.mat');
Stimdat = importdata([pn_stim fn_stim]);
data = matfile('DFF.mat');

repeats = Stimdat.repeats;
on_time = Stimdat.on_time*Stimdat.presentations %* Stimdat.fade_rate; %replace 60 with presentations
pre_time = 2;
off_time = 3;%Stimdat.off_time - pre_time;

fs=10;
on_frames = on_time*fs;
pre_frames = pre_time*fs;
off_frames = off_time*fs;

rep_frames = on_frames+off_frames+pre_frames;


% %% removing entire_frame_activity
% whole_frame_trace = squeeze(mean(mean(data.DFF,1),2));
% 
% for frame = 1:size(data,'DFF',3)
%     DFF_subtracted(:,:,frame) = data.DFF(:,:,frame) - whole_frame_trace(frame);
% end

% extracting coherence values
RespVec = zeros(size(data,'DFF',1),size(data,'DFF',2),on_frames,repeats,'single');

for rep = 1:repeats
    curr_frame = (rep-1)*rep_frames + pre_frames;
    off_resp = mean(data.DFF(:,:,curr_frame+1 : curr_frame+off_frames),3);
    on_resp = data.DFF(:,:,curr_frame+off_frames+1 : curr_frame+off_frames+on_frames);
    RespVec(:,:,:,rep) = on_resp - off_resp; % let's keep this off for now...
end

coherence = Stimdat.dot_parameters.coherence;
direction = Stimdat.dot_parameters.coherence_direction;

unique_directions = unique(direction);
%coh = coherence(1:length(coherence)/length(unique_directions));

for x = 1:size(RespVec,1)
    for y = 1:size(RespVec,2)
        disp([num2str(y) ',' num2str(x)])
%         for rep = 1:size(RespVec,4)
%              reliability(y,x,rep) = corr(squeeze(RespVec(y,x,:,rep)),squeeze(mean(RespVec(y,x,:,1:end ~=rep),4)));
%              coherenceIndCC(y,x,rep) = corr(squeeze(RespVec(y,x,:,rep)),coherence'); 
%         end
        coherenceMeanCC(y,x) = corr(squeeze(mean(RespVec(y,x,:,:),4)),coherence');
    end
end


save maps.mat  coherenceMeanCC
% 
% RespVec_meaned = mean(RespVec,4);
% 
% 
% % sorting into directions....
% 
% %% with the first couple recordings, it won't be perfect because we have gradually changing directions... we'll fix this later...
% direction_tabulated = tabulate(direction);
% [~,max_idx] = maxk(direction_tabulated(:,3),8); % 4 directions... find a better way for this later please
% unique_directions = direction_tabulated(max_idx,1);
% direction_frames = length(direction)/length(unique_directions); % this is a temporary fix fo the issue of 950 instead of 960 coherence...
% %direction_frames = 110;
% for ii = 1:length(unique_directions)
%     curr_frame = (ii-1)*120;
%     RespVec_directions(:,:,:,ii,:) = RespVec(:,:,curr_frame+1:curr_frame+direction_frames,:);
% end
% 
% 
% RespVec2 = mean(RespVec_directions,5);
% 
% 
% 
% %% sorting the directions now...
% for ii = 1:length(unique_directions)
%     direction_plot = double(direction == unique_directions(ii));
%     direction_plot_padded = padarray(direction_plot,[0,5]);
% 
%     [~,peak_loc(ii)] = findpeaks(direction_plot_padded,'MinPeakWidth',20); 
% end
% [~,sorting_vector] = sort(peak_loc);
% 
% RespVec_sorted = RespVec2(:,:,:,sorting_vector);
% % done
% 
% % finding the offset value for each pixel... 07Oct2018
% 
% for y = 1:size(RespVec,1)
%     for x = 1:size(RespVec,2)
%         for dir = 1:length(unique_directions)
%         cross_corr = xcorr(squeeze(RespVec_sorted(y,x,:,dir))',coherence,10);
%         [~,peak_corr(y,x,dir)] = max(cross_corr);
%         end
%     end
% end
% delay_map = peak_corr;
%  offset =0;% mode(peak_corr(:)) - 10;
% 
% 
% for y = 1:size(RespVec,1)
%     for x = 1:size(RespVec,2)
%       for  dir = 1:length(unique_directions)
%           CC = corrcoef(coh,squeeze(RespVec_sorted(y,x,:,dir)));
%        % CC = corrcoef(coherence(1:end-offset),squeeze(RespVec_sorted(y,x,offset+1:end,dir)));
%         coherenceMap(y,x,dir) = CC(1,2);
%       end
%     end
% end
% 
% 
% % figure;
% % for ii = 1:8
% %     subplot(1,3,1)
% %     imagesc(coherenceMap_noSubtraction(:,:,ii))
% %     title('No offresp sub');
% %     axis square
% %     
% %     subplot(1,3,2)
% %     imagesc(coherenceMap_subtracted(:,:,ii))
% %     title('wholeframesub');
% %     axis square
% %     
% %     subplot(1,3,3)
% %     imagesc(coherenceMap_original(:,:,ii))
% %     title('normal');
% %     axis square
% %     
% %     pause
% % end
% 
% % 
% % x_pix = 148;
% % y_pix = 68;
% % 
% % sort_dir = sort(unique_directions);
% % for ii = 1:length(unique_directions)
% %     plot(squeeze(RespVec_sorted(y_pix,x_pix,:,ii)));
% %     hold on
% %     plot(coherence)
% %     hold off
% %     title([num2str(coherenceMap(y_pix,x_pix,ii)) ' | ' num2str(sort_dir(ii)) 'deg'])
% %     pause
% % end
% 
% coherenceMap = rot90(coherenceMap);
% delay_map = rot90(delay_map);
% save coherenceMap.mat coherenceMap delay_map;
% 
% % 
% % 
% % unique_coherences = unique(coherence);
% % unique_directions = unique(direction);
% % 
% % 
% % RespVec_sorted = zeros(400,400,20,5,4);
% % for coh = 1:length(unique_coherences)
% %     for dir = 1:length(unique_directions)
% %         isDirection = direction == unique_directions(dir);
% %         isCoherence =coherence == unique_coherences(coh);
% %         
% %         RespVec_sorted(:,:,:,coh,dir) = RespVec(:,:,isDirection & isCoherence);
% %     end
% % end
% % direction_resp = squeeze(mean(RespVec_sorted(:,:,:,5,:),3)); % response at max coherence
% % [~,preferred_direction] = max(direction_resp,[],3);
% % save OT_Map.mat preferred_direction;
% 
% 
% 
% % 
% % % RespVec is like the heatmaps, contains all the stuff, now let's calculate
% % % CC across reps per pixel
% % 
% % CC_allreps = zeros(size(RespVec,1),size(RespVec,2),repeats);
% % % this part kinda takes forever...
% % for y = 1 : size(RespVec,1)
% %     for x = 1:size(RespVec,2)
% %         for rep = 1:repeats
% %         CC = corrcoef(RespVec(y,x,:,rep),mean(RespVec(y,x,:,1:end ~=rep),4));
% %         CC_allreps(y,x,rep) = CC(1,2);
% %         end
% %     end
% % end
% % 
% % reliabilityMap = rot90(mean(CC_allreps,3));
% % save reliabilityMap.mat reliabilityMap
% 
% %% other analyses....
% % 
% % map2 = rot90(coherenceMap,-1);
% % imagesc(map2)
% % temp = impoly;
% % roi = temp.createMask;
% % 
% % for ii = 1:size(RespVec_meaned,3)
% %     curr_frame = RespVec_meaned(:,:,ii);
% %     resp(ii) = mean(curr_frame(roi));
% % end
% % 
% 
% 
% % magnitued stuff, doesn't work right now...
% %{
% unique_coherences = unique(coherence);
% for ii = 1:length(unique_coherences)
%     for x = 1:size(RespVec_meaned,2)
%         for y = 1:size(RespVec_meaned,1)
%             coherence_responses(y,x,ii) = mean(RespVec_meaned(y,x,coherence == unique_coherences(ii)));
%         end
%     end
% end
% whole_window_bias = mean(mean(coherence_responses,1),2);
% coherence_responses = coherence_responses - whole_window_bias;
% 
% save coherence_curve.mat coherence_responses;
% %}
