% Step 2/2: Preprocess the data first, then use this to combine it
clear 

disp('Choose your horizontal (altitude) data...')
[fn_alt,pn_alt] = uigetfile('.mat');
disp('Choose your vertical (azimuth) data...')
[fn_azi,pn_azi] = uigetfile('.mat');

alt_data = importdata([pn_alt,fn_alt]);
azi_data = importdata([pn_azi,fn_azi]);

altitude = alt_data.RespVec;
azimuth = azi_data.RespVec;

isVisuallyResponsive = alt_data.isVisuallyResponsive | azi_data.isVisuallyResponsive;

% Initial pass, let's just mean everything and turn it into a heatmap...

altitude_m = squeeze(mean(mean(altitude,2),4));
azimuth_m = squeeze(mean(mean(azimuth,2),4));

for ii = 1:size(altitude_m,1)
    [azi_fit{ii},azigof] = fit([1:size(azimuth_m,2)]',azimuth_m(ii,:)','gauss1','Upper',[Inf,size(azimuth_m,2),Inf],'Lower',[0,1,0]);
    azi_pref(ii) = azi_fit{ii}.b1;
    azi_width(ii) = azi_fit{ii}.c1;
    azi_rsquare(ii) = azigof.rsquare;
    
    [alt_fit{ii},altgof] = fit([1:size(altitude_m,2)]',altitude_m(ii,:)','gauss1','Upper',[Inf,size(altitude_m,2),Inf],'Lower',[0,1,0]);
    alt_pref(ii) = alt_fit{ii}.b1;
    alt_width(ii) = alt_fit{ii}.c1;
    alt_rsquare(ii) = altgof.rsquare;
end

for ii = 1:size(altitude,1)
    alt_p(ii) = anova1(squeeze(mean(altitude(ii,:,:,:),2))',[],'off');
    azi_p(ii) = anova1(squeeze(mean(azimuth(ii,:,:,:),2))',[],'off');
end


%edge discarding

isEdgeAzi = false(1,size(azimuth_m,1));
isEdgeAlt = false(1,size(altitude_m,1));
    for ii = 1:size(azimuth_m,1)
        if round(azi_fit{ii}.b1,3) == 1 || round(azi_fit{ii}.b1,3) == size(azimuth_m,2)
        isEdgeAzi(ii) = true;
        end
        
        if round(alt_fit{ii}.b1,3) == 1 || round(alt_fit{ii}.b1,3) == size(altitude_m,2)
        isEdgeAlt(ii) = true;
        end
    end

    
isTooNarrow = (azi_width<1) | (alt_width<1);
isEdge = isEdgeAzi | isEdgeAlt;
isTun = (alt_p<0.05) | (azi_p<0.05);

isFit = (alt_rsquare>0.5) | (azi_rsquare>0.5);

isSpatiallyTuned = isFit & ~isTooNarrow;

disp(['Fraction spatially tuned: ' num2str(mean(isSpatiallyTuned))])


%% extracting centroids
%load your data file...
numCells = size(azi_data.cellMasks,2);
for c = 1:numCells
    currCellMask = azi_data.cellMasks{c};
   temp = regionprops(poly2mask(currCellMask(:,1),currCellMask(:,2),760,760),'centroid');
   roi_centroids(c,:) = temp.Centroid;
end
   
isANOVA = alt_p<0.01 | azi_p<0.01;

isSpat = isSpatiallyTuned;% & isANOVA;
x_locations = roi_centroids(:,1);
y_locations = roi_centroids(:,2);

alt_corr = corr(alt_pref(isSpat)',x_locations(isSpat));
azi_corr = corr(azi_pref(isSpat)',y_locations(isSpat));


save RFmapping_results.mat azi_pref alt_pref altitude azimuth isSpatiallyTuned alt_p azi_p azi_fit alt_fit roi_centroids alt_corr azi_corr


% %%testing flag
% if testing_flag
%     %for altitude
%     monitor_width = 81; % in degrees
%     monitor_height = 67; % in degrees    % assuming 30 deg RF
% 
%   bar_deg = monitor_height/size(altitude_m,2); % each bar center is at... 
%   top_edge =  size(altitude_m,2)/2;
%   top_bar = round(top_edge/bar_deg);
%   bot_edge = monitor_height-size(altitude_m,2)/2;
%   bot_bar = round(bot_edge/bar_deg);
%   scaling_vector_alt = [linspace(2,1,8), ones(1,14),linspace(1,2,8)];
%   altitude_m = altitude_m(:,top_bar:bot_bar);
% 
%   bar_deg = monitor_width/size(azimuth_m,2); % each bar center is at... 
%   left_edge =  size(azimuth_m,2)/2;
%   left_bar = round(left_edge/bar_deg);
%   right_edge = monitor_width-size(azimuth_m,2)/2;
%   right_bar = round(right_edge/bar_deg);
%   scaling_vector_azi = [linspace(2,1,8), ones(1,24),linspace(1,2,8)];
% azimuth_m = azimuth_m(:,left_bar:right_bar);  
%   
% end


% 
% for c = 1 :427
%     plot(altitude_m(c,:));
%     hold on
%     plot(alt_fit{c})
%     hold off
%     title(num2str(alt_rsquare(c)))
%     pause
% end



% 
% 
% % %%%%%%%%%%%%%%%%%%
% % For visualization only
% for ii = 1:length(isVisuallyResponsive)
%     horz_temp = altitude_m(ii,:);
%     vert_temp = azimuth_m(ii,:);
%     [~,pref_vert(ii)] = max(vert_temp);
%     [~,pref_horz(ii)] = max(horz_temp);
%     RF_matrix(:,:,ii) = bsxfun(@plus,horz_temp,vert_temp');
%     
%     % RF_undistorted(:,:,ii) = imageUnDistorter(RF_matrix(:,:,ii));
% end
% 
% % Visualize
% for ii = 1:length(isVisuallyResponsive)
%     if isSpatiallyTuned(ii)
%         imagesc(imgaussfilt(RF_matrix(:,:,ii),2)) % ligth smoothing for beauty reasons
%         pbaspect([800 600 1])
%         % subplot(1,2,2)
%         % imagesc(imgaussfilt(RF_undistorted(:,:,ii),2))
%         pause
%     end
% end
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % playin' around with the data
% 
% for ii = 1:size(altitude_m,1)
%     azi_fit = fit([1:size(azimuth_m,2)]',azimuth_m(ii,:)','gauss1','Upper',[10,size(azimuth_m,2),10],'Lower',[-10,1,-10]);
%     azi_CC(:,ii) = azi_fit(1:size(azimuth_m,2)) * mean(CC_ind_coherence(ii,:),2);
%     
%   alt_fit = fit([1:size(altitude_m,2)]',altitude_m(ii,:)','gauss1','Upper',[10,size(altitude_m,2),10],'Lower',[-10,1,-10]);
%     alt_CC(:,ii) = alt_fit(1:size(altitude_m,2)) * mean(CC_ind_coherence(ii,:),2);
%     
% end     
%      
     
% so the idea is to like bin values and then mean the CCs of neurons within
% each bin, then see wat we get ya???


     
     
% 
% corrcoef(azi_pref(isCoherenceResponsive),nanmean(CC_ind_coherence(isCoherenceResponsive,:),2))
% corrcoef(alt_pref(isCoherenceResponsive),nanmean(CC_ind_coherence(isCoherenceResponsive,:),2))
% 
