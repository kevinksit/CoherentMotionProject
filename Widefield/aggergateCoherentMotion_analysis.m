% for processing widefield coherent motion recordings
load(uigetfile);

for rec = 1:length(coherentMotion_aggregateData)
   v1_roi = coherentMotion_aggregateData(rec).rois.V1;
   pm_roi = coherentMotion_aggregateData(rec).rois.PM;
   temp = regionprops(pm_roi,'centroid');
   pm_centroid = temp.Centroid;
 
   for y = 1:400
       for x = 1:400
           distance_map(y,x) = sqrt((y-pm_centroid(2))^2 + (x - pm_centroid(1))^2);
       end
   end
   curr_coherence = coherentMotion_aggregateData(rec).coherence_map;
   curr_vertret = coherentMotion_aggregateData(rec).vertical_retinotopy;
   curr_horzret = coherentMotion_aggregateData(rec).horizontal_retinotopy;
   curr_magmap = rot90(mean(coherentMotion_aggregateData(rec).coherence_slope, 3));
   temp = corrcoef(distance_map(v1_roi),curr_vertret(v1_roi));
   distance_vertCC(rec) = temp(1,2);
   
   
   temp = corrcoef(distance_map(v1_roi),curr_coherence(v1_roi));
   distanceCC(rec) = temp(1,2);
   
   temp = corrcoef(curr_vertret(v1_roi),curr_coherence(v1_roi));
   vertCC(rec) = temp(1,2);
   temp = corrcoef(curr_horzret(v1_roi),curr_coherence(v1_roi));
   horzCC(rec) = temp(1,2);
   temp = corrcoef(curr_magmap(v1_roi),curr_coherence(v1_roi));
   magCC(rec) = temp(1,2);
   
   distances{rec} = distance_map(v1_roi);
   retinotopies{rec} = curr_vertret(v1_roi);
   coherences{rec} = curr_coherence(v1_roi);
end

%% Comparing correlation to the coherence trace
for rec = 1:length(coherentMotion_aggregateData)
    current_map = coherentMotion_aggregateData(rec).coherence_map;
    % extract the info
    v1_resp_pix{rec} = current_map(coherentMotion_aggregateData(rec).rois.V1);
    lm_resp_pix{rec} = current_map(coherentMotion_aggregateData(rec).rois.LM);
    al_resp_pix{rec} = current_map(coherentMotion_aggregateData(rec).rois.AL);
    pm_resp_pix{rec} = current_map(coherentMotion_aggregateData(rec).rois.PM);
end

    
%some kind of "homogeneity" standardized response... mean / std, works...
%but let's get something betta?
v1_resp = cellfun(@(x) mean(x),v1_resp_pix);
lm_resp = cellfun(@(x) mean(x),lm_resp_pix);
al_resp = cellfun(@(x) mean(x),al_resp_pix);
pm_resp = cellfun(@(x) mean(x),pm_resp_pix);


boxplot([v1_resp; lm_resp; al_resp; pm_resp]')
hold on
for ii = 1:length(v1_resp)
     plot([v1_resp(ii) lm_resp(ii) al_resp(ii) pm_resp(ii)],'--','Color',[0.7 0.7 0.7],'LineWidth',1)
end
xticklabels({'V1' 'LM' 'AL' 'PM'})
ylabel('Correlation to Coherence Trace')
scatter(v1_resp,al_resp)

%% Comparing magnitude of response as a function of coherence
for rec = 1:length(coherentMotion_aggregateData)
    for coh = 1:5
    current_map = coherentMotion_aggregateData(rec).coherence_curve(:,:,coh);
    % extract the info
    v1_curve_pix{rec,coh} = current_map(coherentMotion_aggregateData(rec).rois.V1);
    lm_curve_pix{rec,coh} = current_map(coherentMotion_aggregateData(rec).rois.LM);
    al_curve_pix{rec,coh} = current_map(coherentMotion_aggregateData(rec).rois.AL);
    pm_curve_pix{rec,coh} = current_map(coherentMotion_aggregateData(rec).rois.PM);
    end
end

v1_curve = cellfun(@mean,v1_curve_pix);
lm_curve = cellfun(@mean,lm_curve_pix);
al_curve = cellfun(@mean,al_curve_pix);
pm_curve = cellfun(@mean,pm_curve_pix);

all_curves = cat(2,v1_curve,lm_curve,al_curve,pm_curve);


for ii = 1:size(all_curves,1)
    all_curves_norm(ii,:) = rescale(all_curves(ii,:));
end
all_curves_reshaped = reshape(all_curves_norm,10,5,4);

all_curves_meaned = squeeze(mean(all_curves_reshaped,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% %% transform onto the CCF

for rec = 1:length(coherentMotion_aggregateData)
    current_map = coherentMotion_aggregateData(rec).coherence_map;
    curr_map2 = coherentMotion_aggregateData(rec).cortical_map;
    curr_map3 = coherentMotion_aggregateData(rec).vertical_retinotopy;
    curr_map4 = coherentMotion_aggregateData(rec).horizontal_retinotopy;
    tform = coherentMotion_aggregateData(rec).transformation_parameters.tform;
    Rfixed = coherentMotion_aggregateData(rec).transformation_parameters.Rfixed;
    warped_maps(:,:,rec) =  (imwarp(current_map,tform,'OutputView',Rfixed));
    warped_signmaps(:,:,rec) = (imwarp(curr_map2,tform,'OutputView',Rfixed));
    warped_horzret(:,:,rec) = imwarp(curr_map3,tform,'OutputView',Rfixed);
    warped_vertret(:,:,rec) = imwarp(curr_map4,tform,'OutputView',Rfixed);
end

% get the CCF
outline = logical(importdata('D:\Temporary Code Folder\reference_map.mat'));

mean_warped_maps = mean(warped_maps,3);
mean_warped_maps = imgaussfilt(mean_warped_maps,10);
mean_warped_maps(outline) = max(mean_warped_maps(:))*1.1;

imagesc(mean_warped_maps)



scatter(pm_resp,v1_resp,[],[1 0.4 0.4],'filled')
axis square
axis([0 0.7 0 0.7])
hline = refline([1,0]);
hline.Color = [0.5 0.5 0.5];
hline.LineStyle = '--';




