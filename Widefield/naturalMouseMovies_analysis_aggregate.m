%toeFaded_aggregateAnalysis3

for ii = 1:length(naturalMouse_aggregateData)
    cc_map= rot90((cat(3,naturalMouse_aggregateData(ii).coherence_map_v2)));
    horz_ret = naturalMouse_aggregateData(ii).HorizontalRetinotopy;
    vert_ret = naturalMouse_aggregateData(ii).VerticalRetinotopy;

    temp = struct2cell(naturalMouse_aggregateData(ii).rois);
    rois = cat(3,temp{:});
    
    for ct = 1:size(cc_map,3)
        curr_map = cc_map(:,:,ct);
        
        for r = 1:size(rois,3)
            
            %bin_vert = discretize(vert_ret(rois(:,:,r)),);
            %bin_horz = discretize(horz_ret(rois(:,:,r)),3);
            
           %vert_cc(ii,ct,r) = corr(curr_map(rois(:,:,r)),bin_vert);
            %horz_cc(ii,ct,r) = corr(curr_map(rois(:,:,r)),bin_horz);
            
              vert_cc(ii,ct,r) = corr(curr_map(rois(:,:,r)),vert_ret(rois(:,:,r)));
               horz_cc(ii,ct,r) = corr(curr_map(rois(:,:,r)),horz_ret(rois(:,:,r)));
            
               
                    vert_ret_comb{ii,r} = vert_ret(rois(:,:,r));
               horz_ret_comb{ii,r} = horz_ret(rois(:,:,r));
               
            vert_cc_comb{ii,r} = zscore(curr_map(rois(:,:,r)));
            horz_cc_comb{ii,r} = zscore(curr_map(rois(:,:,r)));
            
        end
       v1_CC(ii,ct) = mean(curr_map(rois(:,:,1)));
       
       all_vert_cc(ii,ct) = corr(curr_map(:),vert_ret(:));
       all_horz_cc(ii,ct) = corr(curr_map(:),horz_ret(:));
 end
end




isJunk = v1_CC < 0.1;


for ii = 1:4
    curr_vert = vert_cc(:,:,ii);
    vert_cc_filt(:,ii) = curr_vert(~isJunk);
    
    curr_horz = horz_cc(:,:,ii);
    horz_cc_filt(:,ii) = curr_horz(~isJunk);
    
    
    
    
end
    

v1_vert = vert_cc(:,:,1);
v1_horz = horz_cc(:,:,1);

vert_v1 = v1_vert(~isJunk);
horz_v1 = v1_horz(~isJunk);


vert_cc_all = reshape(vert_cc,[],size(vert_cc,3));
horz_cc_all = reshape(horz_cc,[],size(horz_cc,3));





boxplot([vert_v1,horz_v1],[ones(1,length(vert_v1)),2*ones(1,length(horz_v1))])

hold on
for ii = 1:length(vert_v1)
    plot([1 2],[vert_v1(ii) horz_v1(ii)],'Color',[0.7 0.7 0.7])
end


% 
%  for ii = 1:length(naturalMouse_aggregateData)
%     cc_map= rot90((cat(3,naturalMouse_aggregateData(ii).coherence_map)));
%     horz_ret = naturalMouse_aggregateData(ii).HorizontalRetinotopy;
%     vert_ret = naturalMouse_aggregateData(ii).VerticalRetinotopy;
%     
%     for j = 1:3
%         imagesc(cc_map(:,:,j))
%         title(['Vert: ' num2str(vert_cc(ii,j)) ' Horz: ' num2str(horz_cc(ii,j))])
%         pause
%     end
% end

vert_cc2 = vert_cc(:);

horz_cc2 = horz_cc(:);

mean(vert_cc2)
mean(horz_cc2)

for ii = 1:length(naturalMouse_aggregateData)
    cc_map = rot90(mean(cat(3,naturalMouse_aggregateData(ii).coherence_map),3));
    cc_map(isnan(cc_map)) = 0;
    
    vert = naturalMouse_aggregateData(ii).VerticalRetinotopy;
    horz = naturalMouse_aggregateData(ii).HorizontalRetinotopy;

    tform = naturalMouse_aggregateData(ii).transformation_parameters.tform;
    Rfixed = naturalMouse_aggregateData(ii).transformation_parameters.Rfixed;
    
    tform_map(:,:,ii) = rescale(imwarp(cc_map,tform,'OutputView',Rfixed));
    
    tform_vert(:,:,ii) = imwarp(vert,tform,'OutputView',Rfixed);
    tform_horz(:,:,ii) = imwarp(horz,tform,'OutputView',Rfixed);
end    
    
mean_warped_map = mean(tform_map,3);
outline = logical(importdata('D:\Temporary Code Folder\reference_map.mat'));

mean_warped_map(outline) = max(mean_warped_map(:))*1.1;
imagesc(mean_warped_map)
