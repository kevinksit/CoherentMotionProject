function ReceptiveField_processingCombined
%Step 1/2: Process each direction using this first, then we'll combine it
%later

% Get the information from the Stimdata
[fn_Stimdat,pn_Stimdat] = uigetfile('.mat');
[fn_data,pn_data] = uigetfile('.mat');

Stimdata = importdata([pn_Stimdat,fn_Stimdat]);
data     = importdata([pn_data, fn_data]);

if ~isfield(data,'spikes')
    data = spikeInference([pn_data '\' fn_data],'Yes');
end

% Get the data that you need

zThresh = 2;

on_time = Stimdata.on_time;
off_time = Stimdata.off_time;
repeats = Stimdata.repeats;
bar_centers = Stimdata.centers;
num_locations = size(bar_centers,2);
fs = 10;

% Convert to frames
on_frames = on_time *fs;
off_frames = off_time*fs;

% Calculate additional times
rep_frames = off_frames + (on_frames * num_locations);

% Sort the data appropriately (reps and shuffled) (cells, on_frames, locs,reps)
RespVec = zeros(size(data.DFF,1),on_frames,num_locations,repeats);

for rep = 1:repeats
    curr_frame = (rep-1)*rep_frames;
    off_resp(:,:,rep) = data.spikes(:,curr_frame+1 : curr_frame+off_frames);
    for loc = 1:num_locations
        curr_frame = (rep-1)*rep_frames + (loc-1)*on_frames + off_frames;
        on_resp(:,:,loc,rep) = data.spikes(:,curr_frame+1:curr_frame+on_frames);
        RespVec(:,:,loc,rep) = on_resp(:,:,loc,rep) - mean(off_resp(:,:,rep),2);
    end
end


% Unscramble the data using Stimdat

for rep = 1:repeats
    [~,sorting_vector] = sort(bar_centers(rep,:));
    on_resp(:,:,:,rep) = on_resp(:,:,sorting_vector,rep);
end

% Check for visual responsiveness
off_dev = std(mean(off_resp,3),[],2); % mean across frames and rep
mean_on_resp = squeeze(mean(mean(on_resp,4),2)); % mean across frames and rep
max_on_resp = max(mean_on_resp,[],2); % preferred location resp
zScoreVec = max_on_resp ./ off_dev;
isVisuallyResponsive = zScoreVec > zThresh;


% sorting... 

for rep = 1:repeats
    isAzi = Stimdata.aziOrAlt(rep,:) == 1;
    temp = RespVec(:,:,isAzi,rep);
    [~,azi_sort] = sort(bar_centers(rep,isAzi));
    azimuth(:,:,:,rep) = temp(:,:,azi_sort);
    
    isAlt = Stimdata.aziOrAlt(rep,:) == 2;
    temp = RespVec(:,:,isAlt,rep);
    [~,alt_sort] = sort(bar_centers(rep,isAlt));
    altitude(:,:,:,rep) = temp(:,:,alt_sort);
end

% Initial pass, let's just mean everything and turn it into a heatmap...

altitude_m = squeeze(mean(mean(altitude,2),4));
azimuth_m = squeeze(mean(mean(azimuth,2),4));



for ii = 1:size(altitude_m,1)
    [azi_fit{ii},azigof] = fit([1:size(azimuth_m,2)]',azimuth_m(ii,:)','gauss1','Upper',[20,size(azimuth_m,2),20],'Lower',[0,1,0]);
    azi_pref(ii) = azi_fit{ii}.b1;
    azi_width(ii) = azi_fit{ii}.c1;
    azi_rsquare(ii) = azigof.rsquare;
    
    [alt_fit{ii},altgof] = fit([1:size(altitude_m,2)]',altitude_m(ii,:)','gauss1','Upper',[20,size(altitude_m,2),20],'Lower',[0,1,0]);
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
numCells = size(data.cellMasks,2);
for c = 1:numCells
    currCellMask = data.cellMasks{c};
   temp = regionprops(poly2mask(currCellMask(:,1),currCellMask(:,2),760,760),'centroid');
   roi_centroids(c,:) = temp.Centroid;
end
   
isANOVA = alt_p<0.01 | azi_p<0.01;

isSpat = isSpatiallyTuned & isANOVA;
x_locations = roi_centroids(:,1);
y_locations = roi_centroids(:,2);

alt_corr = corr(alt_pref(isSpat)',x_locations(isSpat));
azi_corr = corr(azi_pref(isSpat)',y_locations(isSpat));

cd ..
save RFmapping_results.mat azi_pref alt_pref altitude azimuth isSpatiallyTuned alt_p azi_p azi_fit alt_fit roi_centroids alt_corr azi_corr

