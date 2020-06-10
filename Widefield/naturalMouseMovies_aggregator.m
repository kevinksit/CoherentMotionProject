%load your aggregate data

addpath(genpath('C:\Users\sit\Dropbox\CodeInBeta_Kevin\Other'))

load(uigetfile)
recording_number = length(naturalMouse_aggregateData)+1;


parent_directory = uigetdir;
first = cd(parent_directory);

temp = cd('..');
load('VFS_meancat.mat')

try
    load('transformation_parameters.mat')
    registered_flag = 1;
catch
    registered_flag = 0;
end
    
cd(temp);

naturalMouse_aggregateData(recording_number).recording_name = parent_directory(4:end);

% choose the folder you're interested in
% add the entire folder to your path
matfiles = dir('**/*.mat');


for files = 1:length(matfiles)
    if strcmp(matfiles(files).name,'multi_movie_coherence_xcorr.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    elseif strcmp(matfiles(files).name,'reliability_map.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    end
end


naturalMouse_aggregateData(recording_number).cortical_map = VFS_raw;
naturalMouse_aggregateData(recording_number).VerticalRetinotopy = maps.VerticalRetinotopy;
naturalMouse_aggregateData(recording_number).HorizontalRetinotopy = maps.HorizontalRetinotopy;
- 
% run the registerer
if ~registered_flag
    [tform, Rfixed] = referenceMap_registration(VFS_raw);
end

    naturalMouse_aggregateData(recording_number).transformation_parameters.tform = tform;
    naturalMouse_aggregateData(recording_number).transformation_parameters.Rfixed = Rfixed;


for ii = 1:length(naturalMouse_aggregateData)-1
    same_rec(ii) = contains(naturalMouse_aggregateData(ii).recording_name,...
        naturalMouse_aggregateData(recording_number).recording_name(1:end-10));
end

if ~any(same_rec)
    % get the ROIs
    disp('V1, LM, AL, PM')
    [v1_roi, lm_roi, al_roi, pm_roi] = HVAroiChooserGUI(VFS_raw);
    disp('RL, AM, LI, S1')
    [rl_roi, am_roi, li_roi, s1_roi] = HVAroiChooserGUI(VFS_raw);
    
    
    naturalMouse_aggregateData(recording_number).rois.V1 = v1_roi;
    naturalMouse_aggregateData(recording_number).rois.LM = lm_roi;
    naturalMouse_aggregateData(recording_number).rois.AL = al_roi;
    naturalMouse_aggregateData(recording_number).rois.PM = pm_roi;
    naturalMouse_aggregateData(recording_number).rois.RL = rl_roi;
    naturalMouse_aggregateData(recording_number).rois.AM = am_roi;
    naturalMouse_aggregateData(recording_number).rois.LI = li_roi;
    naturalMouse_aggregateData(recording_number).rois.S1 = s1_roi;
else
    rec = find(same_rec, 1, 'first');
    naturalMouse_aggregateData(recording_number).rois.V1 = naturalMouse_aggregateData(rec).rois.V1;
    naturalMouse_aggregateData(recording_number).rois.LM = naturalMouse_aggregateData(rec).rois.LM;
    naturalMouse_aggregateData(recording_number).rois.AL = naturalMouse_aggregateData(rec).rois.AL;
    naturalMouse_aggregateData(recording_number).rois.PM = naturalMouse_aggregateData(rec).rois.PM;
    naturalMouse_aggregateData(recording_number).rois.RL = naturalMouse_aggregateData(rec).rois.RL;
    naturalMouse_aggregateData(recording_number).rois.AM = naturalMouse_aggregateData(rec).rois.AM;
    naturalMouse_aggregateData(recording_number).rois.LI = naturalMouse_aggregateData(rec).rois.LI;
    naturalMouse_aggregateData(recording_number).rois.S1 = naturalMouse_aggregateData(rec).rois.S1;
end
%get the coherencemap
naturalMouse_aggregateData(recording_number).coherence_map = coherenceCC;
naturalMouse_aggregateData(recording_number).reliability = reliability;

cd(first)

disp('Saving...')
save(['naturalMouse_aggregateData' date],'naturalMouse_aggregateData');
clear