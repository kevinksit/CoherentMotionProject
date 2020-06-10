%addpath(genpath('C:\Users\sit\Dropbox\CodeInBeta_Kevin\Other'))

%load your aggregate data

load(uigetfile)
recording_number = length(coherentMotion_aggregateData)+1;

parent_directory = uigetdir;
old = cd(parent_directory);
coherentMotion_aggregateData(recording_number).recording_name = parent_directory(4:end);

% choose the folder you're interested in
% add the entire folder to your path
matfiles = dir('**/*.mat');

for files = 1:length(matfiles)
    if strcmp(matfiles(files).name,'coherence_maps.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    elseif strcmp(matfiles(files).name,'VFS_maps.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    elseif strcmp(matfiles(files).name, 'coherence_curve.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    elseif strcmp(matfiles(files).name, 'additional_maps.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    end
end

coherentMotion_aggregateData(recording_number).cortical_map = VFS_raw;
 
% run the registerer
[tform, Rfixed] = referenceMap_registration(VFS_raw);
coherentMotion_aggregateData(recording_number).transformation_parameters.tform = tform;
coherentMotion_aggregateData(recording_number).transformation_parameters.Rfixed = Rfixed;


% get the ROIs-
disp('Choose V1, LM, AL, PM...')
[v1_roi, lm_roi, al_roi, pm_roi] = HVAroiChooserGUI(VFS_raw);
disp('Choose AM, RL, LI, S1...')
[am_roi, rl_roi, li_roi, s1_roi] = HVAroiChooserGUI(VFS_raw);
coherentMotion_aggregateData(recording_number).rois.V1 = v1_roi;
coherentMotion_aggregateData(recording_number).rois.LM = lm_roi;
coherentMotion_aggregateData(recording_number).rois.AL = al_roi;
coherentMotion_aggregateData(recording_number).rois.PM = pm_roi;
coherentMotion_aggregateData(recording_number).rois.AM = am_roi;
coherentMotion_aggregateData(recording_number).rois.RL = rl_roi;
coherentMotion_aggregateData(recording_number).rois.LI = li_roi;
coherentMotion_aggregateData(recording_number).rois.S1 = s1_roi;

%get the coherencemap
coherentMotion_aggregateData(recording_number).coherence_map = coherenceCC;

coherentMotion_aggregateData(recording_number).vertical_retinotopy = maps.VerticalRetinotopy;
coherentMotion_aggregateData(recording_number).horizontal_retinotopy = maps.HorizontalRetinotopy;
% Return
cd(old)

save(['coherentMotion_aggregateData_' date],'coherentMotion_aggregateData');
