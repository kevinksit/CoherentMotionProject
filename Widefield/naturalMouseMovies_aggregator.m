%load your aggregate data
load(uigetfile)
recording_number = length(naturalMouse_aggregateData)+1;

parent_directory = uigetdir;
cd(parent_directory)
naturalMouse_aggregateData(recording_number).recording_name = parent_directory(4:end);

% choose the folder you're interested in
% add the entire folder to your path
matfiles = dir('**/*.mat');


for files = 1:length(matfiles)
    if strcmp(matfiles(files).name,'multi_movie_coherence.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    elseif strcmp(matfiles(files).name,'VFS_meancat.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    end
end


naturalMouse_aggregateData(recording_number).cortical_map = VFS_raw;
naturalMouse_aggregateData(recording_number).VerticalRetinotopy = maps.VerticalRetinotopy;
naturalMouse_aggregateData(recording_number).HorizontalRetinotopy = maps.HorizontalRetinotopy;

% run the registerer
[tform, Rfixed] = referenceMap_registration(VFS_raw);
naturalMouse_aggregateData(recording_number).transformation_parameters.tform = tform;
naturalMouse_aggregateData(recording_number).transformation_parameters.Rfixed = Rfixed;


% get the ROIs
[v1_roi, lm_roi, al_roi, pm_roi] = HVAroiChooserGUI(VFS_raw);

naturalMouse_aggregateData(recording_number).rois.V1 = v1_roi;
naturalMouse_aggregateData(recording_number).rois.LM = lm_roi;
naturalMouse_aggregateData(recording_number).rois.AL = al_roi;
naturalMouse_aggregateData(recording_number).rois.PM = pm_roi;

%get the coherencemap
naturalMouse_aggregateData(recording_number).coherence_map_v2 = coherenceCC;
naturalMouse_aggregateData(recording_number).coherence_map = coherenceCC;

cd ..

save(['naturalMouse_aggregateData' date],'naturalMouse_aggregateData');
