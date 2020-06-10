addpath(genpath('C:\Users\sit\Dropbox\CodeInBeta_Kevin\Other'))

%load your aggregate data

load(uigetfile)
recording_number = length(naturalMouse_aggregateData)+1;

parent_directory = uigetdir;
old = cd(parent_directory);
naturalMouse_aggregateData(recording_number).recording_name = parent_directory(4:end);

% choose the folder you're interested in
% add the entire folder to your path
matfiles = dir('**/*.mat');

for files = 1:length(matfiles)
    if strcmp(matfiles(files).name,'multi_movie_coherence_xcorr.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    elseif strcmp(matfiles(files).name,'VFS_maps.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    elseif strcmp(matfiles(files).name, 'reliability_map.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    elseif strcmp(matfiles(files).name, 'additional_maps.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    elseif strcmp(matfiles(files).name, 'magnitude_map.mat')
        load([matfiles(files).folder '\' matfiles(files).name])
    end
end

cd ..
load('VFS_meancat');

naturalMouse_aggregateData(recording_number).cortical_map = VFS_raw;
 
% run the registerer
[tform, Rfixed] = referenceMap_registration(VFS_raw);
naturalMouse_aggregateData(recording_number).transformation_parameters.tform = tform;
naturalMouse_aggregateData(recording_number).transformation_parameters.Rfixed = Rfixed;


%get the coherencemap
naturalMouse_aggregateData(recording_number).coherence_map = rot90(mean(coherenceCC, 3));
naturalMouse_aggregateData(recording_number).vertical_retinotopy = maps.VerticalRetinotopy;
naturalMouse_aggregateData(recording_number).horizontal_retinotopy = maps.HorizontalRetinotopy;
naturalMouse_aggregateData(recording_number).magnitude_map = rot90(magnitude_map);
naturalMouse_aggregateData(recording_number).reliability_map = reliability;

% Return
cd(old)

save(['naturalMouse_aggregateData_' date],'naturalMouse_aggregateData');