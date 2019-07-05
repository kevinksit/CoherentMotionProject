function [] = signMapRegistration()
%% Sign Map Registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written 15Dec2017
% Last Updated:
%
% Uses a "master sign map" and its corresponding reference image to
% register the "master recording" to your current recording. Manually
% choose control points to line up the "master recording" to your current
% recording, then apply the transformation to the "master sign map" so that
% it aligns with your current recording window. The major output is a set
% of sign maps (processed and raw) which are correctly aligned to your
% current recording window, preventing the necessity of completing
% retinotopic mapping per recording.
%
% "Moving" refers to the master recording (which will be altered to match
% your current recording).
%
% "Fixed" refers to your current recording, which will not change.
%
%%% Necessary Subfunctions %%%
%
%%% Inputs %%%
%
%%% Outputs %%%
% registered_sign_maps.mat            VFS_raw and VFS_processed saved
%                                     after image registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data selection and import
%Use the reference images's to get the transformation matrix
disp('Choose your sign map...')
[fn1, pn1] = uigetfile('.mat');
disp('Choose the common coordinate map...') 
[fn2, pn2] = uigetfile('.mat');


common_coordinates = importdata([pn2 fn2]);
VFS_struct = importdata([pn1 fn1]);
VFS_p = rot90(VFS_struct.VFS_processed,-1);
VFS_r = rot90(VFS_struct.VFS_raw,-1);
VFS_b = rot90(VFS_struct.VFS_boundaries,-1);


%% Image preprocessing


% resize the ref_img's to 400x400 to match the sign maps
common_coordinates = imresize(common_coordinates, [400 400]);

%% Control point selection and curation
% from https://www.mathworks.com/help/images/control-point-registration.html
% control point selection and checking to ensure good registration
a=0;
while a == 0
    [movingPoints, fixedPoints] = cpselect(common_coordinates,abs(VFS_p),'Wait',true); % choose at least 3 noncollinear points
    
    tform = fitgeotrans(movingPoints,fixedPoints,'affine'); %use the control points to define a transformation matrix
    Rfixed = imref2d(size(common_coordinates)); % define a "world axis" for the images to reference onto
    
    registered_F0 = imwarp(VFS_p,tform,'OutputView',Rfixed); % warp the moving F0 to check performance
    
    % compare unregistered vs registered to check performance, if looks good, then continue
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1)
    imshowpair(common_coordinates,fixed_F0);
    title('Non-registered')
    subplot(1,2,2)
    imshowpair(registered_F0,fixed_F0);
    title('Registered')
    
    % prompt for good registration manual check
    good_registration = questdlg('Does the registration look good?','Registration Performance','Yes','No','Yes');
    
    close all;
    switch good_registration
        case 'Yes'
            a = 1; % break the loop and continue
        case 'No'
            disp('Probably an issue with your control points... let''s choose them again...');
    end
end

%% Applying transformations to sign maps
VFS_processed = rot90(imwarp(VFS_p_moving,tform,'OutputView',Rfixed)); %apply the transformation to VFS_processed
VFS_raw = rot90(imwarp(VFS_r_moving, tform,'OutputView',Rfixed)); % apply the transformation to VFS_raw
VFS_boundaries = rot90(imwarp(VFS_b_moving, tform, 'OutputView', Rfixed)); % apply the transformation to VFS_boundaries

%% Save warped sign maps
save registered_sign_maps.mat VFS_processed VFS_raw VFS_boundaries 
