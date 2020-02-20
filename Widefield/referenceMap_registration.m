function [tform,Rfixed] = referenceMap_registration(VFS)

% recording registration maker doer, thing...
%% first deterimining all the CoMs of the visual areas using the ccf stuff
% outline = importdata('E:\_Code\CoherentMotionProject\Widefield\reference_map.mat');
outline = imread('E:\_Code\CoherentMotionProject\Widefield\ABI_ctx_outline.tif');
%%% REGISTERER %%%

%% Control point selection and curation
% from https://www.mathworks.com/help/images/control-point-registration.html
% control point selection and checking to ensure good registration
% just cause i'm too lazy to change the actual variable names oops. The
% *128 is to scale it into a "seeable range" by image
moving_F0 = ind2rgb(gray2ind(mat2gray(VFS),64),jet); 
fixed_F0 = outline;

a=0;
while a == 0
    [movingPoints, fixedPoints] = cpselect(moving_F0,fixed_F0,'Wait',true); % choose at least 3 noncollinear points
    
    tform = fitgeotrans(movingPoints,fixedPoints,'affine'); %use the control points to define a transformation matrix
    Rfixed = imref2d(size(fixed_F0)); % define a "world axis" for the images to reference onto
    
    registered_F0 = imwarp(moving_F0,tform,'OutputView',Rfixed); % warp the moving F0 to check performance
    
    % compare unregistered vs registered to check performance, if looks good, then continue
    figure('units','normalized','outerposition',[0 0 1 1]);
%     subplot(1,2,1)
%     imshowpair(moving_F0,fixed_F0);
%     title('Non-registered')
%     subplot(1,2,2)
    registered_F0(~outline) = max(registered_F0(:))*1.1;
    imagesc(registered_F0)
    axis image
%     imshowpair(registered_F0,fixed_F0);
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

save transformation_parameters.mat tform Rfixed