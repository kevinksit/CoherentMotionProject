
files = dir('*_nat.mat');
% If you need to read...
v = VideoReader('dierction0.mp4');
f = 1;

raw_movie = importdata('direction0.mp4');

for ii = 1:size(movie, 4)
    gray_movie(:, :, ii) = imresize(rgb2gray(movie(:, :, :, ii)), 0.5);
end

opflow = opticalFlowHS();
    


for mov = 1:8;%length(files)
    
   disp(num2str(mov))
%movie = importdata(files(mov).name);
movie = squeeze(rdk_movie(mov, :, :, :));
for ii = 1:size(movie,3)
disp(ii)
flo = opflow.estimateFlow(movie(:,:,ii));
% 
% mag(ii) = sum(sum(flo.Magnitude));
% vx(:,:,ii) = ((flo.Vx));
% vy(:,:,ii) = ((flo.Vy));
% 
% %frame_motion(:,:,ii) = sqrt((vx(:,:,ii).^2 + mean(vy(:,:,ii).^2)));
% 
 frame_motion(:, :, ii) = flo.Magnitude;
end

sum_vx = sum(sum(vx,1),2);
sum_vy = sum(sum(vy,1),2);

calc = sqrt(sum_vx.^2 + sum_vy.^2);

magnitude_calculated(mov,:) = calc;    

mean_screen_motion(:,:,mov) = mean(frame_motion,3);

reset(opflow);
end

save new_mean_screen_motion.mat mean_screen_motion
save combined_motion_information.mat magnitude_calculated mean_screen_motion
