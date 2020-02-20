% coherentMotionAnalysis
clear;


%

data_type = 'new'; % or old
%% is something wrong with the sorting? idk why we have V1 resps that are sucky..
% load the data
fprintf('Choose your stimulus data \n')
[fn_stim,pn_stim] = uigetfile('.mat');
Stimdat = importdata([pn_stim fn_stim]);
data = matfile('DFF.mat');

repeats = Stimdat.repeats;
on_time    = Stimdat.on_time;
pre_time = 2;
off_time = Stimdat.off_time - pre_time;

num_movies = 3;


fs=10;
on_frames = on_time*fs;
pre_frames = pre_time*fs;
off_frames = off_time*fs;

mov_frames = on_frames + pre_frames + off_frames;

rep_frames = mov_frames * num_movies;

RespVec = zeros(size(data,'DFF',1),size(data,'DFF',2),on_frames,num_movies,repeats,'single');

for rep = 1:repeats
    for mov = 1:num_movies
        curr_frame = (rep-1)*rep_frames + (mov-1)*mov_frames + pre_frames;
        off_resp = mean(data.DFF(:,:,curr_frame+1 : curr_frame+off_frames),3);
        on_resp = data.DFF(:,:,curr_frame+off_frames+1 : curr_frame+off_frames+on_frames);
    RespVec(:,:,:,mov,rep) = on_resp - off_resp; % let's keep this off for now...
end
end

switch data_type
case 'old'
    movie1 = 1;
    movie2 = 2;
    movie3 = 3;
case 'new'
    movie1 = Stimdat.movieID(1);
    movie2 = Stimdat.movieID(2);
    movie3 = Stimdat.movieID(3);
end

load(...
    'C:\Users\sit\Dropbox\CodeInBeta_Kevin\Mouse Movies Stimulus\MouseMovies2014\mouse_movie_matfiles\combined_motion_information.mat');
magnitude_calc = double(magnitude_calculated);

mov_calc(:,1) = (resample(magnitude_calc(movie1,1:300),100,300)); %11 14 16?
mov_calc(:,2) = (resample(magnitude_calc(movie2,1:300),100,300));
mov_calc(:,3) = (resample(magnitude_calc(movie3,1:300),100,300));

% load('D:\Dropbox\CodeInBeta_Kevin\Mouse Movies Stimulus\MouseMovies2014\average_magnitude_information.mat');
% magnitude_meaned = double(magnitude_meaned);
% mov_uncorrected(:,1) = smooth(resample(magnitude_meaned(1:300,movie1),100,300)); %11 14 16?
% mov_uncorrected(:,2) = smooth(resample(magnitude_meaned(1:300,movie2),100,300));
% mov_uncorrected(:,3) = smooth(resample(magnitude_meaned(1:300,movie3),100,300));
% 
% 
% load('D:\Dropbox\CodeInBeta_Kevin\Mouse Movies Stimulus\MouseMovies2014\average_luminance_information.mat');
% mov_luminance(:,1) = smooth(resample(luminance_meaned(1:300,movie1),100,300)); %11 14 16?
% mov_luminance(:,2) = smooth(resample(luminance_meaned(1:300,movie2),100,300));
% mov_luminance(:,3) = smooth(resample(luminance_meaned(1:300,movie3),100,300));

trim = 10; % how much we trim to get rid of onset issues

coherenceCC = zeros(400,400,3);
coherenceCC_old = zeros(size(coherenceCC));
for mov = 1:num_movies
    for y = 1:size(RespVec,1)
        for x = 1:size(RespVec,2)
            disp([num2str(y) ' , ' num2str(x) ' , ' num2str(mov)])
            cell_trace = detrend(squeeze(mean(RespVec(y,x,trim+1:end,mov,:),5)));
            [temp,lags] = xcorr(cell_trace,mov_calc(1:end-trim,mov),5,'coeff');
            [coherenceCC(y,x,mov),idx] = max(temp);
            lag_map(y,x,mov) = lags(idx);
            coherenceCC_old(y,x,mov)  = corr(squeeze(mean(RespVec(y,x,trim:end-trim,mov,:),5)),mov_calc(trim:end-trim,mov));
        end
    end
end


save multi_movie_coherence_xcorr.mat coherenceCC coherenceCC_old
% 
% figure;
% subplot(3,2,1)
% imagesc(coherenceCC(:,:,1))
% title('New')
% axis square
% subplot(3,2,2)
% imagesc(coherenceOld(:,:,1))
% title('Old')
% axis square
% subplot(3,2,3)
% imagesc(coherenceCC(:,:,2))
% axis square
% subplot(3,2,4)
% imagesc(coherenceOld(:,:,2))
% axis square
% subplot(3,2,5)
% imagesc(coherenceCC(:,:,3))
% axis square
% subplot(3,2,6)
% imagesc(coherenceOld(:,:,3))
% axis square
% 
% 
% % 
% % %%
% % for ii = 1:3
% %     temp = rot90(motion_partial(:,:,ii));
% %     temp(VFS_boundaries) = max(temp(:))*1.1;
% %     motion_partial(:,:,ii) = temp;
% % end
% % 
% 
