% Motion Analysis

v = VideoReader('RDKSTIM.mp4');
ct=1;
while hasFrame(v)
    video(:,:,ct) = double(rgb2gray(readFrame(v)));
    ct=ct+1;
end

opticFlow = opticalFlowHS;

for ii = 1:size(movie,3);
    flow(ii) = estimateFlow(opticFlow,movie(:,:,ii));
end

for ii = 1:length(flow)
    RDK_magnitude(ii) = mean(flow(ii).Magnitude(:));
    
    % considering 1 deg = 0.017rad, i think rounding to the hundredth's is
    % good
    ori = round(flow(ii).Orientation,);
  RDK_coherence(ii) = mean(mean(ori == mode(ori(:))));
end


for ii = 1:449
    imagesc(flow(ii).Orientation)
    pause
end
%%%%%%%%%%%55

movie_idx = dir('*_nat.mpg');
for ii = 1:length(movie_idx)
    disp(num2str(ii))
    opticFlow = opticalFlowHS;
    
    v = VideoReader(movie_idx(ii).name);
    ct=1;
    while hasFrame(v)
       temp(ct) = estimateFlow(opticFlow,double(rgb2gray(readFrame(v))));
        ct=ct+1;
    end
    motion_magnitude{ii} = cat(3,temp.Magnitude);
    motion_orientation{ii} = cat(3,temp.Orientation);
    % store variables
    
    clear opticFlow
end


for ii = 1:size

avg_motion = mean(cat(3,flow.Magnitude),3);

imagesc(avg_motion)
axis image
box off
axis off