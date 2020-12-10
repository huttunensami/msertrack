%% Set input file
if ~exist('video', 'var')
   video = '.\videos\video.avi';
end

%% Load configuration
if ~exist('config', 'var')
   load('config.mat');
end

%% Run MSER segmentation
calcMSER(video, config );

%% Run Kalman tracking
emkf2_MSER(video, config );

%% Export trajectories
[ trackData ] = exportKalmanTracks(video, config );
