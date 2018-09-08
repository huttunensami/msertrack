%% Set input file
video = '.\videos\video.avi';

%% Load configuration
load('config.mat');

%% Run MSER segmentation
calcMSER(video, config );

%% Run Kalman tracking
emkf2_MSER(video, config );

%% Export trajectories
[ trackData ] = exportKalmanTracks(video, config );
