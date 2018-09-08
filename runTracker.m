%% Set input file
video = '.\C216-ph_A1\C216-ph_A1_raw.avi';

%% Load configuration
load('config.mat');

%% Run MSER segmentation
calcMSER(video, config );

%% Run Kalman tracking
emkf2_MSER(video, config );

%% Export trajectories
[ trackData ] = exportKalmanTracks(video, config );