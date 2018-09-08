function [ maskPath, strOutput, strMATLAB, strKf, videoFile, strMSER ] = getTrackingPaths( config )
%GETTRACKINGPATHS Summary of this function goes here
%   Detailed explanation goes here

[strPath,strFile] = fileparts(config.sequence);

strMSER = ['_MSERd' num2str(config.Delta)];

if ~config.AreaLimit
   strMSER = [strMSER '_nomax'];
end

if config.MinDiversity < 0.9
   strMSER = [strMSER '_nested'];
end

if config.Preprocess
   strMSER = [strMSER '_filtered'];
end

maskPath = [strPath filesep strrep(strFile, ' ', '_') strMSER];

if config.freq==1
   strKf = ['_emkfq' num2str(config.q) 'r' num2str(config.r)];
else
   strKf = [num2str(freq) '_emkfq' num2str(config.q) 'r' num2str(config.r)];
end

if config.initNew
   strOutput = [strPath filesep strrep(strFile, ' ', '_') strKf strMSER];
else
   strOutput = [strPath filesep strrep(strFile, ' ', '_') strKf strMSER '_noinit'];
end

strMATLAB = [strOutput filesep 'MATLAB'];

%% AVI PARAMETERS
if config.initNew
   videoFile = fullfile(strOutput,[strrep(strFile, ' ', '_') strKf strMSER '.avi']);
else
   videoFile = fullfile(strOutput,[strrep(strFile, ' ', '_') strKf strMSER '_noinit.avi']);
end

end

