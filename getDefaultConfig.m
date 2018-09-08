function [ config ] = getDefaultConfig( )
%GETDEFAULTCONFIG Summary of this function goes here
%   Detailed explanation goes here

   config.MinDiversity = 0.95;
   %config.MinDiversity = 0.2;
   config.MaxVariation = 0.3;
   %config.MaxVariation = 0.5;
   config.Delta = 1;
   config.DarkOnBright = 1;
   config.BrightOnDark= 1-config.DarkOnBright;
   %config.MaxArea = 0.0002;
   config.MaxArea = 0.5;
   config.MinArea = 0.00005;
   config.Verbose = 'Verbose';
   config.Dead = false;
   config.UseGreen = false;
   config.Preprocess = 0;
   if config.Preprocess
      %config.filter = @(I) uint8(bfilter2(double(I)./255,3,[1 3]) * 255);
      config.filter = @(I) medfilt2(I);
   else
      config.filter = @(I) I;
   end

   config.q = 1;
   config.r = 3;
   config.freq = 1;
   config.initNew = 0;
   
   config.saveAVI = 1;
   config.nested = 0;
   config.minArea = 5;
   %config.initFrame = 50;

   % Thresholds
   config.weightThreshold = 0.1;
   config.initThreshold = 2;
end

