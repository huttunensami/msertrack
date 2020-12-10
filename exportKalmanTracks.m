function [ trackData ] = exportKalmanTracks( seq, config )
%ANALYZETRACKS Summary of this function goes here
%   Detailed explanation goes here

if ispc
   saveAVI = config.saveAVI;
else
   saveAVI = 1;
end

%select = [];
%select = [212,119,944,477,201,120,1466,1234,1230,1464,1177,1005,876,761,1349,1042,866];

config.sequence = seq;
[ maskPath, strOutput, strMATLAB, strKf, videoFile ] = getTrackingPaths( config );
[~,strVideofile] = fileparts(videoFile);

%% Open video file
%startFrame = 1;
%endFrame = startFrame + 25;

%% Open ground truth sequence
obj = VideoReader(seq)
NumberOfFrames = get(obj, 'NumberOfFrames');
%NumberOfFrames = min(NumberOfFrames, endFrame-startFrame);

%% AVI PARAMETERS
if saveAVI
   fileName = fullfile([strOutput filesep strVideofile '_trajectory.avi']);
   
   writer = VideoWriter(fileName, 'Uncompressed AVI');
   writer.FrameRate = 1;
   open(writer);
   %mov = avifile(fileName, 'compression', 'None', 'fps', 1);
end

%% VISUALIZATION INIT
handlePosFig = figure('Name','Positions');
I = read(obj, 1);
imshow(I);

%% Allocate structures
results = dir(fullfile(strMATLAB, 'emkf_MSER_*.mat'));
if config.initNew
   matFile = fullfile(strMATLAB, results(end).name);
else
   matFile = fullfile(strMATLAB, results(1).name);
   %matFile = fullfile(strMATLAB, results(end).name);
end
matObj = matfile(matFile);
objects = matObj.objects;

if isempty(objects)
   nObjects = 1;
else
   nObjects = max(objects(:));
end

idx = length(results)*config.freq;
NumberOfFrames = min(NumberOfFrames, idx);

nData = floor( (NumberOfFrames-1) / config.freq ) + 1;
trackData(1:nObjects) = struct('x', nan(nData,1), 'y', nan(nData,1), ...
                               'u', nan(nData,1), 'v', nan(nData,1));
                            
mserCount = zeros(NumberOfFrames, 1);
objectCount = zeros(NumberOfFrames, 1);
trackCount = zeros(NumberOfFrames, 1);

for i=1:config.freq:NumberOfFrames
   disp('**************************************')
   disp([' FRAME: ' int2str(i)])
   
   matFile = fullfile(strMATLAB, sprintf('emkf_MSER_%04d.mat', i));
   matObj = matfile(matFile);

   maskFile = fullfile(maskPath, sprintf('MSER_mask_%04d.png', i));
   %gtFile = fullfile(gtPath, [filesep strrep(strFile, ' ', '_') '_tracking_gt_color' sprintf('%04d.png', i-1)]);

   mask = imread(maskFile);
   maskc = imclearborder(mask);
   maskCC = bwconncomp(mask, 8);
   maskcCC = bwconncomp(maskc, 8);
   
   mserCount(i) = maskCC.NumObjects;
   objectCount(i) = maskcCC.NumObjects;
   
   %% Open tracking data
   if exist(matFile)
      x = matObj.x;
      % targetBox = matObj.targetBox;
      objects = matObj.objects;
      objColors = matObj.objColors;

      xc = x(1,:);
      yc = x(3,:);
      u = x(2,:);
      v = x(4,:);

   %    xc = max(uint32(xc),1);
   %    yc = max(uint32(yc),1);

      NumberOfObjects = size(x,2);
   else
      % targetBox = [];
      objects = [];
      objColors = [];

      xc = [];
      yc = [];
      u = [];
      v = [];

   %    xc = max(uint32(xc),1);
   %    yc = max(uint32(yc),1);
      NumberOfObjects = 0;
   end

   trackCount(i) = NumberOfObjects;
   disp([' Number of tracked cells: ' num2str(NumberOfObjects)]);
   
   I = read(obj, i);
   imshow(I);
   hold on;

   %% Find matches
   for j=1:NumberOfObjects
      
      idx = (i-1)/config.freq + 1;
      
      if objects(j)>length(trackData)
         trackData(end:objects(j)) = struct('x', nan(nData,1), 'y', nan(nData,1), ...
                                            'u', nan(nData,1), 'v', nan(nData,1));
      end
      
      trackData(objects(j)).x(idx) = xc(j);
      trackData(objects(j)).y(idx) = yc(j);
      trackData(objects(j)).u(idx) = u(j);
      trackData(objects(j)).v(idx) = v(j);
      
      if saveAVI %&& ismember(objects(j), select)
         color = objColors( ( mod(objects(j)-1,size(objColors,1)) + 1 ),: );
         plot( xc(j), yc(j), '+', 'MarkerSize', 5, 'Color', color );
      
         xs = [trackData(objects(j)).x(1:idx)];
         ys = [trackData(objects(j)).y(1:idx)];

         line(xs, ys, 'Color', color, 'LineWidth', 2);
         
         tol = 15;
         top = min( size(mask,1) - tol, yc(j) - 3 );
         left = min( size(mask,2) - tol, xc(j) - 5 );
         
         text( left,top, ...
            num2str(objects(j)),...
            'Color', 'w', 'VerticalAlignment','bottom',...
            'FontSize', 8,...
            'BackgroundColor', [0.33 0.33 0.33],...
            'Margin', 1);
      end
         
   end
   
   hold off;
   %% SAVE AVI
   drawnow;
   
   if saveAVI
      %FRAME = getframe(handlePosFig);
      %writeVideo(writer, FRAME);
      set(handlePosFig,'PaperPositionMode','auto');
      %cdata = hardcopy(handlePosFig, '-Dzbuffer', '-r0');
      cdata = print('-RGBImage');
      writeVideo(writer, cdata);
   end   
end

   %% MATLAB data file
   %if initNew
      %matFile = fullfile(strOutput,[strrep(strFile, ' ', '_') strKf '_MSERd' num2str(delta) '_trajectory.mat']);
      matFile = [strOutput filesep strVideofile '_trajectory.mat'];
   %else
   %   matFile = fullfile(strOutput,[strrep(strFile, ' ', '_') strKf '_MSERd' num2str(delta) '_noinit_trajectory.mat']);
   %end
   save(matFile, 'trackData', 'config', 'mserCount', 'objectCount', 'trackCount', '-v7');

   %% Excel
   if ispc
      NumberOfObjects = length(trackData);
      
      excelFile = [strOutput filesep strVideofile '_trajectory.xlsx'];
      %excelFile = [strOutput filesep 'trajectory.xlsx'];
%       if initNew
%          excelFile = fullfile(strOutput,[strrep(strFile, ' ', '_') strKf '_MSERd' num2str(delta) '_trajectory.xlsx']);
%       else
%          excelFile = fullfile(strOutput,[strrep(strFile, ' ', '_') strKf '_MSERd' num2str(delta) '_noinit_trajectory.xlsx']);
%       end
      delete(excelFile);

      f = {'Frame'};
      c = {'Count'};
      o = { 'Tracks' };
      fs = 1:NumberOfFrames;
      os = nan(7*NumberOfObjects, 1, 'single');
      os(1:7:end) = 1:NumberOfObjects;
      osc = num2cell(os);
      osc(isnan(os)) ={''};

      count = {'MSER','active','track'}';
      header = {'x','y','u','v','d','dt',' '}';
      header = repmat(header, [NumberOfObjects 1]);

      xlswrite(excelFile, f, 1, 'C1');
      xlswrite(excelFile, fs, 1, 'C2');

      mserCount = mserCount(1:config.freq:end)';
      objectCount = objectCount(1:config.freq:end)';
      trackCount = trackCount(1:config.freq:end)';
      xlswrite(excelFile, c, 1, 'A2');
      xlswrite(excelFile, count, 1, 'B3');
      xlswrite(excelFile, mserCount, 1, 'C3');
      xlswrite(excelFile, objectCount, 1, 'C4');
      xlswrite(excelFile, trackCount, 1, 'C5');

      xlswrite(excelFile, o, 1, 'A6');
      xlswrite(excelFile, osc, 1, 'A7');

      xlswrite(excelFile, header, 1, 'B7');
      data = nan(nData, 7*NumberOfObjects, 'single');

      for k = 1:NumberOfObjects
         xs = [trackData(k).x];
         ys = [trackData(k).y];
         us = [trackData(k).u];
         vs = [trackData(k).v];
         
         valid = ~isnan(xs);
         d = nan(length(xs),1);
         dt = d;
         xv = xs(valid);
         yv = ys(valid);
         valid(1:find(~isnan(xs),1)) = false;
         
         d(valid) = sqrt( ( xv(2:end)-xv(1:end-1)).^2 + ( yv(2:end)-yv(1:end-1)).^2);
         dt(valid) = cumsum(d(valid));

         idx = 1 + (k-1) * 7;

         data(1:nData, idx) = xs;
         data(1:nData, idx + 1) = ys;
         data(1:nData, idx + 2) = us;
         data(1:nData, idx + 3) = vs;
         data(1:nData, idx + 4) = d;
         data(1:nData, idx + 5) = dt;
      end

      data = data';

      datac = num2cell(data);
      datac(isnan(data)) ={''};
      xlswrite(excelFile, datac, 1, 'C7');
   end

end

