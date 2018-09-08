function calcMSER(videoFile, config)

run('./vlfeat-0.9.21/toolbox/vl_setup.m');

[strPath,strFile] = fileparts(videoFile);
config.sequence = videoFile;

if nargin<2
   config.MinDiversity = 0.95;
   %config.MinDiversity = 0.2;
   config.MaxVariation = 0.3;
   %config.MaxVariation = 0.5;
   config.Delta = 1;
   config.DarkOnBright = 1;
   config.BrightOnDark= 1-config.DarkOnBright;
   config.MaxArea = 0.25;
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
else
  if ~config.Preprocess
      config.filter = @(I) I;
   end 
end

%% Open video
obj = VideoReader(videoFile);

if isempty(obj.NumberOfFrames)
   [~] = read(obj, inf);
end

% if obj.FrameRate ~= 0
%    frameRate = obj.FrameRate;
% else
%    frameRate = obj.NumberOfFrames/obj.Duration;
% end

numberOfFrames = obj.NumberOfFrames;
%numberOfFrames = 1;

disp(obj);

% Read first frame
I = read(obj, 1);

if config.UseGreen
   Ig = I(:,:,2);
elseif ndims(I)>2
   Ig = rgb2gray(I);
else
   Ig = I;
end

Ig = config.filter(Ig);

if config.AreaLimit
   config.MaxArea = 0.5;
   [r,~] = vl_mser(Ig, 'MinDiversity', config.MinDiversity, 'MaxVariation', config.MaxVariation, ...
   'Delta', config.Delta, 'DarkOnBright', config.DarkOnBright, ...
   'MaxArea', config.MaxArea, 'MinArea', config.MinArea, 'BrightOnDark', config.BrightOnDark, ...
   config.Verbose);

   M = zeros(size(Ig));
   
   for x=r'
      s = vl_erfill(Ig,x);
      M(s) = M(s) + 1;
   end
   
   mask = M>0;

   cells = regionprops(mask, 'Area');
   a = cat(1,cells.Area);
   
   meanArea = mean(a);
   imageArea = size(M,1) * size(M,2);
   
   config.MaxArea = 1.25 * meanArea/imageArea;
   config.MinArea = 0.25 * meanArea/imageArea;
end

% strOutput = [strPath filesep strrep(strFile, ' ', '_') '_MSERd' num2str(config.Delta) ];
% if ~config.AreaLimit
%    strOutput = [strOutput '_nomax'];
% end
% 
% if config.MinDiversity < 0.9
%    strOutput = [strOutput '_nested'];
% end
% 
% if config.Preprocess
%    strOutput = [strOutput '_filtered'];
% end
strOutput = getTrackingPaths( config );

%% Create output directories
strMATLAB = [strOutput filesep 'MATLAB'];
strLabel = [strOutput filesep 'label'];
strRaw = [strOutput filesep 'raw'];
[~,~] = mkdir(strOutput);
[~,~] = mkdir(strMATLAB);
[~,~] = mkdir(strLabel);
[~,~] = mkdir(strRaw);

%% Update symbolic links
if isunix
   cmd = ['ln -s "' strOutput '" ./running/'];
   system(cmd);
end

%% Visualization
fig = figure('Visible','on', 'MenuBar','none', 'ToolBar','none', 'Renderer', 'zbuffer');
imagesc(I);
axis off, colormap gray, axis image;
set(gcf,'Units','pixels','Position',[0 0 size(I,2) size(I,1)]);
set(gca,'Units','normalized','Position',[0 0 1 1]);

%set(fig,'PaperUnits','points');
%set(fig,'PaperSize',[size(I,2) size(I,1)]);
set(fig,'PaperPositionMode','auto');

imageFile = fullfile(strOutput, sprintf('MSER_%04d.png', 1));
%print(fig,'-dpng',imageFile,'-r96');

%% Initialize export
%cdata = hardcopy(fig, '-Dzbuffer', '-r0');
%imwrite(cdata,imageFile);
print(fig,'-dpng',imageFile,'-r0');
delete(imageFile);

%numberOfFrames = 1;
for i=1:numberOfFrames

   imageFile = fullfile(strRaw, sprintf('MSER_raw_%04d.png', i));
   matFile = fullfile(strMATLAB, sprintf('MSER_%04d.mat', i));
   
   maskMSER = fullfile(strOutput, sprintf('MSER_mask_%04d.png', i));
   labelMSER = fullfile(strLabel, sprintf('MSER_label_%04d.png', i));
   
%    if(exist(matFile, 'file'))
%       continue;
%    end
      
   I = read(obj, i);
   
   disp('**************************');
   
   tic
   
   if config.UseGreen
      Ig = I(:,:,2);
   elseif ndims(I)>2
      Ig = rgb2gray(I);
   else
      Ig = I;
   end
   
   Ig = config.filter(Ig);
   
%    [r,f] = vl_mser(Ig, 'MinDiversity', 0.7, 'MaxVariation',0.2, 'Delta',5, 'DarkOnBright',1, ...
%                        'MinArea',0.00005, 'BrightOnDark',0, 'Verbose');
%                        %Settings 1
   [r,f] = vl_mser(Ig, 'MinDiversity', config.MinDiversity, 'MaxVariation', config.MaxVariation, ...
                       'Delta', config.Delta, 'DarkOnBright', config.DarkOnBright, ...
                       'MaxArea', config.MaxArea, 'MinArea', config.MinArea, 'BrightOnDark', config.BrightOnDark, ...
                       config.Verbose);

   f = vl_ertr(f); %vl_plotframe(f);
   M = zeros(size(Ig));

   % Fill MSER regions
%    for x=r'
%       s = vl_erfill(Ig,x);
%       M(s) = M(s) + 1;
%    end

   borderMask = false(size(Ig));
   regionMask = borderMask;
   %tmpMask = regionMask;
   
   %se = strel('disk',1);
   for x=1:length(r)
      regionMask(:) = 0;
      s = vl_erfill(Ig,r(x));
      %M(s) = M(s) + 1;
      M(s) = x;
      
      regionMask(s) = 1;
      regionMask = bwmorph(regionMask, 'remove');
      borderMask(regionMask) = 1;
      
      %regionMask(:) = 0;
      %regionMask(s) = true;
      %tmpMask(:) = imerode(regionMask, se);
      
      %borderMask(:) = borderMask | (regionMask & ~tmpMask);
   end
   
   toc
   
   %% Save label and mask images
   mask = M>0;
   
   if ~any(mask)
      imwrite(Ig,imageFile);
      imwrite(mask, maskMSER, 'png');
      continue;
   end
   
   RGB = label2rgb(M, hsv(max(M(:))), 'k', 'shuffle');
   imwrite(mask, maskMSER, 'png');
   imwrite(RGB, labelMSER, 'png');
   
   %% Plot
   set(0, 'CurrentFigure', fig);
     
   %imagesc(Ig, [0 255]);
   I(:,:,1:3) = Ig(:,:,[1 1 1]);
   redMask = cat( 3, borderMask, false(size(borderMask)), false(size(borderMask)) );
   I(redMask) = 255;
   imagesc(I, [0 255]);
   
   hold on;
   % MSER
   %[~,h]=contour(M,(0:max(M(:)))+.5);
   %set(h,'Color','r','Linewidth',1);
   
   plot(f(1,:),f(2,:),'go', 'MarkerSize', 5);
   
   % Dead cells
   if config.Dead
      [stats, metric] = findDead( Ig );
      
      idxDead = metric > 0.8 | [stats.Area]' > config.MaxArea;
      deadCells = {stats.ConvexHull};
      deadCells = deadCells(idxDead);
      funDead = @(boundary)( plot(boundary(:,1), boundary(:,2), 'r', 'LineWidth', 2) );
      cellfun(funDead, deadCells);
   end
   
   hold off;
   drawnow;
   
   %% Save data
   axis off, colormap gray, axis image;
   set(gcf,'Units','pixels','Position',[0 0 size(Ig,2) size(Ig,1)]);
   set(gca,'Units','normalized','Position',[0 0 1 1]);

   %set(fig,'PaperUnits','points');
   %set(fig,'PaperSize',[size(Ig,2) size(Ig,1)]);
   set(fig,'PaperPositionMode','auto');
   
   print(fig,'-dpng',imageFile,'-r0'); 
   %cdata = hardcopy(fig, '-Dzbuffer', '-r0');
   %imwrite(cdata,imageFile);
   
   if config.Dead
      save(matFile, 'r', 'f', 'config', 'M', 'stats', 'metric', '-v7');
   else
      save(matFile, 'r', 'f', 'config', 'M', '-v7');
   end

end
close(fig);

%% Update symbolic links
if isunix
   paths = strfind(strOutput,filesep);
   seq = strOutput(paths(end)+1:end);
   cmd = ['unlink ./running/' seq];
   system(cmd);
   
   cmd = ['ln -s "' strOutput '" ./completed/'];
   system(cmd);
end
end
