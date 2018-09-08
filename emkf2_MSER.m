function [] = emkf2_MSER(sequence,config)

%sequence = '/research/imag/development/biocenter/Videos/woundhealingcontrol/woundhealingcontrol.avi';
%freq = 3;
%[strPath,strFile] = fileparts(sequence);

%% Parameters
if nargin<2
   config = getDefaultConfig();
end

config.sequence = sequence;

% KALMAN FILTERING PARAMETERS
InitialStateVar = 6;
SystemVar = config.q;
r = [config.r 0;
     0 config.r];

%% Open files
[ maskPath, strOutput, strMATLAB, ~, fileName ] = getTrackingPaths( config );

%% Old
% maskPath = [strPath filesep strrep(strFile, ' ', '_') '_MSERd' num2str(config.delta)];
% 
% if config.nomax
%    maskPath = [maskPath '_nomax'];
% end
% 
% if config.nested
%    maskPath = [maskPath '_nested'];
% end
% 
% if config.filtered
%    maskPath = [maskPath '_filtered'];
% end
% 
% if config.freq==1
%    strKf = ['_emkfq' num2str(SystemVar) 'r' num2str(r(1,1))];
% else
%    strKf = [num2str(freq) '_emkfq' num2str(SystemVar) 'r' num2str(r(1,1))];
% end
% 
% if config.initNew
%    strOutput = [strPath filesep strrep(strFile, ' ', '_') strKf '_MSERd' num2str(delta)];
% else
%    strOutput = [strPath filesep strrep(strFile, ' ', '_') strKf '_MSERd' num2str(delta) '_noinit'];
% end
% 
% if config.nested
%    strOutput = [strOutput '_nested'];
% end
% 
% if config.filtered
%    strOutput = [strOutput '_filtered'];
% end
% 
% strMATLAB = [strOutput filesep 'MATLAB'];

%% Create directories
[~,~] = mkdir(strOutput);
[~,~] = mkdir(strMATLAB);

%% Update symbolic links
if isunix
   cmd = ['ln -s "' strOutput '" /research/biovis/development/biocenter/tracking/running/'];
   system(cmd);
end

%% AVI PARAMETERS
if config.saveAVI
   writer = VideoWriter(fileName, 'Uncompressed AVI');
   writer.FrameRate = 1;
   open(writer);
   %mov = avifile(fileName, 'compression', 'None', 'fps', 1);
end

%% LOAD IMAGE SEQUENCE
startFrame = 1;
endFrame = startFrame+10000;
endFrame = endFrame*config.freq;

%% Old objects
x_old = [];
P_old = [];
Q_old = [];
R_old = [];
objects_old = [];
targetBox_old = [];

%% Removed objects
x_remove = [];
P_remove = [];
Q_remove = [];
R_remove = [];
objects_remove = [];
targetBox_remove = [];
ws_remove = [];

%% Open video
obj = VideoReader(sequence)

IMG = read(obj, 1);

NumberOfFrames = get(obj, 'NumberOfFrames');
NumberOfFrames = min(NumberOfFrames, endFrame-startFrame);

height=size(IMG,1);
width=size(IMG,2);

disp(['Number of frames: ' int2str(NumberOfFrames)])

%% VISUALIZATION INITIALIZATION
scale = 1;

%handleBBFig = figure('Name','Detections');
%imshow(IMG(:,:,:,1),'Border','Tight');
%set(handleBBFig,'Position',[25 655 scale*width scale*height]);

handlePosFig = figure('Name','Positions', 'Renderer', 'zbuffer');
%posFig = get(handlePosFig,'Position');
%set(handlePosFig,'Position',[posFig(1) posFig(2) scale*2*width scale*height]);
%iptwindowalign(handleBBFig, 'top', handlePosFig, 'top');
%iptwindowalign(handleBBFig, 'right', handlePosFig, 'left');
if width<1.75*height
   imshow([IMG(:,:,:,1) IMG(:,:,:,1)],'Border','Tight');
else
   imshow([IMG(:,:,:,1); IMG(:,:,:,1)],'Border','Tight');
end
%handlePos = gca;

if width<1.75*height
   handlePos = subplot(1,2,1);
   set(gca,'Position',[0 0 0.5 1]);
else
   handlePos = subplot(2,1,1);
   set(gca,'Position',[0 0.5 1 0.5]);
end

%handleWeightFig = figure('Name','Weights');
%weightFig = get(handlePosFig,'Position');
%imshow(IMG(:,:,:,1),'Border','Tight');
%handleWeight = gca;

if width<1.75*height
   handleWeight = subplot(1,2,2);
   set(gca,'Position',[0.5 0 0.5 1]);
else
   handleWeight = subplot(2,1,2);
   set(gca,'Position',[0 0 1 0.5]);
end

%figure(handleBBFig);

%% SELECT OBJECTS
disp('--------------------------------------------')
disp('SELECT OBJECTS')
disp('--------------------------------------------')


%% Detect objects
se = strel('disk',1);

mask = imread([maskPath filesep 'MSER_mask_0001.png']);
%mask = computeRoundness(mask);

mask = imclearborder(mask);
mask = imfill(mask, 'holes');
mask = imopen(mask,se);
mask = bwareaopen(mask, config.minArea);
[ ~,~,bboxs,nIter ] = CreateMeasurements( mask );

%% Process boxes
%figure(handleBBFig);
%set(0,'CurrentFigure',handleBBFig);
set(0,'CurrentFigure',handlePosFig);

subplot(handleWeight);
imshow(mask,'Border','tight');
   
subplot(handlePos);
imshow(IMG,'Border','tight');

subplot(handleWeight);
hold on
[ newBoxes, numberOfMatches, covar ] = process_bboxes( bboxs );
%% Only reliable observations are considered
indexes = numberOfMatches > config.initThreshold;
newBoxes = newBoxes(indexes,:);
covar = covar(:,:,indexes);
hold off

NumberOfObjects = size(newBoxes, 1);
objects = 1:NumberOfObjects;
targetBox = repmat(struct('xmin',[],'xmax',[],'ymin',[],'ymax',[]),1,NumberOfObjects);

%% CREATE COLOR LABELS
RGB = label2rgb(1:NumberOfObjects, 'hsv','r','shuffle');
RGB=double(RGB)./255;
objColors = reshape(RGB,[NumberOfObjects 3]);

%% INITIALIZE STATES
x = zeros(4,NumberOfObjects);
R = zeros(2,2,NumberOfObjects);

subplot(handlePos);
hold on
for k=1:NumberOfObjects
   %% Create new object
   w = newBoxes(k,3);
   h = newBoxes(k,4);
   xc = newBoxes(k,1) + w/2;
   yc = newBoxes(k,2) + h/2;
   targetBox(k) = get_bbox(xc,yc,w,h);
   
   disp(['Object #' int2str(k) ': xc =' int2str(xc) ': yc =' int2str(yc) ])
   x(:,k) = [xc;0;yc;0];
   %R(:,:,k) = r;
   var = covar(:,:,k);
   varx = max( r(1,1), min(w/4,var(1,1)) );
   vary = max( r(2,2), min(h/4,var(2,2)) );
   
   R(:,:,k) = diag([varx vary]);
   
   rectBox = get_rect(targetBox(k));
   color = objColors(objects(k),:);
   rectangle('Position',rectBox, 'LineWidth',1,'LineStyle','-','EdgeColor',color);
   text(targetBox(k).xmin,targetBox(k).ymin,num2str(objects(k)),'Color','w','VerticalAlignment','bottom','FontSize',6);
end
hold off

drawnow;

if config.saveAVI
   %FRAME = getframe(handlePosFig);
   %mov = addframe(mov,FRAME);
   set(handlePosFig,'PaperPositionMode','auto');
   %cdata = hardcopy(handlePosFig, '-Dzbuffer', '-r0');
   cdata = print('-RGBImage');
   writeVideo(writer, cdata);
end

disp(['Number of objects: ' int2str(NumberOfObjects)])
objectsTotal = NumberOfObjects;
%% INITIALIZE FILTER
[Q,P,pi,F,H] = emkf2_init(repmat(InitialStateVar,1,NumberOfObjects),repmat(SystemVar,1,NumberOfObjects));
pi = ones(NumberOfObjects,1);%/NumberOfObjects;

disp('--------------------------------------------')
disp('INIT FILTERING')
disp('--------------------------------------------')
%% GO THROUGH ALL FRAMES
%pause;

matFile = fullfile(strMATLAB, sprintf('emkf_MSER_%04d.mat', startFrame));
save(matFile, 'x', 'R', 'targetBox', 'objColors', 'objects', 'config');

for i=(startFrame+config.freq):config.freq:startFrame+NumberOfFrames-1
   disp('****************************************************************************************')
   disp([' FRAME: ' int2str(i)])
   
   clear IMG;
   
%    if i==initFrame
%       initNew = 1;
%    else
%       initNew = 0;      
%    end
      
   IMG = read(obj, i);
   
   %% Predict position
   
   [xp,Pp] = emkf2_predict(x,F,P,Q);
   
   %    for k=1:NumberOfObjects
   %       targetBox(k) = update_bbox(targetBox(k),xp(1,k),xp(3,k));
   %    end
   
   arrayfun(@update_bbox,targetBox,xp(1,:),xp(3,:));
   
   %% GET MEASUREMENTS
   disp('GET MEASUREMENTS')
   
   mask = imread([maskPath filesep sprintf('MSER_mask_%04d.png', i)]);
   %mask = computeRoundness(mask);
   
   mask = imclearborder(mask);
   mask = imfill(mask, 'holes');
   mask = imopen(mask,se);
   mask = bwareaopen(mask, config.minArea);
   
   [ xs,ys,bboxs,boxInit,ccInit ] = CreateMeasurements( mask );
   NumberOfDetections = length(xs);
   
   z=[xs ys]';
   
   disp(['Number of Measurements: ' int2str(size(z,2))])
   
   %% UPDATE FILTER
   disp('--------------------------------------------')
   disp('UPDATE FILTER')
   
   clear ws;
   %[mu,C,ws] = emkf2_getWeights(xp,Pp,H,z,pi,R);
   [x(:),P(:),pi,ws] = emkf2_update(xp,Pp,H,z,pi,R,config.weightThreshold);
   
   rectScore = zeros(NumberOfDetections, NumberOfObjects);
   bboxAreas = bboxs(:,3).*bboxs(:,4);
   
   allMatches = false(NumberOfDetections,1);
   alpha = 0.2;
   foundObjects = false(1,length(objects));
   
   %% Update objects
   for k=1:NumberOfObjects
      
      targetBox(k) = update_bbox(targetBox(k),x(1,k),x(3,k));
      
      %% Use weights
      rectBox = get_rect(targetBox(k));
      rectArea = repmat(rectBox(3)*rectBox(4),NumberOfDetections,1);
      areaUnion = rectunion(rectBox,bboxs)';
      areaIntersect = rectint(rectBox,bboxs)';
      rectScore(:,k) = areaIntersect./areaUnion;
      
      %  Find matches that are inside an object box
      intersects = (areaIntersect==bboxAreas) | (areaIntersect==rectArea);
            
      matches = ws(:,k)>config.weightThreshold;
      weights = ws(matches,k);
      matchBoxes = bboxs(matches,:);
      allMatches = allMatches | matches | intersects | rectScore(:,k)>0.2;
      NumberOfMatches = size(matchBoxes,1);
            
      % Adjust scale
      if NumberOfMatches>0
         dims = weights'*matchBoxes;
         dims = dims./sum(weights);
         
         w = alpha*dims(3)+(1-alpha)*rectBox(3);
         h = alpha*dims(4)+(1-alpha)*rectBox(4);
         targetBox(k) = get_bbox(x(1,k),x(3,k),w,h);
         
         foundObjects(k) = 1;
      end
      %
      %% Use rectangles
%       weightMatches = ws(:,k)>weightThreshold ;
%       
%       rectBox = get_rect(targetBox(k));
%       rectArea = repmat(rectBox(3)*rectBox(4),NumberOfDetections,1);
%       areaUnion = rectunion(rectBox,bboxs)';
%       areaIntersect = rectint(rectBox,bboxs)';
%       rectScore(:,k) = areaIntersect./areaUnion;
%       
%       %Find matches that are inside an object box
%       bboxAreas = bboxs(:,3).*bboxs(:,4);
%       intersects = (areaIntersect==bboxAreas) | (areaIntersect==rectArea);
%       
%       matches =  weightMatches | intersects | rectScore(:,k)>0.2;
%       allMatches = allMatches | matches;
%       matchBoxes = bboxs(matches,:);
%       NumberOfMatches = size(matchBoxes,1);
%       matchScores = rectScore(matches,k);
%       
%       % Adjust scale
%       if NumberOfMatches>0 && sum(matchScores)>0
%          dimNew = matchScores'*matchBoxes(:,3:4)/sum(matchScores);
%          
%          w = alpha*dimNew(1)+(1-alpha)*rectBox(3);
%          h = alpha*dimNew(2)+(1-alpha)*rectBox(4);
%          targetBox(k) = get_bbox(x(1,k),x(3,k),w,h);
%          
%          if isnan(targetBox(k).xmin)
%             disp('Error');
%          end
%          
%          foundObjects(k) = 1;
%       end
      
      %% Display object status
      disp(['Object #' int2str(k) ' State'])
      disp(['Pi: ' num2str(pi(k)*100) ' %'])
      disp('Predicted State:  Corrected State:')
      disp([xp(:,k)  x(:,k)])
      disp('-----------------')
   end
   
   %% Remove objects
   removeObjects = ~foundObjects;
   if i==startFrame+config.freq
      x_old = x(:,removeObjects);
      P_old = P(:,:,removeObjects);
      Q_old = Q(:,:,removeObjects);
      R_old = R(:,:,removeObjects);
      objects_old = objects(removeObjects);
      targetBox_old = targetBox(removeObjects);
      ws_old = ws(:,removeObjects);
   else
      x_remove = x(:,removeObjects);
      P_remove = P(:,:,removeObjects);
      Q_remove = Q(:,:,removeObjects);
      R_remove = R(:,:,removeObjects);
      objects_remove = objects(removeObjects);
      targetBox_remove = targetBox(removeObjects);
      ws_remove = ws(:,removeObjects);
   end
   
   x = x(:,foundObjects);
   P = P(:,:,foundObjects);
   Q = Q(:,:,foundObjects);
   R = R(:,:,foundObjects);
   pi = pi(foundObjects);
   objects = objects(foundObjects);
   NumberOfObjects = length(objects);
   targetBox = targetBox(foundObjects);
   ws = ws(:,foundObjects);
   
   %% Show images
   %figure(handlePosFig)
   set(0,'CurrentFigure',handlePosFig);
   subplot(handleWeight);
   imshow(mask,'Border','tight');
   
   subplot(handlePos);
   imshow(IMG,'Border','tight');
   
   %% Group unmatched boxes
   %figure(handlePosFig)
   set(0,'CurrentFigure',handlePosFig);
   %imshow(IMG,'Border','tight');
   
   %% Find not allocated cells
   L = uint32(labelmatrix(ccInit));
   
   idxz = uint32( sub2ind(size(L),uint32(ys(allMatches)),uint32(xs(allMatches))) );
   Lm = unique(L(idxz));
   Li = ~ismember(L,Lm);
   Li = Li & mask;
   Li = imclearborder(Li);
   
   stats = regionprops(Li,'basic');
   notMatched = cat(1,stats.BoundingBox);
   
   subplot(handleWeight);
   hold on
   %notMatched = bboxs(~allMatches&boxInit,:);
   %[ newBoxes, numberOfMatches, covar ] = process_bboxes( notMatched );
   newBoxes = notMatched;
   
   for k = 1:size(newBoxes,1)
      rectangle('Position', newBoxes(k,:), 'LineWidth',1,'LineStyle','-','EdgeColor','r');
   end
   hold off;
   
   %% Only reliable observations are considered
   %indexes = numberOfMatches>0;
   %newBoxAreas = newBoxes(:,3).*newBoxes(:,4);
   
   %indexes = indexes & newBoxAreas>minArea;
   %newBoxes = newBoxes(indexes,:);
   
   %covar = covar(:,:,indexes);
   %numberOfMatches = numberOfMatches(indexes);
   
   numberOfNewBoxes = size(newBoxes, 1);
   
   %% Match or add objects
   matched = false(1,length(objects_old));
   
   %% Predict position
   %[x_old,P_old] = emkf2_predict(x_old,F,P_old,Q_old);
   
   for k=1:numberOfNewBoxes
      
      mw = newBoxes(k,3);
      mh = newBoxes(k,4);
      mx = newBoxes(k,1) + mw/2;
      my = newBoxes(k,2) + mh/2;
      
      % Expand rectangle
      xmin = newBoxes(k,1) - newBoxes(k,3);
      xmax = newBoxes(k,1) + 2*newBoxes(k,3);
      ymin = newBoxes(k,2) - newBoxes(k,4);
      ymax = newBoxes(k,2) + 2*newBoxes(k,4);
%       xmin = newBoxes(k,1);
%       xmax = newBoxes(k,1) + newBoxes(k,3);
%       ymin = newBoxes(k,2);
%       ymax = newBoxes(k,2) + newBoxes(k,4);
      
      %% Check if an old object is matched
      pointInRect = 0;
      numberOfOldObjects = size(x_old,2);
      
      %% Go through all old objects
      for l=1:numberOfOldObjects
         
         %% Go to the next one if the current object has been already matched
         if matched(l)
            continue;
         end
         
         %% Make matching area a bit larger
         bw = targetBox_old(l).xmax-targetBox_old(l).xmin;
         bh = targetBox_old(l).ymax-targetBox_old(l).ymin;
                
         ox = x_old(1,l);
         oy = x_old(3,l);
         
         pointInRect = ox>xmin && ox<xmax && oy>ymin && oy<ymax;
         %pointInRect = 0;
         if pointInRect
            %% Resume old state
            %             Pn = P_old(:,:,l);
            %             Qn = Q_old(:,:,l);
            %             Rn = R_old(:,:,l);
            %             idx = objects_old(l);            
            matched(l) = 1;
            
            %% Update position and dimensions
            w = alpha*mw + (1-alpha)*bw;
            h = alpha*mh + (1-alpha)*bh;
            %xc = alpha*mx + (1-alpha)*x_old(1,l);
            %yc = alpha*my + (1-alpha)*x_old(3,l);
            xc = mx;
            yc = my;
            
            vx = xc - x_old(1,l);
            vy = yc - x_old(3,l);
            
            x_old(:,l) = [xc;vx;yc;vy];
            
            %[x_old(:,l),P_old(:,:,l),~,~] = emkf2_update(x_old(:,l),P_old(:,:,l),H,[mx;my],1,R,1,0);
            
            targetBox_old(l) = get_bbox(x_old(1,l),x_old(3,l),w,h);

            %% Old
            %xn = [xc;vx;yc;vy];
            %             x = cat(2, x, xn);
            %             P = cat(3, P, Pn);
            %             Q = cat(3, Q, Qn);
            %             R = cat(3, R, Rn);
            %             NumberOfObjects = NumberOfObjects+1;
            %             if isempty(pi)
            %                pi = 1;
            %             else
            %                %pi = cat(1,pi,1/NumberOfObjects);
            %                pi = (NumberOfObjects-1)/NumberOfObjects*pi;
            %                pi = cat(1,pi,1/NumberOfObjects);
            %             end
            
            %targetBox(NumberOfObjects) = get_bbox(xc,yc,w,h);
            %objects = cat(2,objects,idx);
            
            break;
         end
      end
      
      %% Check whether there are enough detections or not
      if ~config.initNew || pointInRect
         continue;
      end
      
      %% No match found, initialize entirely new object
      [Qn,Pn,pin,Fn,Hn] = emkf2_init(InitialStateVar,SystemVar);
      vx = 0;
      vy = 0;
      %Rn = r;
      var = covar(:,:,k);
      %       varx = min(mw/4,var(1,1));
      %       vary = min(mh/4,var(2,2));
      varx = max( r(1,1), min(mw/4,var(1,1)) );
      vary = max( r(2,2), min(mh/4,var(2,2)) );
      
      Rn = diag([varx vary]);
      %Rn = covar(:,:,k);
      objectsTotal = objectsTotal + 1;
      idx = objectsTotal;
      disp(['New object: ' 'ID: ' num2str(idx) ' R: ']);
      disp(Rn);
      
      %% Create new object
      w = newBoxes(k,3);
      h = newBoxes(k,4);
      xc = newBoxes(k,1) + w/2;
      yc = newBoxes(k,2) + h/2;
      newTargetBox = get_bbox(xc,yc,w,h);
      xn = [xc;vx;yc;vy];
      
      x = cat(2, x, xn);
      P = cat(3, P, Pn);
      Q = cat(3, Q, Qn);
      R = cat(3, R, Rn);
      NumberOfObjects = NumberOfObjects+1;
      if isempty(pi)
         pi = 1;
      else
         pi = (NumberOfObjects-1)/NumberOfObjects*pi;
         pi = cat(1,pi,1/NumberOfObjects);
      end
      targetBox(NumberOfObjects) = newTargetBox;
      objects = cat(2,objects,idx);
   end
   
   %% Add matched objects and remove not matched
   if ~isempty(matched)
      prevNumberOfObjects = length(objects);
      x = cat(2, x, x_old(:,matched));
      P = cat(3, P, P_old(:,:,matched));
      Q = cat(3, Q, Q_old(:,:,matched));
      R = cat(3, R, R_old(:,:,matched));
      objects = cat(2, objects, objects_old(matched));
      targetBox = cat(2, targetBox, targetBox_old(matched));
      
      x_old = x_old(:,~matched);
      P_old = P_old(:,:,~matched);
      Q_old = Q_old(:,:,~matched);
      R_old = R_old(:,:,~matched);
      objects_old = objects_old(~matched);
      targetBox_old = targetBox_old(~matched);
      
      NumberOfObjects = length(objects);
      numberOfNewObjects = NumberOfObjects - prevNumberOfObjects;
      
      if numberOfNewObjects > 0
         if prevNumberOfObjects > 0
            pi = prevNumberOfObjects/NumberOfObjects*pi;
         end
         pi = cat(1,pi,repmat(1/NumberOfObjects, numberOfNewObjects, 1));
      end
   end

%% Keep older objects
%    x_old = cat(2, x_old, x_remove);
%    P_old = cat(3, P_old, P_remove);
%    Q_old = cat(3, Q_old, Q_remove);
%    R_old = cat(3, R_old, R_remove);
%    objects_old = cat(2,objects_old, objects_remove);
%    targetBox_old = cat(2, targetBox_old, targetBox_remove);
%% Discard older objects
if i>startFrame+config.freq
   x_old = x_remove;
   P_old = P_remove;
   Q_old = Q_remove;
   R_old = R_remove;
   objects_old = objects_remove;
   targetBox_old = targetBox_remove;
end

   %% Calculate pi
   NumberOfObjects = length(objects);
   pi = ones(NumberOfObjects,1);%/NumberOfObjects;
   
   %% VISUALIZATION
  set(0,'CurrentFigure',handlePosFig)
  subplot(handlePos);
         
  hold on
   
  text(600,10,['Number of objects: ' num2str(NumberOfObjects)],'Margin',1,'Color','w','FontSize',8,'BackgroundColor',[0.1 0.1 0.1]);
   
   for k=1:NumberOfObjects
      %% Position
      %figure(handlePosFig);
      
      color = objColors(( mod(objects(k)-1,size(objColors,1)) + 1 ),:);
      
      %if ~mod(objects(k),10)
         set(0,'CurrentFigure',handlePosFig)
         subplot(handlePos);
         
         hold on
         
         plot(x(1,k),x(3,k),'+','MarkerSize',3,'Color',color);
         text(targetBox(k).xmin,targetBox(k).ymin,num2str(objects(k)),'Margin',1,'Color','w','VerticalAlignment','bottom','FontSize',8,'BackgroundColor',[0.1 0.1 0.1]);
         rectBox = get_rect(targetBox(k));

         rectangle('Position', rectBox, 'LineWidth',1,'LineStyle','-','EdgeColor',color);

         hold off
      %end
      
      %% Show weights
         % figure(handlePosFig)
         
         subplot(handleWeight);
         
         hold on
         
         if k<=size(ws,2)
            ind = ws(:,k)>config.weightThreshold;
            weights = ws(ind,k);
            zw = z(:,ind);
         
            plot(zw(1,:),zw(2,:),'o','MarkerEdgeColor',color);
         else
            plot(x(1,k),x(3,k),'+','MarkerSize',3,'Color',color);
            rectangle('Position', rectBox, 'LineWidth',1,'LineStyle','-','EdgeColor','g');
         end
         
         hold off
   end
   %% Draw old objects
   %    figure(handlePosFig)
   %    subplot(handlePos);
   %    hold on
   %
   %    for k=1:length(objects_old)
   %       color = objColors(objects_old(k),:);
   %       plot(x_old(1,k),x_old(3,k),'+','MarkerSize',12,'Color',color);
   %       text(x_old(1,k)-5,x_old(3,k)-5,num2str(objects_old(k)),'Color','w','VerticalAlignment','bottom');
   %
   %       rectBox = get_rect(targetBox_old(k));
   %       rectangle('Position',rectBox,'LineWidth',1,'LineStyle','--','EdgeColor',color);
   %
   %       x_old(1,k) = x_old(1,k) + x_old(2,k);
   %       x_old(3,k) = x_old(3,k) + x_old(4,k);
   %       w = targetBox_old(k).xmax-targetBox_old(k).xmin;
   %       h = targetBox_old(k).ymax-targetBox_old(k).ymin;
   %       targetBox_old(k) = get_bbox(x_old(1,k),x_old(3,k),w,h);
   %    end
   %
   %    hold off
   drawnow
   
   %% SAVE AVI
   if config.saveAVI
      drawnow;
      
      %FRAME = getframe(handlePosFig);
      %mov = addframe(mov,FRAME);
      %writeVideo(writer, FRAME);
      set(handlePosFig,'PaperPositionMode','auto');
      %cdata = hardcopy(handlePosFig, '-Dzbuffer', '-r0');
      cdata = print('-RGBImage');
      writeVideo(writer, cdata);
   end
   
   matFile = fullfile(strMATLAB, sprintf('emkf_MSER_%04d.mat', i));
   save(matFile, 'x', 'R', 'Q', 'P', 'targetBox', 'objects', ...
                 'x_old', 'R_old', 'Q_old', 'P_old', 'targetBox_old', 'objects_old', ...
                 'x_remove', 'R_remove', 'Q_remove', 'P_remove', 'targetBox_remove', 'objects_remove', ...
                 'z', 'ws', 'ws_remove', 'objColors', 'pi',...
                 'config', '-v7');
              
   if NumberOfObjects == 0
      break;
   end
   %pause(0.1)
end %END FRAMES


%% CLOSE AVI
if config.saveAVI
   %mov = close(mov);
   close(writer);
end

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
