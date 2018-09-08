function [ x,y,bbox,boxInit, ccInit ] = CreateMeasurements( imgProb )
%CREATEMEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here

EnableVisualization = 0;

binImg = imgProb>0;
width = size(binImg,2);
height = size(binImg,1);

se = strel('disk',1);
nIter = 0;
maxIterErode = 15;
maxIterDilate = 0;
MaxElementSize = 1;

%IM=imclose(binImg,se);
IM = binImg;
IMp = bwpack(IM);
%[L,NUM] = bwlabeln(IM);
%stats = regionprops(L,'basic');
cc = bwconncomp(IM,4);
ccInit = cc;
NUM = cc.NumObjects;
stats = regionprops(cc,'basic');

if EnableVisualization
   figure(5);clf;
   rect = get(gcf,'Position');
   %set(gcf,'Position',[rect(1) rect(2) width height]);
   set(gca,'Position',[0 0 1 1]);
   figure(6);clf;
   rect = get(gcf,'Position');
   %set(gcf,'Position',[rect(1) rect(2) width height]);
   set(gca,'Position',[0 0 1 1]);
   figure(7);clf;
   rect = get(gcf,'Position');
   %set(gcf,'Position',[rect(1) rect(2) width height]);
   set(gca,'Position',[0 0 1 1]);
   figure(8);clf;
   rect = get(gcf,'Position');
   %set(gcf,'Position',[rect(1) rect(2) width height]);
   set(gca,'Position',[0 0 1 1]);

   figure(5);clf;
   imshow(IM,'Border','tight');
   
   figure(6);clf;
   imshow(IM,'Border','tight');
   
   figure(7);clf;
   imshow(IM,'Border','tight');
   
   figure(8);clf;
   imshow(IM,'Border','tight');
end

xy=cat(1,stats.Centroid);
   
x=xy(:,1);
y=xy(:,2);

BoundingBoxes=cat(1,stats.BoundingBox);
bbox=BoundingBoxes;
bbox(:,5)=0.0;
% hold on;
% plot(xy(:,1),xy(:,2),'x')
% hold off

%% DILATION
IMDILp = IMp;
disp('Dilation...');
while (NUM>1 && nIter<maxIterDilate)
   ElementSize = ceil(MaxElementSize*rand);
   se = strel('disk',ElementSize);
   IMDILp = imdilate(IMDILp,se,'ispacked');
   
   %% Erode randomly
   if rand < .0
      %disp('	Erode...');
      IMDILp = imerode(IMDILp,se,'ispacked',size(IM,1));
   else
      %disp('	No erosion...');
   end
   
   IMDIL = bwunpack(IMDILp, size(IM,1));

   %[L,NUM] = bwlabeln(IMDIL);
   %stats = regionprops(L,'basic');
   cc = bwconncomp(IMDIL,8);
   NUM = cc.NumObjects;
   stats = regionprops(cc,'basic');
  
   xy = cat(1,stats.Centroid);

   x = cat(1,x,xy(:,1));
   y = cat(1,y,xy(:,2));
   
   BoundingBoxes = cat(1,stats.BoundingBox);
   BoundingBoxes(:,5) = -nIter;
   bbox = cat(1,bbox,BoundingBoxes);

   if EnableVisualization
      L = labelmatrix(cc);
      RGB = label2rgb(L, hsv(max(L(:))), 'k', 'shuffle');

      figure(6);imshow(IMDIL,'Border','tight');
      hold on
      plot(xy(:,1),xy(:,2),'r+')
      hold off

      figure(7);imshow(RGB,'Border','tight');
      pause
   end
   
   nIter = nIter+1;
end

%% EROSION
nIter = 1;
IMERp = IMp;

disp('Erosion...');
while (NUM>0 && nIter<maxIterErode)
   ElementSize = ceil(MaxElementSize*rand);
   se = strel('disk',ElementSize);
   IMERp = imerode(IMERp,se,'ispacked',size(IM,1));
   
   %% Dilate randomly
   if rand < .0
      %disp('	Dilate...');
      IMERp = imdilate(IMERp,se,'ispacked');
   else
      %disp('	No dilation...');
   end
   
   IMER = bwunpack(IMERp, size(IM,1));
   
   %[L,NUM] = bwlabeln(IMER);
   cc = bwconncomp(IMER,8);
   NUM = cc.NumObjects;
   
   if NUM==0
      break;
   end
   
   stats = regionprops(cc,'basic');
   
   %% Centroid
   xy=cat(1,stats.Centroid);
   
%    tmp=[x;xy(:,1)];
%    x=tmp;
%    tmp=[y;xy(:,2)];
%    y=tmp;
   
   x = cat(1,x,xy(:,1));
   y = cat(1,y,xy(:,2));
   
   %% Bounding box
   BoundingBoxes = cat(1,stats.BoundingBox);
   BoundingBoxes(:,5) = nIter;
   %tmp=[bbox;BoundingBoxes];
   %bbox=tmp;
   bbox = cat(1,bbox,BoundingBoxes);
   
   if EnableVisualization
      L = double(labelmatrix(cc));
      RGB = label2rgb(L, hsv(max(L(:))), 'k', 'shuffle');
   
      figure(6);imshow(IMER,'Border','tight');
   
      hold on
      plot(xy(:,1),xy(:,2),'r+')
      hold off

      figure(7);imshow(RGB,'Border','tight');
      pause
   end
   
   nIter = nIter+1;
end

boxInit = abs(bbox(:,5))>1;
bbox = [bbox(:,1)-bbox(:,5) bbox(:,2)-bbox(:,5) bbox(:,3)+2*bbox(:,5) bbox(:,4)+2*bbox(:,5)]; 

if EnableVisualization
   figure(8);
   imshow(imgProb,'Border','tight');
   hold on
   plot(x,y,'r+')
   hold off
   drawnow;
   pause(0.5)
end