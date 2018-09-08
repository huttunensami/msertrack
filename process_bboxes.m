function [ newBoxes, numberOfMatches, covar ] = process_bboxes( bboxs )
%process_bboxes Post processes bounding boxes
%   Detailed explanation goes here


newBoxes = [];
numberOfMatches = [];
covar = [];

if size(bboxs,1)==0
   bContinue = 0;
else
   bContinue = 1;
end

remainingBoxes = bboxs;
intersectLimit = 0.95;

while(bContinue)
   %% Calculate score
   rectBox = remainingBoxes(1,:);
   rectArea = repmat(rectBox(3)*rectBox(4),size(remainingBoxes,1),1);
   areaUnion = rectunion(rectBox,remainingBoxes)';
   areaIntersect = rectint(rectBox,remainingBoxes)';
   rectScore = areaIntersect./areaUnion;
   
   %% Find matches
   bboxAreas = remainingBoxes(:,3).*remainingBoxes(:,4);
   intersects = (areaIntersect>=intersectLimit*bboxAreas) | (areaIntersect>=intersectLimit*rectArea);
   %intersects = ( (areaIntersect>=intersectLimit*bboxAreas) & (areaIntersect<bboxAreas) )...
   %             | ( (areaIntersect>=intersectLimit*rectArea) & (areaIntersect<rectArea) );
   matches = rectScore>0.25 | intersects;
   matchBoxes = remainingBoxes(matches,:);
   %rectScore = rectScore(matches);
   nMatches = size(matchBoxes,1);
   
   %% Calculate statistics
   x = matchBoxes(:,1) + matchBoxes(:,3)/2;
   y = matchBoxes(:,2) + matchBoxes(:,4)/2;
   
   newCovar = [var(x)   0;
               0        var(y)];
   
   %% Calculate and add a new box to the list
   %newBox = rectScore'*matchBoxes/sum(rectScore);
   xc = sum(x)/nMatches;
   yc = sum(y)/nMatches;
   dims = sum(matchBoxes(:,3:4),1)/nMatches;
   newBox = [xc-dims(1)/2 yc-dims(2)/2 dims];
   remainingBoxes = remainingBoxes(~matches,:);
   
   newBoxes = cat(1, newBoxes, newBox);
   numberOfMatches = cat(1, numberOfMatches, nMatches);
   covar = cat(3, covar, newCovar);
   
   if size(remainingBoxes,1) < 1
      bContinue = 0;
   end
end

for i = 1:size(bboxs,1)
   rectangle('Position', bboxs(i,:), 'LineWidth',1,'LineStyle','-','EdgeColor','g');
end

for i = 1:size(newBoxes,1)
   rectangle('Position', newBoxes(i,:), 'LineWidth',1,'LineStyle','-','EdgeColor','r');
end
end