function [ rect,xmin,ymin,w,h ] = get_rect( boundingBox )

w = boundingBox.xmax-boundingBox.xmin;
h = boundingBox.ymax-boundingBox.ymin;
xmin = boundingBox.xmin;
ymin = boundingBox.ymin;
rect = [xmin,ymin,w,h];

end