function [ boundingBox ] = get_bbox( xc, yc, width, height )

boundingBox.xmin = floor(xc-width/2);
boundingBox.xmax = floor(xc+width/2);
boundingBox.ymin = floor(yc-height/2);
boundingBox.ymax = floor(yc+height/2);

end