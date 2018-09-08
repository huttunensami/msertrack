function [ w,h ] = get_dims( boundingBox )
%GET_DIMS Summary of this function goes here

w = boundingBox.xmax-boundingBox.xmin;
h = boundingBox.ymax-boundingBox.ymin;

end

