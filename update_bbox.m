function [ boundingBox ] = update_bbox( boundingBox, xc_new, yc_new )
%UPDATE_BBOX Summary of this function goes here

[w,h] = get_dims(boundingBox);
boundingBox = get_bbox(xc_new,yc_new,w,h);

end

