% this function is used for crop frames needed in video and transfer them
% into gray images.
% By Rui LUO 2017/12/27
function [imageGrayCrop,imageColorCrop] = cropImage(VideoToProcess,...
                                            indexToProcessFrame,CropVector)
imageColorOriginal = readframe(VideoToProcess,indexToProcessFrame);
imageGrayOrignial = rgb2gray(imageColorOriginal);
imageColorCrop = imcrop(imageColorOrignial,CropVector);
imageGrayCrop = imcrop(imageGrayOrignial,CropVector);
end