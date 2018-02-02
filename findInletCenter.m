% this function is used for finding inlet position in image.
% By Rui Luo 2017/12/25
function [rowCenter,columnCenter,inletImageRadius]=findInletCenter(video,...
                                  FirstFrameIndex,CropVector,DataDirectory)

medianFilterThreshold = [5,5];
pixelUpperLimit=100;
cannyThreshold=[0.7,0.8];
cannySigma= 10;
ImageSizeVector = [0,0,CropVector(3),CropVector(4)];
GrayImageWithEdgeDirectory = fullfile(DataDirectory,'Gray Image with Edge');
ColorImageDirectory = fullfile(DataDirectory,'Color Image');

imageColorOriginal = read(video,FirstFrameIndex);
imageGrayOrignial = rgb2gray(imageColorOriginal);
imageColorCrop = imcrop(imageColorOriginal,CropVector);
imageGrayCrop = imcrop(imageGrayOrignial,CropVector);
imageMedianFiltered = medfilt2(imageGrayCrop,medianFilterThreshold);
imageCannyEdgeDetected = edge(imageMedianFiltered,'canny',cannyThreshold,cannySigma);
imageRemovedSmallObject = bwareaopen(imageCannyEdgeDetected,pixelUpperLimit);
[rowEdge,columnEdge] = find(imageRemovedSmallObject);
rowCenter = mean(rowEdge);
columnCenter = mean(columnEdge);
temp = (rowEdge-rowCenter).^2+(columnEdge-columnCenter).^2;
inletImageRadius = mean(sqrt(temp));

imshow(imageGrayCrop,'border','tight','initialmagnification','fit');
hold on
plot(columnEdge,rowEdge,'b*','markersize',1)
set(gcf,'PaperPosition',ImageSizeVector/100)

fname = sprintf('0.png');
GrayImageFullDirectory = fullfile(GrayImageWithEdgeDirectory,fname);
print('-dpng',GrayImageFullDirectory,'-r100'); 

hold off
figure
imshow(imageColorCrop,'border','tight','initialmagnification','fit');
set(gcf,'PaperPosition',ImageSizeVector/100)
fname = sprintf('0.png');
ColorImageFullDirectory = fullfile(ColorImageDirectory,fname);
print('-dpng',ColorImageFullDirectory,'-r100');
end