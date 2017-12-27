function [rowCenter,columnCenter,inletImageRadius]=findInletCenter(video,FirstFrameIndex,cropVector,GreyImageWithEdgeDirectory,ColorImageDirectory)
medianFilterThreshold = [15,15];
pixelUpperLimit=100;
cannyThreshold=[0.4,0.5];
cannySigma=5;

imageColorOriginal = readframe(video,FirstFrameIndex-5);
imageGrayOrignial = rgb2gray(imageColorOriginal);
imageColorCrop = imcrop(imageColorOrignial,cropVector);
imageGrayCrop = imcrop(imageGrayOrignial,cropVector);
imageMedianFiltered = medfilt2(imageGrayCrop,medianFilterThreshold);
imageCannyEdgeDetected = edge(imageMedianFiltered,'canny',cannyThreshold,cannySigma);
imageRemovedSmallObject = bwareaopen(imageCannyEdgeDetected,pixelUpperLimit);
[rowEdge,columnEdge] = find(imageRemovedSmallObject);
rowCenter = mean(rowEdge);
columnCenter = mean(columnEdge);
temp = (rowEdge-rowCenter).^2+(columnEdge-columnCenter).^2;
inletImageRadius = mean(sqrt(temp));

imshow(imageGrayCrop);
hold on
plot(columnEdge,rowEdge,'b*','markersize',3)
%set(gcf,'PaperUnits','inches','PaperPosition',SaveVector)
fname = sprintf('finger_0.png');
fullfileName=fullfile(GreyImageWithEdgeDirectory,fname);
print('-dpng',fullfileName,'-r100'); 

end