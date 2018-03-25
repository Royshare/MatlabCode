% this function is used for labeling clusters and obtaining their numbers and area.
% By Rui Luo 2017/12/29

function [imageCluster,numberCluster,areaCluster] = findCluster(imageIntensity,imageConcentration,phiInitial)
   medianFilterThreshold = [3,3];
pixelUpperLimit=100;
cannyThreshold=[0.1,0.7];
cannySigma= 4;
% 
%    figure;imshow(imageIntensity)
   
   indexCluster = find(imageConcentration < phiInitial*1.05); %1.3 * phiinitial
   imageCluster = rgb2gray(imageIntensity);
   imageCluster(indexCluster) = 255;
   figure;imshow(imageCluster) 
   binaryThreshold = graythresh(imageCluster);
   imageClusterBinary = imbinarize(imageCluster,binaryThreshold);
   con = imageConcentration;
   con(indexCluster) = 0;
   mean2(con(con~=0))
   max(max(con))
%    figure;imshow(imageClusterBinary)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     imageedge = bwperim(imageClusterBinary);
%     figure
%     imshow(imageedge)
%    imageEdgeFilter = medfilt2(imageClusterBinary,medianFilterThreshold);
%    figure
%    imshow(imageEdgeFilter)
%    imageEdge = edge(imageEdgeFilter,'canny',cannyThreshold,cannySigma);
%    figure
%    imshow(imageEdge)
%    imageEdgeRemovedSmallObject = bwareaopen(imageEdge,pixelUpperLimit);
%    figure 
%    imshow(imageEdgeRemovedSmallObject)
%    edgeLength = find(imageEdgeRemovedSmallObject > 0);
%    size(edgeLength)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   imageClusterBinary = ~imageClusterBinary;
   bw1 = imfill(imageClusterBinary,'holes');
%    figure
%    imshow(bw1)
   se1 = strel('disk',3);
%    bw2 = imerode(bw1,se1);
%    figure
%    imshow(bw2)
   bw2 = imdilate(bw1,se1);

%    bw3 = imclose(bw1,se1);
%    figure
%    imshow(bw3)
%    imageEdge = edge(bw3,'canny',cannyThreshold,cannySigma);
%    figure
%    imshow(imageEdge)
%    imageEdgeRemovedSmallObject = bwareaopen(imageEdge,pixelUpperLimit);
%    figure 
%    imshow(imageEdgeRemovedSmallObject)
   se2 = strel('disk',10);
   bw3 = imopen(bw2,se2);
%    figure
%    imshow(bw3)
   [Location,clusterNumber] = bwlabel(bw3,4);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    figure;imshow(imageCluster)
%    [fx,fy] = gradient(imageCluster);
%    gradientTotal = sqrt(fx.^2+fy.^2);
%    figure;imshow(gradientTotal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    hold on
%    for labelNum = 1:clusterNumber
%        [row,column] = find(Location == labelNum);
%        rowCenter = mean(row);
%        columnCenter = mean(column);
%        plot(columnCenter,rowCenter,'marker','*','markeredgecolor','r','markersize',10);
%        text(columnCenter+20,rowCenter+20,num2str(labelNum));
%    end
%    hold off
   numberCluster = clusterNumber;
   areaCluster = 1;
end