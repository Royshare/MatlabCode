% this function is used for get concentration information from gray image.
% By Rui Luo 2017/12/27
function [concentrationCartesian,concentrationPolar,imageConcentration] = getPixelConcentration(imageIntensity,imageReferenceAverageValue,inletRowPosition,inletColumnPosition,phiInitial)
pixelUpperLimitDelete = 40; % image inside this circle will be removed.
pixelLowerLimitDelete = 440;
[imageHeight,imageWidth,~] = size(imageIntensity); 

imageGrayIntensity = rgb2gray(imageIntensity);
imageBinarySuspension = imbinarize(imageGrayIntensity);

imageGrayIntensity(imageBinarySuspension) = 255;

% figure
% imshow(imageBinarySuspension)
% figure
% imshow(imageGrayIntensity)
%     
% obtain minimum intensity randomly
rowRandom = fix([0.02*imageHeight+0.15*imageHeight*rand(1),...
             0.83*imageHeight+0.15*imageHeight*rand(1),...
             0.02*imageHeight+0.15*imageHeight*rand(1),...
             0.83*imageHeight+0.15*imageHeight*rand(1)]);
columnRandom = fix([0.02*imageWidth+0.15*imageWidth*rand(1),...
                0.02*imageWidth+0.15*imageWidth*rand(1),...
                0.83*imageWidth+0.15*imageWidth*rand(1),...
                0.83*imageWidth+0.15*imageWidth*rand(1)]);
intensityRandom = zeros(1,4);
for iRandom = 1:4
    intensityRandom(iRandom) = imageGrayIntensity(rowRandom(iRandom),...
                               columnRandom(iRandom));
end
intensityMinimum = mean(intensityRandom);

imageIntensityModified = double(imageGrayIntensity)./double(imageReferenceAverageValue);
imageIntensityModified = uint8(imageIntensityModified*180);
% figure
% imshow(imageIntensityModified)

imageLog=double(imageIntensityModified);
imageLog=log(imageLog/intensityMinimum);
imageLog(imageBinarySuspension) = 0;
% delete data near inlet.
kNum=1;
for rowNum=1:imageHeight
    for columnNum=1:imageWidth
        distanceToCenter = sqrt((rowNum-inletRowPosition)^2+...
                           (columnNum-inletColumnPosition)^2);
        if distanceToCenter <= pixelUpperLimitDelete || distanceToCenter > pixelLowerLimitDelete
            indexDelete(1,kNum)=rowNum;
            indexDelete(2,kNum)=columnNum;
            kNum=kNum+1;
        end
    end
end
indexRowDelete=indexDelete(1,:);
indexColumnDelete=indexDelete(2,:);
indexDeleteArea=sub2ind(size(imageLog),indexRowDelete,indexColumnDelete);

imageLog(indexDeleteArea)=0;
% temp = imageGrayIntensity;
% temp(indexDeleteArea) = 255;
% figure
% imshow(temp)
% calculate area after deleting region near inlet.
indexCaluateSuspension=find(imageLog ~= 0);
imageCaltemp=zeros(size(imageLog));
imageCaltemp(indexCaluateSuspension)=1;
areaTemp=sum(imageCaltemp(:));
% calculate coefficient
k = phiInitial*areaTemp/sum(sum(imageLog(:)));
imageConcentration=k*imageLog;


kNum=1;
for rowNum=1:imageHeight
    for columnNum=1:imageWidth
        rowLocation(kNum)=(rowNum-inletRowPosition);  % unit: pixel
        columnLocation(kNum)=(columnNum-inletColumnPosition);
        concentration(kNum)=imageConcentration(rowNum,columnNum);
        kNum=kNum+1;
    end
end
[thetaLocation,rhoLocation]=cart2pol(columnLocation,rowLocation);
%%%
    indexNoSuspension = concentration == 0;
    concentrationSuspension = concentration(~indexNoSuspension);
    rhoSuspensionLocation = rhoLocation(~indexNoSuspension);
    thetaSuspensionLocation = thetaLocation(~indexNoSuspension);
    rowSuspensionLocation = rowLocation(~indexNoSuspension);
    columnSuspensionLocation = columnLocation(~indexNoSuspension);
    concentrationCartesian = [columnSuspensionLocation',rowSuspensionLocation',concentrationSuspension'];
    concentrationPolar = [rhoSuspensionLocation',thetaSuspensionLocation',concentrationSuspension'];
%%%

% concentrationCartesian = [columnLocation',rowLocation',concentration'];
% concentrationPolar = [rhoLocation',thetaLocation',concentration'];

% [conPlotx,conPloty]=meshgrid(1:imageWidth,1:imageHeight);
% figure
% pcolor(conPlotx,conPloty,imageConcentration)
% axis([400,750,300,650])
end