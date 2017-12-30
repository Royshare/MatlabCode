% Main code for obtain information from processed frame(image).
% By Rui Luo 2017/12/26
clear
clc
close all

% main functions' switch, '1' means 'on' while '0' means 'off'
doCalculateRingAverageConcentration = 1;
doCalculateShearViscosity = 1;
doCalculateNormalViscosity = 0;
doCalculateRingAverageShearViscosity = 1;
doCalculateRingAverageNormalViscosity = 0;
doFindClusterNumber = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables to define.
indexFrameToCal = [120];
medianFilterThreshold = [3,3];
pixelUpperLimit=100;
cannyThreshold=[0.4,0.5];
cannySigma= sqrt(2);
denom = 3;
phiInitial = 0.31;
ringWidth = 5;  % unit: pixel
DataDirectory = 'C:\Users\lr546\Desktop\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlsDirectory = [DataDirectory,'phi',num2str(phiInitial*100),'\data.xls'];
imageDirectory =[DataDirectory,'phi',num2str(phiInitial*100),'\Gray Image\'];
Ratio = xlsread(xlsDirectory,1,'A2');
inletRowPosition = xlsread(xlsDirectory,1,'B2');
inletColumnPosition = xlsread(xlsDirectory,1,'C2');
[~,numTotalFrame]=size(indexFrameToCal);


% 1. obtain reference image.
imageCalDirectory = [imageDirectory,'1.png'];
imageIntensity = imread(imageCalDirectory);
imageReferenceAverageValue = getReferenceImage(imageIntensity,DataDirectory,...
                                    inletRowPosition,inletColumnPosition);
% temp = uint8(imageReferenceAverageValue);
% figure
% imshow(temp)
for mCal = 1:numTotalFrame
    imageCalDirectory = [imageDirectory,num2str(indexFrameToCal(mCal)),'.png'];
    imageIntensity = imread(imageCalDirectory);
%     figure
%     imshow(imageIntensity)

% 2. calculate concentration at each pixel.
    [concentrationCartesian,concentrationPolar,imageConcentration] = getPixelConcentration(imageIntensity,imageReferenceAverageValue,inletRowPosition,inletColumnPosition,phiInitial);
    
% 3.1 calculate shear viscosity at each pixel.    
    [shearViscosityCartesian,shearViscosityPolar] = getShearViscosity(concentrationCartesian,concentrationPolar);
% 3.2 calculate normal viscosity at each pixel.
    [normalViscosityCartesian,normalViscosityPolar] = getNormalViscosity(concentrationCartesian,concentrationPolar);
% 4.1 calculate ring average shear viscosity.
    shearViscosityRingAverage = getRingAverageValue(shearViscosityPolar,ringWidth);
% 4.2 calculate ring average normal viscosity.
    normalViscosityRingAverage = getRingAverageValue(normalViscosityPolar,ringWidth);
% 4.3 calculate ring average concentration.
    concentrationRingAverage = getRingAverageValue(concentrationPolar,ringWidth);
% 5. calculate mobility ratio
    shearMobilityRatioRing = getMobilityRatio(shearViscosityRingAverage);
    normalMobilityRatioRing = getMobilityRatio(normalViscosityRingAverage);
% 6. find cluster number and calculate its area.
    [imageCluster,numberCluster,areaCluster] = findCluster(imageIntensity,imageConcentration,phiInitial)
% 7. find cluster position.

end


