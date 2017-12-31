% Main code for obtain information from processed frame(image).
% By Rui Luo 2017/12/26
clear
clc
close all

% main functions' switch, '1' means 'on' while '0' means 'off'
doCalculateRingAverageConcentration = 0;
doCalculateShearViscosity = 1;
doCalculateNormalViscosity = 0;
doCalculateRingAverageShearViscosity = 1;
doCalculateRingAverageNormalViscosity = 0;
doCalculateShearMobilityRatio = 1;
doCalculateNormalMobilityRatio = 0;
doFindCluster = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables to define.
phiInitial = 0.22;
indexFrameToCal = 30:15:450;

medianFilterThreshold = [3,3];
pixelUpperLimit=100;
cannyThreshold=[0.4,0.5];
cannySigma= sqrt(2);
denom = 3;
ringWidth = 5;  % unit: pixel
DataDirectory = 'C:\Users\lr546\Desktop\large particle\';
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
mobilityRatioMax = zeros(1,numTotalFrame);
for mCal = 1:numTotalFrame
    imageCalDirectory = [imageDirectory,num2str(indexFrameToCal(mCal)),'.png'];
    imageIntensity = imread(imageCalDirectory);
%     figure
%     imshow(imageIntensity)

% 2. calculate concentration at each pixel.
    [concentrationCartesian,concentrationPolar,imageConcentration] = getPixelConcentration(imageIntensity,imageReferenceAverageValue,inletRowPosition,inletColumnPosition,phiInitial);
    
% 3.1 calculate shear viscosity at each pixel.
if doCalculateShearViscosity == 1
    [shearViscosityCartesian,shearViscosityPolar] = getShearViscosity(concentrationCartesian,concentrationPolar);
end
% 3.2 calculate normal viscosity at each pixel.
if doCalculateNormalViscosity == 1
    [normalViscosityCartesian,normalViscosityPolar] = getNormalViscosity(concentrationCartesian,concentrationPolar);
end
% 4.1 calculate ring average shear viscosity.
if doCalculateRingAverageShearViscosity == 1
    shearViscosityRingAverage = getRingAverageValue(shearViscosityPolar,ringWidth);
end
% 4.2 calculate ring average normal viscosity.
if doCalculateRingAverageNormalViscosity == 1
    normalViscosityRingAverage = getRingAverageValue(normalViscosityPolar,ringWidth);
end
% 4.3 calculate ring average concentration.
if doCalculateRingAverageConcentration == 1
    concentrationRingAverage = getRingAverageValue(concentrationPolar,ringWidth);
end
% 5. calculate mobility ratio
if doCalculateShearMobilityRatio == 1
    shearMobilityRatioRingAverage = getMobilityRatio(shearViscosityRingAverage);
end
mobilityRatioMax(mCal) = max(shearMobilityRatioRingAverage(:,2));

if doCalculateNormalMobilityRatio == 1
    normalMobilityRatioRingAverage = getMobilityRatio(normalViscosityRingAverage); 
end
% 6. find cluster number and calculate its area.
if doFindCluster == 1
    [imageCluster,numberCluster,areaCluster] = findCluster(imageIntensity,imageConcentration,phiInitial);
end
% 7. find cluster position.

end

mobilityRatioDirectory = [DataDirectory,'phi',num2str(phiInitial*100),'\mobility ratio.xls'];
output1 = {'frame','maximum mobility ratio'};
output2 = [indexFrameToCal',mobilityRatioMax'];
xlswrite(mobilityRatioDirectory,output1,1,'A1');
xlswrite(mobilityRatioDirectory,output2,1,'A2');

