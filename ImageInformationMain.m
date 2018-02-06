% Main code for obtain information from processed frame(image).
% By Rui Luo 2017/12/26
clear
clc
close all

% main functions' switch, '1' means 'on' while '0' means 'off'
doCalculatePixelConcentration = 0;
doCalculateRingAverageConcentration = 0;
doCalculateShearViscosity = 0;
doCalculateNormalViscosity = 0;
doCalculateRingAverageShearViscosity = 0;
doCalculateRingAverageNormalViscosity = 0;
doCalculateShearMobilityRatio = 0;
doCalculateNormalMobilityRatio = 0;
doFindCluster = 0;
doCalculateOuterPerimeter = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables to define.
% phiInitialArray = [0.14,0.17,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35];
phiInitialArray = linspace(0.22,0.36,15);
for informationGetIndex = 1:length(phiInitialArray)
phiInitial = phiInitialArray(informationGetIndex);
indexFrameToCal = 30:1:600;

ringWidth = 5;  % unit: pixel
DataDirectory = 'C:\Users\lr546\Desktop\125particle 1.397gap\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameterDirectory = [DataDirectory,'phi',num2str(phiInitial*100),'\data.xls'];
imageDirectory =[DataDirectory,'phi',num2str(phiInitial*100),'\Gray Image\'];
Ratio = xlsread(parameterDirectory,1,'A2'); % unit: pixel/cm
inletRowPosition = xlsread(parameterDirectory,1,'B2');
inletColumnPosition = xlsread(parameterDirectory,1,'C2');
[~,numTotalFrame]=size(indexFrameToCal);


% 1. obtain reference image.
if doCalculatePixelConcentration == 1
      imageCalDirectory = [imageDirectory,'1.png'];
      imageIntensity = imread(imageCalDirectory);
      imageReferenceAverageValue = getReferenceImage(imageIntensity,DataDirectory,...
                                    inletRowPosition,inletColumnPosition);
      % temp = uint8(imageReferenceAverageValue);
      % figure
      % imshow(temp)
end
shearMobilityRatioMax = zeros(1,numTotalFrame);
shearMobilityRatioMaxLocation = zeros(1,numTotalFrame);
shearViscosityContrastMax = zeros(1,numTotalFrame);
numberCluster = zeros(1,numTotalFrame);
for mCal = 1:numTotalFrame
    imageCalDirectory = [imageDirectory,num2str(indexFrameToCal(mCal)),'.png'];
    imageIntensity = imread(imageCalDirectory);
%     figure
%     imshow(imageIntensity)

% 2. calculate concentration at each pixel.
if doCalculatePixelConcentration == 1
    [concentrationCartesian,concentrationPolar,imageConcentration] = getPixelConcentration(imageIntensity,imageReferenceAverageValue,inletRowPosition,inletColumnPosition,phiInitial);
end    
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
    [shearMobilityRatioMax(mCal),shearMobilityRatioMaxIndex] = max(shearMobilityRatioRingAverage(:,2));
    shearMobilityRatioMaxLocation(mCal) = shearMobilityRatioRingAverage(shearMobilityRatioMaxIndex,1);
end

if doCalculateNormalMobilityRatio == 1
    normalMobilityRatioRingAverage = getMobilityRatio(normalViscosityRingAverage); 
end
% 5.2 calculate viscosity contrast (Homsy)
% if doCalculateShearViscosityContrast ==1
%    viscosityContrastRingAverage = getViscosityContrast(shearViscosityRingAverage);
%    shearViscosityContrastMax(mCal) = max(viscosityContrastRingAverage(:,2));
% end
% 6. find cluster number and calculate its area.
if doFindCluster == 1
    [imageCluster,numberCluster,areaCluster] = findCluster(imageIntensity,imageConcentration,phiInitial);
%      set(gcf,'PaperPosition',ImageSizeVector/100);
    ClusterImageFullDirectory = [ClusterImageDirectory,'\',num2str(indexFrameToCal(mCal)),'.png'];
    print('-dpng',ClusterImageFullDirectory,'-r100');
end
% 7. find outer perimeter.
if doCalculateOuterPerimeter == 1
    interfaceOuter = getOuterInterfacePosition(imageIntensity,inletRowPosition,inletColumnPosition);
%     interfaceDirectory = [DataDirectory,'phi',num2str(phiInitial*100),'\interface all.xls'];
    interfaceDirectory = [DataDirectory,'phi',num2str(phiInitial*100),'\interface data'];
    mkdir(interfaceDirectory);
    interfaceFileName = [interfaceDirectory,'\',num2str(indexFrameToCal(mCal)),'.csv'];
    output = array2table(interfaceOuter,'VariableNames',{'theta','rho'});
    writetable(output,interfaceFileName)
%     sheetName=['Sheet',num2str(mCal)];
%     output1 = {'frame',num2str(indexFrameToCal(mCal))};
%     output2 = {'theta','rho'};
%     output3 = interfaceOuter;
%     xlswrite(interfaceDirectory,output1,sheetName,'A1');
%     xlswrite(interfaceDirectory,output2,sheetName,'A2');
%     xlswrite(interfaceDirectory,output3,sheetName,'A3');
    
end
      
end
if doCalculateShearMobilityRatio == 1
      tempSmoothed = smooth(shearMobilityRatioMax);
      mobilityRatioDirectory = [DataDirectory,'phi',num2str(phiInitial*100),'\mobility ratio.xls'];
      output1 = {'frame','maximum mobility ratio','location/pixel','maximum mobility ratio cleaned'};
      output2 = [indexFrameToCal',shearMobilityRatioMax',shearMobilityRatioMaxLocation',tempSmoothed];
      xlswrite(mobilityRatioDirectory,output1,1,'A1');
      xlswrite(mobilityRatioDirectory,output2,1,'A2');
end

% Save cluster number result data.
if doFindCluster == 1
      clusterNumberDirectory = [DataDirectory,'phi',num2str(phiInitial*100),'\cluster number.xls'];
      output1 = {'frame','cluster number'};
      output2 = [indexFrameToCal',numberCluster'];
      xlswrite(clusterNumberDirectory,output1,1,'A1');
      xlswrite(clusterNumberDirectory,output2,1,'A2');
end
end 
