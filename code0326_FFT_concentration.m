% for testing FFT result of ring concentration.
% by Rui Luo 2018/3/26

clear
clc
close all

particleSize = 125;  % unit: um
gapThickness = 1.397; % unit: mm
FlowRate = 150;      % unit: ml/min
phiInitial = [0.31];
timeFrame = [600];
ringLocationMin = 410;
ringLocationMax = 415;
mainDirectory = 'C:\Users\rui\Desktop\';
dataDirectory = [mainDirectory,num2str(particleSize),'particle ',...
      num2str(gapThickness),'gap\phi',num2str(phiInitial*100)];
parameterDirectory = [dataDirectory,'\data.xls'];
imageDirectory = [dataDirectory,'\Gray Image\'];
imageFileDirectory = [imageDirectory,num2str(timeFrame),'.png'];
inletRowPosition = xlsread(parameterDirectory,1,'B2');
inletColumnPosition = xlsread(parameterDirectory,1,'C2');
ratio = xlsread(parameterDirectory,1,'A2');

imageOriginal = imread(imageFileDirectory);
figure;imshow(imageOriginal)
imageGray = rgb2gray(imageOriginal);
imageHisteq = histeq(imageGray);
figure;imshow(imageHisteq)
[imageHeight,imageWidth] = size(imageHisteq);
kNum=1;
for rowNum=1:imageHeight
    for columnNum=1:imageWidth
        rowLocation(kNum)=(rowNum-inletRowPosition);  % unit: pixel
        columnLocation(kNum)=(columnNum-inletColumnPosition);
        intensity(kNum)=double(imageHisteq(rowNum,columnNum));
        kNum=kNum+1;
    end
end
[thetaLocation,rhoLocation]=cart2pol(columnLocation,rowLocation);

inRingCondition = logical(rhoLocation > ringLocationMin &...
                                rhoLocation <= ringLocationMax);
indexInRing = find(inRingCondition);
rhoRing = rhoLocation(indexInRing);
thetaRing = thetaLocation(indexInRing);
intensityRing = intensity(indexInRing);
[thetaUse,indexUse] = unique(thetaRing);
intensityUse = intensityRing(indexUse);
figure;plot(thetaUse,intensityUse)
test = medfilt1(resample(intensityUse,2,5),3);
test = smooth(test);
figure; plot(test)
fs = 50;
N = 400;
n = 0:N-1;
t = n/fs;
f = n*fs/N;
intensityFT = fft(test,N);
mag = abs(intensityFT);
figure;plot(f(1:N/2),mag(1:N/2),'*-')