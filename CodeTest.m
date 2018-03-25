%%
clear
clc
close all
x = linspace(0,200,100);
y = log(x)-1000*x-50*x.^-0.5;
z = -x+100*log(x);
figure;plot(x,y)
figure;plot(x,z)
%% 
clear
clc
close all
x = -20:0.2:20;
k = 2;
y = sin(k*x)/pi./x;
figure; plot(x,y)
k = 20;
y = sin(k*x)/pi./x;
figure; plot(x,y)
%% 
clear
clc
close all

nx=1000;xmin=-15;xmax=15;
ny=1000;ymin=-15;ymax=15;
[x,y]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
a = 3; h = 1;
Q = 10;

xx = x.^4-6*x.^2*y.^2+y.^4+2*(x.^2-y.^2)*(h^2-a^2)+a^4+2*a^2*h^2+h^4;
yy = 4*x.*y*(x.^2-y.^2+h^2-a^2);
potentialfun = Q/4/pi*log(xx.^2+yy.^2);
% potentialfun = Q/4/pi*log((x.^2-a^2-h^2-y.^2+2*h*y)^2+4*x.^2*(y-h)^2);
levmax=max(max(potentialfun));
levmin=min(min(potentialfun));
lev=linspace(levmin,levmax,100);
figure;
contour(x,y,potentialfun,lev)
title('Equipotential line of two source');

streamfun = Q/2/pi*atan2(yy,xx);
% streamfun = Q/2/pi*atan2(2*x.*(y-h),(x.^2-a^2-h^2-y^2+2*h*y));
lev1max=max(max(streamfun));
lev1min=min(min(streamfun));
lev1=linspace(lev1min,lev1max,100);
figure;
contour(x,y,streamfun,lev1)
title('Steamline of two source');


%% find interface inside
clear
clc
close all
particleSize = 125;  % unit: um
gapThickness = 1.397; % unit: mm
phiInitial = 0.2;
FlowRate = 150;      % unit: ml/min
pixelUpperLimitDelete = 40; % image inside this circle will be removed.
pixelLowerLimitDelete = 440;

MainDirectory = 'C:\Users\lr546\Desktop\';
DataDirectory = [MainDirectory,num2str(particleSize),'particle ',...
      num2str(gapThickness),'gap\phi',num2str(phiInitial*100)];
parameterDirectory = [DataDirectory,'\data.xls'];
imageDirectory =[DataDirectory,'\Gray Image\'];
Ratio = xlsread(parameterDirectory,1,'A2'); % unit: pixel/cm
inletRowPosition = xlsread(parameterDirectory,1,'B2');
inletColumnPosition = xlsread(parameterDirectory,1,'C2');

timeFrame = [450];
indexTotal = length(timeFrame);
for indexAnalysis = 1:length(timeFrame)
imageCalDirectory = [imageDirectory,num2str(timeFrame(indexAnalysis)),'.png'];
imageIntensity = imread(imageCalDirectory);

[imageHeight,imageWidth,~] = size(imageIntensity); 
imageGrayIntensity = rgb2gray(imageIntensity);
imageAnalysis = ezsmoothn(double(imageGrayIntensity));


image1 = imageAnalysis;
[fx,fy] = gradient(image1);

figure; imshow(image1)

imageOut = edge(image1,'canny');
figure; imshow(imageOut)
figure;imshow(fx)
figure;imshow(fy)
figure; imshow(fx.^2+fy.^2)
end
%% test perimeter of interface
clear
clc
close all


MainDirectory = 'C:\Users\lr546\Desktop\';
phiInitialArray = [0.17,0.2,0.25,0.28,0.31];
for informationGetIndex = 1:length(phiInitialArray)
      particleSize = 125;  % unit: um
      gapThickness = 1.397; % unit: mm
      phiInitial = phiInitialArray(informationGetIndex);
      FlowRate = 150;      % unit: ml/min
      timeFrame = linspace(30,600,571);
      % timeFrame = linspace(6,20,15).^2;
      DataDirectory = [MainDirectory,num2str(particleSize),'particle ',...
      num2str(gapThickness),'gap\phi',num2str(phiInitial*100)];
      interfaceDirectory = [DataDirectory,'\interface data'];
      parameterDirectory = [DataDirectory,'\data.xls'];
      ratio = xlsread(parameterDirectory,1,'A2')-3.7; % unit: pixel/cm

      indexTotal = length(timeFrame);
      s = zeros(1,indexTotal);
for indexAnalysis = 1:indexTotal
    fileDirectory = [interfaceDirectory,'\',num2str(timeFrame(indexAnalysis)),'.csv'];
    interface = csvread(fileDirectory,1,0);
    theta = interface(:,1);
    rho = interface(:,2);
    theta1 = theta;
    theta1(1,:) = [];
    rho1 = rho;
    rho1(1,:) = [];
    pointNumber = length(theta);
    theta(pointNumber) = [];
    rho(pointNumber) = []; 
    radiusCircle1 = ratio*sqrt(FlowRate/60 * timeFrame(indexAnalysis)/30/pi/(gapThickness/10));
    s(indexAnalysis) = sum(abs((theta-theta1).*(rho+rho1)/2));
    s1 = s;
    s2 = s;
    s1(:,1) = [];
    s2(:,indexTotal) = [];
    sdiff = s2./s1;
    rhoMean(indexAnalysis) = mean(rho1);
end
time = timeFrame;
time(:,indexTotal) = [];
% semilogx(time,sdiff,'-*')
% plot(sqrt(timeFrame),rhoMean);
loglog(timeFrame,s)
hold on
end
hold off
%% test smooth setting
clear
clc
close all

particleSize = 125;  % unit: um
gapThickness = 1.397; % unit: mm
phiInitial = 0.17;
FlowRate = 150;      % unit: ml/min
MainDirectory = 'C:\Users\lr546\Desktop\';


DataDirectory = [MainDirectory,num2str(particleSize),'particle ',...
      num2str(gapThickness),'gap\phi',num2str(phiInitial*100)];
interfaceDirectory = [DataDirectory,'\interface data'];
parameterDirectory = [DataDirectory,'\data.xls'];
% ratio = xlsread(parameterDirectory,1,'A2')-3.6; % unit: pixel/cm
% timeFrame = linspace(30,600,571);
timeFrame = 400;

fileDirectory = [interfaceDirectory,'\',num2str(timeFrame),'.csv'];
      interface1 = csvread(fileDirectory,1,0);
theta1 = interface1(:,1);
rho1 = interface1(:,2);
rho1S = smooth(rho1,'sgolay');
rho1S = smooth(rho1S,'sgolay');
rho1Sm = ezsmoothn(rho1);
      figure;plot(theta1,rho1)
      figure;plot(theta1,rho1S)
      figure;plot(theta1,rho1Sm)
%% test output csv
clear
clc
close all
output2 = {'theta','rho','rho_nomalized'};
output3 = [1,23,23;2,34,34];
output = array2table(output3,'VariableNames',{'theta','rho','rho_nomalized'});
fileName = 'test.csv';
writetable(output,fileName);
% fid = fopen(fileName,'w');
% fprintf(fid,'%s,%s,%s\n',output2{:});
% dlmwrite(fileName,output,',',1,0);
% fclose(fid);
%% using hausdorff distance to determine similarity between curve
clear
clc
close all

particleSize = 125;  % unit: um
gapThickness = 1.397; % unit: mm
phiInitial = 0.17;
MainDirectory = 'C:\Users\lr546\Desktop\';
DataDirectory = [MainDirectory,num2str(particleSize),'particle ',...
      num2str(gapThickness),'gap\phi',num2str(phiInitial*100)];
InterfaceDirectory = [DataDirectory,'\interface1.xls'];

comparisonMax = sheetNumber-1;
hd = zeros(1,comparisonMax);
for indexAnalysis = 1:indexTotal-1
      fileDirectory = [interfaceDirectory,'\',num2str(timeFrame),'.csv'];
      interfaceSheetFinal = xlsread(InterfaceDirectory,sheetNumber);
      interfaceSheetFinal(1:2,:) = [];
      interfaceFinalSmooth = ezsmoothn(interfaceSheetFinal(:,3));
      interfaceFinal = interp1(interfaceSheetFinal(:,1),interfaceFinalSmooth,thetaInterpolation);
%       interfaceFinal = interp1(interfaceSheetFinal(:,1),interfaceSheetFinal(:,3),thetaInterpolation);
      sample1 = [thetaInterpolation',interfaceFinal'];
%       interfaceFinal = [interfaceSheetFinal(:,1)',interfaceSheetFinal(:,3)'];
      interfaceSheetCompare = xlsread(InterfaceDirectory,comparisonIndex);
      interfaceSheetCompare(1:2,:)=[];
      interfaceCompareSmooth = ezsmoothn(interfaceSheetCompare(:,3));
      interfaceCompare = interp1(interfaceSheetCompare(:,1),interfaceCompareSmooth,thetaInterpolation);
%       interfaceCompare = interp1(interfaceSheetCompare(:,1),interfaceSheetCompare(:,3),thetaInterpolation);
      sample2 = [thetaInterpolation',interfaceCompare'];
      %       interfaceCompared = [interfaceSheetCompare(:,1)',interfaceSheetFinal(:,3)'];
      [hd(comparisonIndex),csq] = HausdorffDist(sample1,sample2);
      figure;plot(thetaInterpolation,interfaceFinal)
      hold on
      plot(thetaInterpolation,interfaceCompare)
end
%% USING Frechet Distance to determine similarity between curve
clear
clc
close all
thetaInterpolation = linspace(-3.13,3.13,1800);
particleSize = 125;  % unit: um
gapThickness = 1.27; % unit: mm
phiInitial = 0.35;
MainDirectory = 'C:\Users\lr546\Desktop\';
DataDirectory = [MainDirectory,num2str(particleSize),'particle ',...
      num2str(gapThickness),'gap\phi',num2str(phiInitial*100)];
InterfaceDirectory = [DataDirectory,'\interface1.xls'];
[~,sheetName] = xlsfinfo(InterfaceDirectory);
[~,sheetNumber] = size(sheetName);
comparisonMax = sheetNumber-1;
hd = zeros(1,comparisonMax);
for comparisonIndex = comparisonMax:-1:70
      interfaceSheetFinal = xlsread(InterfaceDirectory,sheetNumber);
      interfaceSheetFinal(1:2,:) = [];
      interfaceFinalSmooth = smooth(interfaceSheetFinal(:,3));
      interfaceFinal = interp1(interfaceSheetFinal(:,1),interfaceFinalSmooth,thetaInterpolation);
%       interfaceFinal = interp1(interfaceSheetFinal(:,1),interfaceSheetFinal(:,3),thetaInterpolation);
      sample1 = [thetaInterpolation',interfaceFinal'];
%       interfaceFinal = [interfaceSheetFinal(:,1)',interfaceSheetFinal(:,3)'];
      interfaceSheetCompare = xlsread(InterfaceDirectory,comparisonIndex);
      interfaceSheetCompare(1:2,:)=[];
      interfaceCompareSmooth = smooth(interfaceSheetCompare(:,3));
      interfaceCompare = interp1(interfaceSheetCompare(:,1),interfaceCompareSmooth,thetaInterpolation);
%       interfaceCompare = interp1(interfaceSheetCompare(:,1),interfaceSheetCompare(:,3),thetaInterpolation);
      sample2 = [thetaInterpolation',interfaceCompare'];
      %       interfaceCompared = [interfaceSheetCompare(:,1)',interfaceSheetFinal(:,3)'];
      [hd(comparisonIndex),csq] = DiscreteFrechetDist(sample1,sample2);
      figure;plot(thetaInterpolation,interfaceFinal)
      hold on
      plot(thetaInterpolation,interfaceCompare)
end

%% test Frechet Distance
clear
clc
close all

t = 0:pi/8:2*pi;
y = linspace(1,3,6);
P = [(2:7)' y']+0.3.*randn(6,2);
Q = [t' sin(t')]+2+0.3.*randn(length(t),2);
P1 = [1,1;2,2;3,3;4,4;5,5];
Q1 = [1,2;2,3;3,4;4,5;5,6];
[hd, cSq] = DiscreteFrechetDist(P,Q);
% plot result
figure
plot(Q1(:,1),Q1(:,2),'o-r','linewidth',3,'markerfacecolor','r')
hold on
plot(P1(:,1),P1(:,2),'o-b','linewidth',3,'markerfacecolor','b')
title(['Discrete Frechet Distance of curves P and Q: ' num2str(hd)])
legend('Q','P','location','best')
line([2 hd+2],[0.5 0.5],'color','m','linewidth',2)
text(2,0.4,'dFD length')
for i=1:length(cSq)
  line([P1(cSq(i,1),1) Q1(cSq(i,2),1)],...
       [P1(cSq(i,1),2) Q1(cSq(i,2),2)],...
       'color',[0 0 0]+(i/length(cSq)/1.35));
end
axis equal
% display the coupling sequence along with each distance between points
disp([cSq sqrt(sum((P1(cSq(:,1),:) - Q1(cSq(:,2),:)).^2,2))])
%% test deformation
clear
clc
close all

phiInitial = 0.31;
DataDirectory = 'C:\Users\lr546\Desktop\125particle 1.397gap\';
interfaceDirectory = [DataDirectory,'phi',num2str(phiInitial*100),'\interface.xls'];
sheetIndex = 15;
interfaceDataFrame = xlsread(interfaceDirectory,sheetIndex);
interfaceDataFrame([1,2],:) = [];
interfaceTheta = interfaceDataFrame(:,1);
rho = interfaceDataFrame(:,2);
interface = interfaceDataFrame(:,3);
 figure;plot(interfaceTheta,interface)
 figure;plot(interfaceTheta,rho)

% interface1 = smooth(interface);
interface2 = smooth(interface,'lowess');
figure;plot(interfaceTheta,interface2)
interface2 = smooth(interface2,'lowess');
figure;plot(interfaceTheta,interface2)

% interface3 = smooth(interface,'loess');
% interface4 = smooth(interface,'sgolay');
% interface5 = smooth(interface,'rlowess');
% interface6 = smooth(interface,'rloess');

% figure;plot(interfaceTheta,interface1)

% figure;plot(interfaceTheta,interface3)
% figure;plot(interfaceTheta,interface4)
% figure;plot(interfaceTheta,interface5)
% figure;plot(interfaceTheta,interface6)

interfacem = medfilt1(interface2,10);
figure;plot(interfaceTheta,interfacem)
interfacem = medfilt1(interfacem,5);
figure;plot(interfaceTheta,interfacem)
test = 3-interfacem;
figure;plot(interfaceTheta,test)
[PKS,LOCS]= findpeaks(test);
interfaceTheta(LOCS)
%% test smooth function use
clear
clc
close all

%% plot z-score rho
clear
clc
close all
 phiInitial = [0.11,0.20,0.23,0.24,0.25,0.27,0.29,0.31];
% phiInitial = 0.31;
indexMaximum = length(phiInitial);
for plotIndex =1:indexMaximum
color = char('b','r','k','g','m','c','y','b','r');
DataDirectory = 'C:\Users\lr546\Desktop\125particle 1.397gap\';
interfaceDirectory = [DataDirectory,'phi',num2str(phiInitial(plotIndex)*100),'\interface.xls'];
rsDirectory = [DataDirectory,'phi',num2str(phiInitial(plotIndex)*100),'\rsPlot5-9.png'];
[~,SheetName,~]=xlsfinfo(interfaceDirectory);
for sheetIndex = 5:1:length(SheetName);
      [rsData,~,~] = xlsread(interfaceDirectory,sheetIndex);
      rsData([1,2],:)=[];
      plot(rsData(:,1),rsData(:,3),color(sheetIndex))
      hold on
end
hold off
print('-dpng',rsDirectory,'-r100')
end
%% test calculate outer perimeter 
clear
clc
close all
inletRowPosition = 459.0964;
inletColumnPosition = 601.386;
imageIntensity=imread(fullfile('C:\Users\lr546\Desktop\325particle 1.397gap\phi27\Color Image\315.png'));
[imageHeight,imageWidth,~] = size(imageIntensity); 
imageIntensity=rgb2gray(imageIntensity);
 figure
 imshow(imageIntensity)

BW1 = imbinarize(imageIntensity);
figure
imshow(BW1);
imageInterface = bwperim(BW1);
figure
imshow(imageInterface)
test = edge(BW1,'canny');
figure
imshow(test)
test2 = bwareaopen(test,100);
figure
imshow(test2)
kNum = 1;
for rowNum=1:imageHeight
    for columnNum=1:imageWidth
        rowLocation(kNum)=(rowNum-inletRowPosition);  % unit: pixel
        columnLocation(kNum)=(columnNum-inletColumnPosition);
        interfaceVector(kNum)=imageInterface(rowNum,columnNum);
        kNum=kNum+1;
    end
end
[thetaLocation,rhoLocation]=cart2pol(columnLocation,rowLocation); 
indexInterface = interfaceVector ==1;

rhoInterfaceLocation = rhoLocation(indexInterface);
thetaInterfaceLocation = thetaLocation(indexInterface);
interfaceOuter = [rhoInterfaceLocation',thetaInterfaceLocation']

%% Test image process.
clear 
clc
close all

intial_x_crop = 450;  %x-coordinate of the top-left-most point of the cropped image on original image 
intial_y_crop =98;    %y-coordinate of the top-left-most point of the cropped image on original image 
x_max = 1160;   
y_max = 932; 
medianFilter = [3 3];   %input parameter for median filter, usually not changed
cannyvec = [0.4 0.5]; %input parameter for canny method in edge function
sigma = sqrt(2);     %input parameter for canny method in edge function
denom = 3;    %the parameter determines the threshold with which the code will remove unwanted points so
imageOriginal=imread('test3.png');
x_pos = x_max/100;
y_pos = y_max/100;
savevec = [0 0 x_pos y_pos];
ratio=110;   % pixel per inch.
phi=0.31;
cropvec = [intial_x_crop,intial_y_crop,x_max,y_max];

imageGray = rgb2gray(imageOriginal);
imageCrop = imcrop(imageGray,cropvec);    
im_med = medfilt2(imageCrop,medianFilter);    
canny=edge(im_med,'canny',cannyvec,sigma);
test=bwareaopen(canny,100);
[row,col]=find(test);    
len2 = length(row);
test=bwareaopen(canny,floor(len2/denom));
clear row;
clear col;
[row,col]=find(test);
matr_cart=[col,row];
[cen_fit, radius_fit, ~]=fitcircle(matr_cart, 'tol', 1e-5);

zz = length(row);
radius = zeros(1,zz);
center_y = sum(col)/zz;
center_x = sum(row)/zz;
% imshow(imageCrop);
% hold on
% plot(col,row,'b*','markersize',3)
% set(gcf,'PaperUnits','inches','PaperPosition',savevec)
% fname = sprintf('finger_1.png');
% fullfileName=fullfile(dir3,fname);
% print('-dpng',fullfileName,'-r100');    
for i = 1:zz
    radius(i) = sqrt((col(i)-center_y )^2+(row(i)-center_x)^2);
end
radius_ave = sum(radius)/zz;
% ratio = dia_ave/(5/32); 


len1 = length(row);
min_col = min(col);
max_col = max(col);

% figure
% imshow(imageCrop)



%% Obtain intensity
clear 
clc
close all

Imin=165.5575;

intial_x_crop1 = 405;
intial_y_crop1 = 90;
x_max1 = 1159;
y_max1 = 931;
cropvec1 = [intial_x_crop1 intial_y_crop1 x_max1 y_max1];
%  imageIntensity=imageCrop;

%%%%
phi=0.31;
ratio=111.8493;
center_x=469.7894;
center_y=589.8303;
imageIntensity=imread(fullfile('C:\Users\lr546\Desktop\concentration\test3.png'));
imageIntensity=rgb2gray(imageIntensity);
 figure
 imshow(imageIntensity)
%%%%

BW1 = imbinarize(imageIntensity);
indexS=find(BW1<1);
 imshow(BW1)
BW1 = ~BW1;
BW2 = bwfill(BW1, 'holes');
 figure 
imshow(BW2)
 index=find(BW1<1);

  imageProcess=imageIntensity;
  imageProcess(index)=255;
   figure
   imshow(imageProcess)

image_ave = IntCrct(cropvec1);
 figure
 imshow(image_ave)

ref_corrected1 = double(imageIntensity)./double(image_ave);
ref_corrected2 = uint8(ref_corrected1*180);
figure
imshow(ref_corrected2)

imageCal=double(ref_corrected2);
imageCal=log(imageCal/Imin);
imageCal(index)=0;
kNum=1;

for iNum=1:933
    for jNum=1:1161
        comp=sqrt((iNum-center_x)^2+(jNum-center_y)^2);
        if comp<= 25
            index1(1,kNum)=iNum;
            index1(2,kNum)=jNum;
            kNum=kNum+1;
        end
    end
end
indxi=index1(1,:);
indxj=index1(2,:);
indexAnalysis=sub2ind(size(imageCal),indxi,indxj);
imageCal(indexAnalysis)=0;
indexTemp=find(imageCal ~= 0);
imageCaltemp=zeros(size(imageCal));
imageCaltemp(indexTemp)=1;
areaTemp=sum(imageCaltemp(:));
pixelArea=(1*0.0254/ratio)^2;
k=phi*areaTemp/sum(sum(imageCal(:)))
imageCon=k*imageCal;

ref_corrected2(indexAnalysis)=255;
figure
imshow(ref_corrected2)

[conPlotx,conPloty]=meshgrid(1:1160,1:932);

figure
pcolor(conPlotx,conPloty,imageCon)
axis([400,750,300,650])

kNum=1;
for iNum=1:932
    for jNum=1:1160
        yLocation(kNum)=(iNum-center_x)/ratio*0.0254;  % unit: m
        xLocation(kNum)=(jNum-center_y)/ratio*0.0254;
        concentration(kNum)=imageCon(iNum,jNum);
        kNum=kNum+1;
    end
end

[thetaLocation,rLocation]=cart2pol(xLocation,yLocation);
limitPixel=30:5:150;
limit=limitPixel/ratio*0.0254;


for logicNum=1:length(limit)-1
      indexLogical=logical(rLocation>limit(logicNum) & rLocation<=limit(logicNum+1));
      concentrationMatrixIndex=find(indexLogical );
      concentrationRing=concentration(concentrationMatrixIndex);
      concentrationAverage(logicNum)=mean(concentrationRing);
      locationReal(logicNum)=(limit(logicNum)+limit(logicNum+1))/2;
end

figure
plot(locationReal,concentrationAverage,'*-')

figure
plot3(xLocation,yLocation,concentration)
% axis([-0.0335,0.0335,-0.0335,0.0335,0.1,0.35])




