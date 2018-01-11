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
indexCal=sub2ind(size(imageCal),indxi,indxj);
imageCal(indexCal)=0;
indexTemp=find(imageCal ~= 0);
imageCaltemp=zeros(size(imageCal));
imageCaltemp(indexTemp)=1;
areaTemp=sum(imageCaltemp(:));
pixelArea=(1*0.0254/ratio)^2;
k=phi*areaTemp/sum(sum(imageCal(:)))
imageCon=k*imageCal;

ref_corrected2(indexCal)=255;
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




