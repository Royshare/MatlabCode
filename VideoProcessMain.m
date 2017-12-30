% Main code for processing video to get cropped image and ratio.
% By Rui Luo 2017/12/24
clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables to define.
% target video should be stored in same folder of as where the code locates.
phiInitial = 0.27;
VideoToProcess = VideoReader('MVI_0643.MOV');   
FirstFrameIndex = 207;
LastFrameIndex = FirstFrameIndex+450;
CroppedImageXPosition = 450;  
CroppedImageYPosition = 96; 
CroppedImageWidth = 1160;    
CroppedImageHeight = 920;

FrameRate = VideoToProcess.FrameRate;
InletRadius= 5/32*2.54;   % unit: cm   
DataDirectory = ['C:\Users\lr546\Desktop\phi',num2str(phiInitial*100)];
GrayImageDirectory = fullfile(DataDirectory,'Gray Image');
GrayImageWithEdgeDirectory = fullfile(DataDirectory,'Gray Image with Edge');
ColorImageDirectory = fullfile(DataDirectory,'Color Image');
mkdir(DataDirectory)
mkdir(GrayImageDirectory)
mkdir(GrayImageWithEdgeDirectory)
mkdir(ColorImageDirectory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CropVector = [CroppedImageXPosition,CroppedImageYPosition,...
                    CroppedImageWidth,CroppedImageHeight];
ImageSizeVector = [0,0,CroppedImageWidth,CroppedImageHeight];

% 1. find ratio between real and image dimension. 
[InletRowPosition,InletColumnPostion,InletImageDiameter] = ...
   findInletCenter(VideoToProcess,FirstFrameIndex,CropVector,DataDirectory);
Ratio = InletImageDiameter/InletRadius;     % unit: pixel/cm
% 2. crop and save every frame
imageProcessIndex=1;
for indexToProcessFrame = FirstFrameIndex:1:LastFrameIndex
    [imageGrayCrop,imageColorCrop] = cropImage(VideoToProcess,...
                              indexToProcessFrame,CropVector);
    imageName = sprintf(num2str(imageProcessIndex),'.png');
    imshow(imageGrayCrop,'border','tight','initialmagnification','fit');
    set(gcf,'PaperPosition',ImageSizeVector/100);
    imageGrayCropDirectory = fullfile(GrayImageDirectory,imageName);
    print('-dpng',imageGrayCropDirectory,'-r100');
    imshow(imageColorCrop,'border','tight','initialmagnification','fit');
    imageColorCropDirectory = fullfile(ColorImageDirectory,imageName);
    set(gcf,'PaperPosition',ImageSizeVector/100);
    print('-dpng',imageColorCropDirectory,'-r100');
    imageProcessIndex = imageProcessIndex+1;
end

% 3. save data useful for following process.
output1 = {'Ratio(pixel/cm)','InletRowPosition(pixel)','InletColumnPosition(pixel)',...
'InletImageDiameter(pixel)','starting frame','ending frame','FrameRate(/s)'};
output2 = [Ratio,InletRowPosition,InletColumnPostion,InletImageDiameter,...
                                  FirstFrameIndex,LastFrameIndex,FrameRate];
xlsDirectory = fullfile(DataDirectory,'data.xls');
xlswrite(xlsDirectory,output1,1,'A1');
xlswrite(xlsDirectory,output2,1,'A2');

