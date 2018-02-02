% Main code for processing video to get cropped image and ratio.
% By Rui Luo 2017/12/24
clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables to define.
% target video should be stored in same folder of as where the code locates.
particleSize = 125;  % unit: um
gapThickness = 1.15; % unit: mm
phiInitial = 0.20;
VideoName = 'MVI_0594.MOV';
FlowRate = 150;      % unit: ml/min

FirstFrameIndex = 195;
numberFrameProcess = 660;  % unit: frame
CroppedImageXPosition = 510;  
CroppedImageYPosition = 90; 
CroppedImageWidth = 1160;    
CroppedImageHeight = 920;
MainDirectory = 'C:\Users\lr546\Desktop\';
DataDirectory = [MainDirectory,num2str(particleSize),'particle ',...
      num2str(gapThickness),'gap\phi',num2str(phiInitial*100)];
GrayImageDirectory = fullfile(DataDirectory,'Gray Image');
GrayImageWithEdgeDirectory = fullfile(DataDirectory,'Gray Image with Edge');
ColorImageDirectory = fullfile(DataDirectory,'Color Image');
VideoDirectory = [DataDirectory,'\',VideoName];
mkdir(DataDirectory)
mkdir(GrayImageDirectory)
mkdir(GrayImageWithEdgeDirectory)
mkdir(ColorImageDirectory)
VideoToProcess = VideoReader(VideoDirectory);
LastFrameIndex = FirstFrameIndex+numberFrameProcess;
FrameRate = VideoToProcess.FrameRate;
InletRadius= 5/32*2.54;   % unit: cm   
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

