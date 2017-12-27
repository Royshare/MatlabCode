% Main code for processing video to get cropped image and ratio.
% By Rui Luo 2017/12/24
clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables to define.
% target video should be stored in same folder of as where the code locates.
VideoToProcess = VideoReader('MVI_0699.MOV');   
FirstFrameIndex = 197;
LastFrameIndex = 780;
InletRadius= 5/32*2.54;   % unit: cm
CroppedImageXPosition = 450;  
CroppedImageYPosition = 98; 
CroppedImageWidth = 1160;    
CroppedImageHeight = 932;    
DataDirectory = 'C:\Users\Rui\Desktop\phi31';
GrayImageDirectory = fullfile(DataDirectory,'Gray Image');
GrayImageWithEdgeDirectory = fullfile(DataDirectory,'Gray Image with Edge');
ColorImageDirectory = fullfile(DataDirectory,'Color Image');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SavedImageWidth = CroppedImageWidth/100;
%SavedImageHeight = CroppedImageHeight/100;
CropVector = [CroppedImageXPosition,CroppedImageYPosition,...
                    CroppedImageWidth,CroppedImageHeight];
%SaveVector = [0,0,SavedImageWidth,SavedImageHeight];

% find ratio between real and image dimension. 
[InletRowPosition,InletColumnPostion,InletImageDiameter] = ...
                    findInletCenter(VideoToProcess,FirstFrameIndex,...
                 CropVector,GrayImageWithEdgeDirectory,ColorImageDirectory);
Ratio = InletImageDiameter/InletRadius;     % unit: pixel/cm
% crop and save every frame
imageProcessIndex=1;
for indexToProcessFrame = FirstFrameIndex:1:LastFrameIndex
    [imageGrayCrop,imageColorCrop] = cropImage(VideoToProcess,...
                              indexToProcessFrame,DataDirectory,CropVector);
    imageName = sprintf(num2str(imageProcessIndex),'.png');
    imshow(imageGrayCrop);
    imageGrayCropDirectory = fullfile(GrayImageDirectory,imageName);
    print('-dpng',imageGrayCropDirectory);
    imageColorCropDirectory = fullfile(ColorImageDirectory,imageName);
    print('-dpng',imageColorCropDirectory);
    imageProcessIndex = imageProcessIndex+1;
end

% save data useful for following process.

output1 = ['Ratio(pixel/cm)','InletXPosition','InletYPosition',...
                     'InletImageDiameter','starting frame','ending frame'];
output2 = [Ratio,InletRowPosition,InletColumnPostion,InletImageDiameter,...
                                            FirstFrameIndex,LastFrameIndex];
xlsDirectory = fullfile(DataDirectory,'data.xls');
xlswrite(xlsDirectory,output1,'A1');
xlswrite(outputLoc,output2,'A2');

