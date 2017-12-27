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
InletRadius= 5/32;   % unit: inch
CroppedImageXPosition = 450;  
CroppedImageYPosition = 98; 
CroppedImageWidth = 1160;    
CroppedImageHeight = 932;    
MedianFilterThreshold = [3 3];  
CannyVectorThreshold = [0.4 0.5]; 
CannyVectorSigma = sqrt(2);     
denom = 3;    %the parameter determines the threshold with which the code will remove unwanted points so as to keep points on edge only.

IntensityDirectory = 'C:\Users\Rui\Desktop\phi22\Intensity Data';
GreyImageDirectory = 'C:\Users\Rui\Desktop\phi22\Grey Image';
GreyImageWithEdgeDirectory = 'C:\Users\Rui\Desktop\phi22\Grey Image with Edge';
ColorImageDirectory = 'C:\Users\Rui\Desktop\phi22\Color Image';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VideoHeight = VideoToProcess.Height; 
VideoWidth = VideoToProcess.Width;

FrameTotalNumber = LastFrameIndex-FirstFrameIndex+1;

%SavedImageWidth = CroppedImageWidth/100;
%SavedImageHeight = CroppedImageHeight/100;
CropVector = [CroppedImageXPosition,CroppedImageYPosition,...
                    CroppedImageWidth,CroppedImageHeight];
%SaveVector = [0,0,SavedImageWidth,SavedImageHeight];

% find ratio between real and image dimension. 
[InletXPosition,InletYPostion,InletImageDiameter] = ...
                    findInletCenter(VideoToProcess,FirstFrameIndex,...
                    CropVector,GreyImageWithEdgeDirectory,ColorImageDirectory);

% process every frame
