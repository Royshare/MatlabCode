



%% reference
clear
clc
close all

ref_cal = VideoReader('reference_correction.MOV');
indexStartRead = 1;
indexEndRead = 50;
indexTotalRead = indexEndRead-indexStartRead+1;

indexImageData = 1;
imageOriginal = zeros(1,indexTotalRead);
for indexRead = indexStartRead:1:indexEndRead
imageOriginal(indexImageData) = read(ref_cal,indexRead);
imageGray = rgb2gray(imageOriginal(indexImageData));
end
imageAverageValue = mean(double(imageGray));
imageAverageBit = uint8(imageAverageValue)



%% temp reference inlet center
clear
clc
close all
video = VideoReader('reference_correction.MOV');
medianFilterThreshold = [5,5];
pixelUpperLimit=100;
cannyThreshold=[0.4,0.5];
cannySigma= 28;

% ImageSizeVector = [0,0,CropVector(3),CropVector(4)];
% GrayImageWithEdgeDirectory = fullfile(DataDirectory,'Gray Image with Edge');
% ColorImageDirectory = fullfile(DataDirectory,'Color Image');

imageColorOriginal = imread('C:\Users\lr546\Desktop\reference\Gray Image\1.png');
imageGrayOrignial = rgb2gray(imageColorOriginal);
% imageColorCrop = imcrop(imageColorOriginal,CropVector);
% imageGrayCrop = imcrop(imageGrayOrignial,CropVector);
imageMedianFiltered = medfilt2(imageGrayOrignial,medianFilterThreshold);
imageCannyEdgeDetected = edge(imageMedianFiltered,'canny',cannyThreshold,cannySigma);
imageRemovedSmallObject = bwareaopen(imageCannyEdgeDetected,pixelUpperLimit);
[rowEdge,columnEdge] = find(imageRemovedSmallObject);
rowCenter = mean(rowEdge);
columnCenter = mean(columnEdge);
temp = (rowEdge-rowCenter).^2+(columnEdge-columnCenter).^2;
inletImageRadius = mean(sqrt(temp));

imshow(imageGrayOrignial,'border','tight','initialmagnification','fit');
hold on
plot(columnEdge,rowEdge,'b*','markersize',1)
% set(gcf,'PaperPosition',ImageSizeVector/100)

% fname = sprintf('0.png');
% GrayImageFullDirectory = fullfile(GrayImageWithEdgeDirectory,fname);
% print('-dpng',GrayImageFullDirectory,'-r100'); 

hold off
figure
imshow(imageColorOriginal,'border','tight','initialmagnification','fit');
% set(gcf,'PaperPosition',ImageSizeVector/100)
% fname = sprintf('0.png');
% ColorImageFullDirectory = fullfile(ColorImageDirectory,fname);
% print('-dpng',ColorImageFullDirectory,'-r100');


%% crop reference image
% Main code for processing video to get cropped image and ratio.
% By Rui Luo 2017/12/24
clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables to define.
% target video should be stored in same folder of as where the code locates.
VideoToProcess = VideoReader('reference_correction.MOV');   
FirstFrameIndex = 14;
LastFrameIndex = FirstFrameIndex+59;
% LastFrameIndex = FirstFrameIndex+5;   % for test to decide first frame
CroppedImageXPosition = 390;  
CroppedImageYPosition = 83; 
CroppedImageWidth = 1200;    
CroppedImageHeight = 942;

FrameRate = VideoToProcess.FrameRate;
InletRadius= 5/32*2.54;   % unit: cm   
DataDirectory = 'C:\Users\lr546\Desktop\reference';
GrayImageDirectory = fullfile(DataDirectory,'Gray Image');
GrayImageWithEdgeDirectory = fullfile(DataDirectory,'Gray Image with Edge');
ColorImageDirectory = fullfile(DataDirectory,'Color Image');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CropVector = [CroppedImageXPosition,CroppedImageYPosition,...
                    CroppedImageWidth,CroppedImageHeight];
ImageSizeVector = [0,0,CroppedImageWidth,CroppedImageHeight];


% 2. crop and save every frame
imageProcessIndex=1;
for indexToProcessFrame = FirstFrameIndex:1:LastFrameIndex
    [imageGrayCrop,imageColorCrop] = cropImage(VideoToProcess,...
                              indexToProcessFrame,CropVector);
    imageName = sprintf(num2str(imageProcessIndex),'.png');
    imshow(imageGrayCrop,'border','tight','initialmagnification','fit');
    set(gcf,'PaperPosition',ImageSizeVector/100);
    imageGrayCropDirectory = fullfile(GrayImageDirectory,imageName);
    print('-dpng',imageGrayCropDirectory);
    imshow(imageColorCrop,'border','tight','initialmagnification','fit');
    imageColorCropDirectory = fullfile(ColorImageDirectory,imageName);
    set(gcf,'PaperPosition',ImageSizeVector/100);
    print('-dpng',imageColorCropDirectory);
    imageProcessIndex = imageProcessIndex+1;
end






