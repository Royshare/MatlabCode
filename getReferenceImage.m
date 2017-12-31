% this function is used for modify original image according to reference.
% By Rui Luo 2017/12/27

function imageReferenceAverageValue = getReferenceImage(imageIntensity,DataDirectory,inletRowPosition,inletColumnPostition)

indexStartRead = 1;
indexEndRead = 30;

% indexTotalRead = indexEndRead-indexStartRead+1;
centerDataDirectory = [DataDirectory,'reference\center.xls'];
inletRowReferencePosition = xlsread(centerDataDirectory,1,'A2');
inletColumnReferencePosition = xlsread(centerDataDirectory,1,'B2');
[ImageSuspensionHeight,ImageSuspensionWidth,~] = size(imageIntensity);

cropVector = [inletColumnReferencePosition-inletColumnPostition,...
                    inletRowReferencePosition-inletRowPosition,...
                        ImageSuspensionWidth-1,ImageSuspensionHeight-1];
for indexRead = indexStartRead:1:indexEndRead
    imageReferenceDirectory = [DataDirectory,'reference\Gray Image\',num2str(indexRead),'.png'];
    imageReferenceOriginal(:,:,:,indexRead) = imread(imageReferenceDirectory);
    imageReferenceCrop(:,:,:,indexRead) = imcrop(imageReferenceOriginal(:,:,:,indexRead),cropVector);
end
imageReferenceAverageValue = mean(double(imageReferenceCrop),4);
imageReferenceAverageValue = uint8(imageReferenceAverageValue);
imageReferenceAverageValue = rgb2gray(imageReferenceAverageValue);
end