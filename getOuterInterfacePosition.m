% this function is used for calculate interfacial perimeter between suspension
% and air.
% By Rui Luo 2018/1/9

function interfaceOuter = getOuterInterfacePosition(imageIntensity,inletRowPosition,inletColumnPosition)

      [imageHeight,imageWidth,~] = size(imageIntensity); 
      imageIntensity=rgb2gray(imageIntensity); 
%       imageIntensity = imadjust(imageIntensity);
% figure, imshow(imageIntensity)
%       imageIntensity = histeq(imageIntensity);
%       imageIntensity = adapthisteq(imageIntensity);
% figure
% imshow(imageIntensity)      
%       imageBinary = imbinarize(imageIntensity);
% figure
% imshow(imageBinary)
      
%       imageEdgeTemp = edge(imageBinary,'canny');
imageBinary = imbinarize(imageIntensity);
      [~,threshold] = edge(imageBinary,'sobel');
      imageEdgeTemp = edge(imageBinary,'sobel',0.5*threshold);
% figure, imshow(imageEdgeTemp)      
%       imageEdge = bwareaopen(imageEdgeTemp,190);
% figure
% imshow(imageEdge)
      se90 = strel('line', 3, 90);
      se0 = strel('line', 3, 0);
      BWsdil = imdilate(imageEdgeTemp, [se90 se0]);
% figure, imshow(BWsdil), title('dilated gradient mask');
      BWdfill = imfill(BWsdil, 'holes');
      BWnobord = imclearborder(BWdfill, 4);
% figure, imshow(BWnobord);
      seD = strel('diamond',1);
BWfinal = imerode(BWnobord,seD);
% BWfinal = imerode(BWfinal,seD);
% figure, imshow(BWfinal), title('segmented image');
BWfinal = bwareaopen(BWfinal,300);
% BWoutline = bwperim(BWfinal);
imageEdge = edge(BWfinal,'canny');
% figure, imshow(imageEdge);
% Segout = imageIntensity; 
% Segout(imageEdge) = 255; 
% figure, imshow(Segout), title('outlined original image');


      kNum = 1;
      rowLocation = zeros(1,imageHeight*imageWidth);
      columnLocation = zeros(1,imageHeight*imageWidth);
      interfaceVector = zeros(1,imageHeight*imageWidth);
      for rowNum=1:imageHeight
            for columnNum=1:imageWidth
                  rowLocation(kNum)=(rowNum-inletRowPosition);  % unit: pixel
                  columnLocation(kNum)=(columnNum-inletColumnPosition);
                  interfaceVector(kNum)=imageEdge(rowNum,columnNum);
                  kNum=kNum+1;
            end
      end
      [thetaLocation,rhoLocation]=cart2pol(columnLocation,rowLocation); 
      indexInterface = interfaceVector ==1;

      rhoInterfaceLocation = rhoLocation(indexInterface);
      thetaInterfaceLocation = thetaLocation(indexInterface);
      rhoNormalized = zscore(rhoInterfaceLocation);
      interfaceOuter = [thetaInterfaceLocation',rhoInterfaceLocation',rhoNormalized'];
      interfaceOuter = sortrows(interfaceOuter,1);
end