% this function is used for calculate interfacial perimeter between suspension
% and air.
% By Rui Luo 2018/1/9

function interfaceOuter = getOuterInterfacePosition(imageIntensity,inletRowPosition,inletColumnPosition)

      [imageHeight,imageWidth,~] = size(imageIntensity); 
      imageIntensity=rgb2gray(imageIntensity);

      imageBinary = imbinarize(imageIntensity);
      imageEdgeTemp = edge(imageBinary,'canny');
      
      imageEdge = bwareaopen(imageEdgeTemp,190);
      figure
      imshow(imageEdge)
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