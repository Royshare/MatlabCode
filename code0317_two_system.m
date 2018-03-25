% test code for abtaining two system miscible fingering 
% by Rui Luo 3/18/2018

clear
clc
close all

particleSize = 325;  % unit: um
gapThickness = 0.406; % unit: mm
FlowRate = 10;      % unit: ml/min

MainDirectory = 'D:\! Backup\Data Using\';

phiInitial = [0.25];
timeFrame = linspace(30,960,32);
% timeFrame = [400];
pixelLowerLimitDelete = 900;
[~,phiIndexTotal] = size(phiInitial);
indexTimeTotal = length(timeFrame);
for phiIndex = 1:phiIndexTotal
      DataDirectory = [MainDirectory,num2str(particleSize),'particle ',...
      num2str(gapThickness),'gap\phi',num2str(phiInitial(phiIndex)*100)];
      pictureDirectory = [DataDirectory,'\Gray Image'];
      fileDirectory = [pictureDirectory,'\1.png'];
      imageGrayOriginal = imread(fileDirectory);
      parameterDirectory = [DataDirectory,'\data.xls'];
      imageFingerDirectory = [DataDirectory,'\Finger Image'];
      mkdir(imageFingerDirectory)
      [imageHeight,imageWidth,~] = size(imageGrayOriginal);
      inletRowPosition = xlsread(parameterDirectory,1,'B2');
      inletColumnPosition = xlsread(parameterDirectory,1,'C2');
      kNum=1;
      for rowNum=1:imageHeight
            for columnNum=1:imageWidth
                  distanceToCenter = sqrt((rowNum-inletRowPosition)^2+...
                           (columnNum-inletColumnPosition)^2);
                  if distanceToCenter > pixelLowerLimitDelete
                        indexDelete(1,kNum)=rowNum;
                        indexDelete(2,kNum)=columnNum;
                        kNum=kNum+1;
                  end
            end
      end
      indexRowDelete=indexDelete(1,:);
      indexColumnDelete=indexDelete(2,:);
      indexDeleteArea=sub2ind([imageHeight,imageWidth],indexRowDelete,indexColumnDelete);
      for indexTime = 1:indexTimeTotal
            fileDirectory = [pictureDirectory,'\',num2str(timeFrame(indexTime)),'.png'];
            imageGrayOriginal = imread(fileDirectory);
%             imfinfo(fileDirectory)
            imageGray = rgb2gray(imageGrayOriginal);
            imageComplement = imcomplement(imageGray);
%             figure; imshow(imageGray)

            %%%%% ????
%             figure; imshow(imageComplement)
%             test = im2uint8(mat2gray(log(1+double(imageComplement))));
%             figure; imshow(test)  

            %%%% ??????
%             figure; imhist(imageGray)
            ylim('auto')
            test = histeq(imageGray);
%             figure, imshow(test)
            %%%% save hist equal image
            [imageSavingHeight,imageSavingWidth] = size(test);
            imageSavingSize =[0,0,imageSavingWidth,imageSavingHeight];
            imageName = sprintf(num2str(timeFrame(indexTime)),'.png');
            imageFingerFileDirectory = fullfile(imageFingerDirectory,imageName);
            imshow(test,'border','tight','initialmagnification','fit');
            set(gcf,'PaperPosition',imageSavingSize/100);
            print('-dpng',imageFingerFileDirectory,'-r100');

            
%             figure, imhist(test)
%             ylim('auto')
            %%%% ??????????
            
            %%%% ????
%             test = medfilt2(test,[1,1]);
            % laplacian ??
%             w = fspecial('laplacian', 0);
%             w = [1,1,1;1,-8,1;1,1,1];
%             g = imfilter(im2double(test), w, 'replicate');
%             figure, imshow(g)
%             figure, imshow(im2double(test)-g)
            
            %%%% ?????

%             fo = imopen(test,se);
%             figure,imshow(fo)
%             foc = imclose(fo,se);
%             figure, imshow(foc)
%             binaryThreshold = graythresh(test);
%             testBW = imbinarize(test,binaryThreshold);
%             figure, imshow(testBW)
%             testBW(indexDeleteArea)=0;
%             imageRemoveBW = bwareaopen(testBW,1000);
%             imageRemoveBW2 = bwareaopen(~imageRemoveBW,1000);
%             figure, imshow(imageRemoveBW)
%             figure, imshow(imageRemoveBW2)
%             imageRemoveBW = imageRemoveBW | ~imageRemoveBW2;
%             figure, imshow(imageRemoveBW)
%             se = strel('diamond',2);
%             gErode = imerode(imageRemoveBW,se);
%             g = imreconstruct(gErode,imageRemoveBW);
%             figure, imshow(g)
%             gfill = imfill(g,'holes');
%             figure, imshow(gfill)
            
%             se1 = strel('square',5);
%             gfillo = imopen(gfill,se1);
%             gfillc = imclose(gfillo,se1);
%             figure, imshow(gfillo)
%             indexTemp = find(gfillo==1);
%             gedge = edge(gfillo,'canny',[0.1,0.15],1.4);
%             figure, imshow(gedge)
%             gedge = bwareaopen(gedge,700);
%             figure, imshow(gedge)
            %%%% ???????
%             gc = ~gfill;
%             D = bwdist(gc);
%             L = watershed(-D);
%             w = L ==0;
%             g2 = gfill & ~w;
%             figure, imshow(g2)
      end
end