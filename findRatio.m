function [xPostion,yPosition,diameter,ratio]=findRatio(video,FirstFrameIndex,inletRadius)

VideoFrameColor(1).cdata=readframe(VideoToProcess,FirstFrameIndex+1);
VideoFrameGray(1).cdata = rgb2gray(VideoFrameColor(1).cdata);
image = VideoFrameGray(1).cdata;
imgcrop = imcrop(image,CropVector);
im_med = medfilt2(imgcrop,[15 15]);
canny=edge(im_med,'canny',[0.4 0.5],5);
test=bwareaopen(canny,100);
[row,col]=find(test);
zz = length(row);
dia = zeros(1,zz);
center_x = sum(col)/zz;
center_y = sum(row)/zz;
imshow(imgcrop);
hold on
plot(col,row,'b*','markersize',3)
set(gcf,'PaperUnits','inches','PaperPosition',SaveVector)
fname = sprintf('finger_1.png');
fullfileName=fullfile(ColorImageDirectory,fname);
print('-dpng',fullfileName,'-r100');    
for i = 1:zz
    dia(i) = sqrt((col(i)-center_x)^2+(row(i)-center_y)^2);
end
dia_ave = sum(dia)/zz;
ratio = dia_ave/(inletRadius); 
end