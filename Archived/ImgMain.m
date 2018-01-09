%% Section I Raw Image Processing
clear all
close all
clc
readerobj = VideoReader('MVI_0295.MOV');   %load the video with its name. target video should be stored in same folder of as where the code locates.
vid_Height = readerobj.Height; 
vid_Width = readerobj.Width;

% parameters need to be changed
starting_frame = 242;
ending_frame = 400;
intial_x_crop = 440;  %x-coordinate of the top-left-most point of the cropped image on original image 
intial_y_crop =102;    %y-coordinate of the top-left-most point of the cropped image on original image 
x_max = 1150;    %size of cropped image (x)
y_max = 920;     %size of cropped image (y)
medvec = [3 3];   %input parameter for median filter, usually not changed
cannyvec = [0.3 0.5]; %input parameter for canny method in edge function
sigma = sqrt(2);     %input parameter for canny method in edge function
denom = 3;    %the parameter determines the threshold with which the code will remove unwanted points so as to keep points on edge only.
dir1 = 'C:\Users\Rui\Desktop\data points';%local path of directory to store data
dir2 = 'C:\Users\Rui\Desktop\pics2';%local path of directory to store images w/o plotted edge
dir3 = 'C:\Users\Rui\Desktop\pics';%local path of directory to store images w/ plotted edge

% input definition
numFrames = ending_frame-starting_frame;
mov(1:numFrames) = struct('cdata',zeros(vid_Height,vid_Width,3,'uint8'),'colormap',[]);
mov_gray(1:numFrames) = struct('cdata',zeros(vid_Height,vid_Width,3,'uint8'),'colormap',[]);
x_pos = x_max/100;
y_pos = y_max/100;
cropvec = [intial_x_crop intial_y_crop x_max y_max];
savevec = [0 0 x_pos y_pos];

% determine the center of the hole and an average diameter
mov(1).cdata=read(readerobj,starting_frame+1);
mov_gray(1).cdata = rgb2gray(mov(1).cdata);
image = mov_gray(1).cdata;
imgcrop = imcrop(image,cropvec);
im_med = medfilt2(imgcrop,[15 15]);
canny=edge(im_med,'canny',[0.1 0.2],5);
test=bwareaopen(canny,100);
[row,col]=find(test);
zz = length(row);
dia = zeros(1,zz);
center_x = sum(col)/zz;
center_y = sum(row)/zz;
imshow(imgcrop);
hold on
plot(col,row,'b*','markersize',3)
set(gcf,'PaperUnits','inches','PaperPosition',savevec)
fname = sprintf('finger_1.png');
fullfileName=fullfile(dir3,fname);
print('-dpng',fullfileName,'-r100');    
for i = 1:zz
    dia(i) = sqrt((col(i)-center_x)^2+(row(i)-center_y)^2);
end
dia_ave = sum(dia)/zz;
ratio = dia_ave/(5/32); 

% whole processing
count = 0;
for z = 1:numFrames  %
    ImgPro(readerobj,starting_frame,vid_Height,vid_Width,z,cropvec,medvec,cannyvec,sigma,dia_ave,ratio,savevec,dir1,dir2,dir3,denom);
    count = count+1;
end
count

%% Section II Concentration Plot
clear all;
clc;
close all;
count = 620;
C=180;
phi = 0.15;
intial_x_crop = 400;
intial_y_crop = 120;
x_max = 1149;
y_max = 919;
% axis_vec = [0 1.2 0.3 0.38];
titl = char('concentration profile at \phi_0 = 32%');
color = char('b','r','k','g','m','c','y','b','r');
dir1 = 'C:\Users\Rui\Desktop\pics2';
dir2 = 'C:\Users\Rui\Desktop\data points';

cropvec = [intial_x_crop intial_y_crop x_max y_max];
step = floor(count/90);
image_ave = IntCrct(cropvec);
ratio = 110;

fname = sprintf('finger_%i.png',1);
ref = imread(fullfile(dir1,fname));
ref=rgb2gray(ref);
ref_corrected = double(ref)./double(image_ave);
ref_corrected = uint8(ref_corrected*C);
figure
% imshow(ref)
imshow(ref_corrected)
[x1 y1] = ginput(4);
ref_crop1 = imcrop(ref_corrected,[x1(1) y1(1) x1(4)-x1(1) y1(4)-y1(1)]);
ref_in1 = mean2(ref_crop1);

[x2 y2] = ginput(4);
ref_crop2 = imcrop(ref_corrected,[x2(1) y2(1) x2(4)-x2(1) y2(4)-y2(1)]);
ref_in2 = mean2(ref_crop2);
clear savdir fname
close all

%fitting process
    fname = sprintf('finger_%i.mat',61);
    fig = load(fullfile(dir2,fname)); 
    
    fname2 = sprintf('finger_%i.png',61);
    I = imread(fullfile(dir1,fname2));
    I = rgb2gray(I);
    I_corrected = double(I)./double(image_ave);
    I_corrected = uint8(I_corrected*C);
    m=1;
    for r = 15/25.4*ratio:10:(fig.radius_fit-10)
        n=1; 
        for rr = r:0.5:min((r+10),fig.radius_fit-10)
            k = 1;
            for j = 0:6:360
                xi(n,k) = floor(rr*cos(j*pi/180)+fig.cen_fit(1));
                yi(n,k) = floor(rr*sin(j*pi/180)+fig.cen_fit(2));
                k = k+1;
            end
            for ii = 1:length(xi)
                c(ii)=I_corrected(yi(n,ii),xi(n,ii));
%                 c(ii)=I(yi(n,ii),xi(n,ii));
            end
            
            ave_c(m,n) = sum((c))/length(xi);
            n = n+1;
        end
        [col row] = size(ave_c);
        ave_cc(m) = sum(ave_c(m,:))/row;
        dis(m) = r+5;
        f_log = inline('log(I/I_min)/log(I_max/I_min)','I','I_min','I_max');
        concen(m) = f_log(ave_cc(m),ref_in2,ref_in1);        
        m = m+1;
    end  
min_int = ave_cc(1);
fitting_int = min_int;

for i = 1:1
    [a b c] = CntrPlot(C,i,image_ave,dir2,dir1,color,ref_in1,ref_in2,fitting_int,phi);
    hold on
    clear x1 y1 I I_corrected
end

% axis(axis_vec)
xlabel('R','FontSize',15)
ylabel('\Phi','FontSize',15)
set(gca,'FontSize',15)
title(titl)

%% Section III Intensity Plot
clear all
close all
clc
count = 610;
C=180;
phi = 0.28;
intial_x_crop = 440;
intial_y_crop = 102;
x_max = 1160;
y_max = 932;
color = char('b','r','k','g','m','c','y','b','r');
dir1 = 'C:\Users\Rui\Desktop\pics2';
dir2 = 'C:\Users\Rui\Desktop\data points';

cropvec = [intial_x_crop intial_y_crop x_max y_max];
step = floor(count/90);
image_ave = IntCrct(cropvec);
ratio = 110;

fname = sprintf('finger_%i.png',600);
ref = imread(fullfile(dir1,fname));
ref=rgb2gray(ref);
ref_corrected = double(ref)./double(image_ave);
ref_corrected = uint8(ref_corrected*C);
figure
imshow(ref_corrected)
[x1 y1] = ginput(4);
ref_crop1 = imcrop(ref,[x1(1) y1(1) x1(4)-x1(1) y1(4)-y1(1)]);
ref_in1 = mean2(ref_crop1);

[x2 y2] = ginput(4);
ref_crop2 = imcrop(ref,[x2(1) y2(1) x2(4)-x2(1) y2(4)-y2(1)]);
ref_in2 = mean2(ref_crop2);
clear savdir fname
close all

%fitting process
    fname = sprintf('finger_%i.mat',90);
    fig = load(fullfile(dir2,fname)); 
    
    fname2 = sprintf('finger_%i.png',90);
    I = imread(fullfile(dir1,fname2));
    I = rgb2gray(I);
    I_corrected = double(I)./double(image_ave);
    I_corrected = uint8(I_corrected*C);
    m=1;
    for r = 20/25.4*ratio:10:(fig.radius_fit-10)
        k=1; 
        for rr = r:0.5:(r+10)
            n = 1;
            for j = 0:6:360
                xi(n,k) = floor(rr*cos(j*pi/180)+fig.cen_fit(1));
                yi(n,k) = floor(rr*sin(j*pi/180)+fig.cen_fit(2));
                k = k+1;
            end
            for ii = 1:length(xi)
%                 c(ii)=I_corrected(yi(n,ii),xi(n,ii));
                c(ii)=I(yi(n,ii),xi(n,ii));
            end
            
            ave_c(m,n) = sum((c))/length(xi);
            n = n+1;
        end
        [col row] = size(ave_c);
        ave_cc(m) = sum(ave_c(m,:))/row;
        dis(m) = r+5;
        f_log = inline('log(I/I_min)/log(I_max/I_min)','I','I_min','I_max');
        concen(m) = f_log(ave_cc(m),ref_in2,ref_in1);        
        m = m+1;
    end
    max_dis = max(dis);
    min_con = concen(1);
    min_int = min(ave_cc);
    fitting_int = min_int;
    
    sum_con = 0;
    for p = 2:length(dis)
        sum_con = concen(p)*((dis(p)-5)^2-(dis(p-1)-5)^2)+sum_con;
    end
    sum_con = sum_con+concen(1)*(dis(1)-5)^2;
    fitting_con = sum_con/(phi*(max_dis-5)^2);


figure;
for i = 1:step
    IntPlot(C,i,image_ave,dir1,dir2,color,ref_in1,ref_in2,fitting_con,fitting_int); 
    hold on    
    clear x1 y1 I I_corrected
end

axis([20 100 20 40])
xlabel('R','FontSize',15)
ylabel('I','FontSize',15)
set(gca,'FontSize',15)
title('intensity profile at \phi_0 = 28%')
% title('Intensity and \Phi profile at 20% initial concentration')
% legend('t=3s','t=6s','t=9s','t=12s','t=15s')


%% calculation for Lambda value
clear all;
clc;
close all;
figure
count = 740;
Q = 150*10^-6/60;
dir1 = 'C:\Users\Rui\Desktop\data points';
LambdaCal(count,Q,dir1)
title('\Lambda of 30%','FontSize',15)
xlabel('tQ/R^2 h','FontSize',15)
ylabel('\Lambda','FontSize',15)
axis([0 5 0 4*10^-4])
set(gca,'FontSize',15)


