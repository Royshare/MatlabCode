%%%%%%%%%%%%%%%%%%%%%%%%
% for test concentration plot

clear all;
clc;
close all;

count = 540;
phi = 0.22;
intial_x_crop = 395;
intial_y_crop = 92;

C=180;
x_max = 1159;
y_max = 931;
% axis_vec = [0.12 0.22 0.31 0.35];
titl = char('concentration profile at \phi_0 = 22%');
color = char('b','r','k','g','m','c','y','b','r');
pictureLoc = 'C:\Users\Rui\Desktop\phi22\pics2';
dataLoc = 'C:\Users\Rui\Desktop\phi22\data points';
outputLoc='C:\Users\Rui\Desktop\phi22.xlsx';


cropvec = [intial_x_crop intial_y_crop x_max y_max];
 step = floor(count/90);
 %step=1;
image_ave = IntCrct(cropvec);

fname = sprintf('finger_%i.png',1);
fnamet=sprintf('finger_%i.mat',1);
ref = imread(fullfile(pictureLoc,fname));
refData=load(fullfile(dataLoc,fnamet));
ref=rgb2gray(ref);
ref_corrected = double(ref)./double(image_ave);
ref_corrected = uint8(ref_corrected*C);
ratio=refData.ratio;
figure
% imshow(ref)
imshow(ref_corrected)
[x1,y1] = ginput(4);
ref_crop1 = imcrop(ref_corrected,[x1(1) y1(1) x1(4)-x1(1) y1(4)-y1(1)]);
ref_in1 = mean2(ref_crop1);

[x2,y2] = ginput(4);
ref_crop2 = imcrop(ref_corrected,[x2(1) y2(1) x2(4)-x2(1) y2(4)-y2(1)]);
ref_in2 = mean2(ref_crop2);
clear savdir fname
close all

%fitting process
    fname = sprintf('finger_%i.mat',30);
    fig = load(fullfile(dataLoc,fname)); 
    
    fname2 = sprintf('finger_%i.png',30);
    I = imread(fullfile(pictureLoc,fname2));
    I = rgb2gray(I);
    I_corrected = double(I)./double(image_ave);
    I_corrected = uint8(I_corrected*C);
    m=1;
    for r = 10/25.4*ratio:10:(fig.radius_fit-10)
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
        [col,row] = size(ave_c);
        ave_cc(m) = sum(ave_c(m,:))/row;
        dis(m) = r+5;
        f_log = inline('log(I/I_min)/log(I_max/I_min)','I','I_min','I_max');
        concen(m) = f_log(ave_cc(m),ref_in2,ref_in1);        
        m = m+1;
    end  
min_int = ave_cc(1);
fitting_int = min_int;

for i = 1:step
    [a,concentration,location,error] = CntrPlot(C,i,image_ave,dataLoc,pictureLoc,color,ref_in1,ref_in2,fitting_int,phi);
    hold on
    clear x1 y1 I I_corrected
    sheetName=['Sheet',num2str(i)];
    xlswrite(outputLoc,location',sheetName,'A1');
    xlswrite(outputLoc,concentration',sheetName,'B1');
    xlswrite(outputLoc,error',sheetName,'C1')
end

% axis(axis_vec)
xlabel('R/m','FontSize',15)
ylabel('\Phi','FontSize',15)
set(gca,'FontSize',15)
title(titl)