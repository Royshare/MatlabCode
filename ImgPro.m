%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to process pictures of specific frames of video 
function ImgPro(OBJ,starting_frame,vid_Height,vid_Width,frmlocator,cropvec,medvec,cannyvec,sigma,dia_ave,ratio,savevec,dir1,dir2,dir3,denom)
    % Initialize variable mov, mov_gray
    mov = struct('cdata',zeros(vid_Height,vid_Width,3,'uint8'),'colormap',[]);
    mov_gray = struct('cdata',zeros(vid_Height,vid_Width,3,'uint8'),'colormap',[]);

    mov.cdata=read(OBJ,starting_frame+frmlocator);
    mov_gray.cdata = rgb2gray(mov.cdata);
    image = mov_gray.cdata; 
    imgcrop = imcrop(image,cropvec);    
    im_med = medfilt2(imgcrop,medvec);
    
    canny=edge(im_med,'canny',cannyvec,sigma);
    test=bwareaopen(canny,100);
    [row,col]=find(test);    
    len2 = length(row);
    test=bwareaopen(canny,floor(len2/denom));
    clear row;
    clear col;
    [row,col]=find(test);
    matr_cart=[col,row];
   
   [cen_fit, radius_fit, ~]=fitcircle(matr_cart, 'tol', 1e-5);
    circ = 2*pi*radius_fit;
     
    len1 = length(row);
    min_col = min(col);
    max_col = max(col);
   
    i = min_col;
    angle1 = zeros(1,len1);
    dist = zeros(1,len1);
    while (i >= min_col & i <= max_col)
        k = 1;
        for j = 1:len1
            if col(j) == i
                col_temp(k) = col(j);
                row_temp(k) = row(j);
                coor_temp(k) = j;
                k = k+1;
            end
        end
        k = k-1;
        tan_temp = zeros(1,k);
        angle_temp = zeros(1,k);
        if i <= cen_fit(1)            
            for n = 1:k
                tan_temp(n) = (row_temp(n)-cen_fit(2))/(col_temp(n)-cen_fit(1));
                angle_temp(n) = atan(tan_temp(n))/pi*180+180;
                angle1(coor_temp(n)) = angle_temp(n);
                dist(coor_temp(n)) = sqrt((row(coor_temp(n))-cen_fit(2))^2+(col(coor_temp(n))-cen_fit(1))^2)/ratio;
                if dist(coor_temp(n)) < (dia_ave+5)/ratio
                    angle1(coor_temp(n)) = 0;
                    dist(coor_temp(n)) = 0;
                end
            end
        else             
            for n = 1:k
                tan_temp(n) = (row_temp(n)-cen_fit(2))/(col_temp(n)-cen_fit(1));
                if row_temp(n) <= cen_fit(2)
                    angle_temp(n) = atan(tan_temp(n))/pi*180+360;
                else
                    angle_temp(n) = atan(tan_temp(n))/pi*180;
                end
                angle1(coor_temp(n)) = angle_temp(n);
                dist(coor_temp(n)) = sqrt((row(coor_temp(n))-cen_fit(2))^2+(col(coor_temp(n))-cen_fit(1))^2)/ratio;
                if dist(coor_temp(n)) < (dia_ave+5)/ratio
                    angle1(coor_temp(n)) = 0;
                    dist(coor_temp(n)) = 0;
                end
            end            
        end
        i = i+1;        
    end
    
    A1=[angle1;dist];
    A2=[angle1;row'];
    A3=[angle1;col'];
    [Y,I]=sort(A1(1,:));
    B1=A1(:,I);
    B2=A2(:,I);
    B3=A3(:,I);
    
    savdir1 = dir1;
    fname = sprintf('finger_%i.mat',frmlocator);
    save(fullfile(savdir1,fname),'B1','B2','B3','cen_fit','radius_fit','ratio')

    h=figure(1);
    set(h,'PaperPositionMode','auto')
    clf
    
    imshow(imgcrop);
    set(gcf,'PaperUnits','inches','PaperPosition',savevec)
    fname = sprintf('finger_%i.png',frmlocator);
    savdir2 = dir2;
    fullfileName=fullfile(savdir2,fname);
    print('-dpng',fullfileName,'-r100');
    
    hold on
    plot(col,row,'b*','markersize',3)
    hold on
    ang = 0:0.01:2*pi;
    xp = radius_fit*cos(ang)+cen_fit(1);
    yp = radius_fit*sin(ang)+cen_fit(2);
    plot(xp,yp,'r-','markersize',3)
    
    
    set(gcf,'PaperUnits','inches','PaperPosition',savevec)
    fname = sprintf('finger_%i.png',frmlocator);
    savdir2 = dir3;
    fullfileName=fullfile(savdir2,fname);
    print('-dpng',fullfileName,'-r100');
    close all
end