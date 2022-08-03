%{
..................Property of Dr. Chong's Lab...............
This code was produced as an algorithm for solving a problem
...........Originally written by Mr. Rifat Zaman............
...............Edited by Md Kamrul Hasan Foysal.............
.............Under Supervision of Dr. Jo Woon Chong.........
..................©Copyright: ALL RIGHTS RESERVED...........
%}
clc;
clearvars;
close all;

addpath(genpath(pwd));

a='1';
A =dir( fullfile(a , '*.png') );

%%
NofFrames=numel(A);
tmfrinit=10;
tmfrmax=NofFrames;
 
tmfr=tmfrinit:tmfrmax;
Fs = 30;
Ts = 1/Fs;
time=0:Ts:(tmfrmax-1)*Ts;
sg=zeros(1,tmfrmax);
sgg=zeros(1,tmfrmax);
count1=0;
count2=0;
count2g=0;
count1g=0;
sgn=zeros(1,tmfrmax);
sggn=zeros(1,tmfrmax);
for ii=tmfrinit:tmfrmax
          files='1\a ('; 
                
    stringname=num2str (ii);
    f1=strcat(files, stringname, ').png');
    
I1 = imread(f1);
I1=imresize(I1, .4);

ii
height = size(I1,1);
width = size(I1,2);

r=height;
c=width;

minarea=width*height/50; % for discarding small white regions during curve detection
    img_reduced_current = I1;

    %% Bit Rearrangment

%     figure(201), imshow( I1);
    Img565 =uint8(zeros(height,width,2));
    
    Img565(:,:,2) = bitand(img_reduced_current(:,:,1),bin2dec('11111000'))...
        +(bitand(img_reduced_current(:,:,2),bin2dec('11100000'))*(2^-5));

    Img565(:,:,1) = bitand(img_reduced_current(:,:,2),bin2dec('00011100'))*(2^3)+...
        (bitand(img_reduced_current(:,:,3),bin2dec('11111000'))*(2^-3));
    
    img_rearranged = uint8(zeros(size(img_reduced_current)));
    
    img_rearranged(:,:,1) = bitand(Img565(:,:,1),bin2dec('11111000'));
    img_rearranged(:,:,2) = bitand(Img565(:,:,1),bin2dec('00000111'))*(2^5)...
        +bitand(Img565(:,:,2),bin2dec('11100000'))*(2^-3);
    img_rearranged(:,:,3) = bitand(Img565(:,:,2),bin2dec('00011111'))*(2^3);  
    figure(200), imshow( img_rearranged);

    %%
    
    x = rgb2gray(img_rearranged);    %convert to greyscale image (not binary)
    
    figure(209), imshow(x);
    G1t = x; % G1t is used later

     figure(202), imshow(x);
     

    % Edge Detection

    mm =30;

    yyy=x;
% 
    for j=2:r-1 % r is no of rows in G1t
        for i=2:c-1  % c is no of columns in G1t
            p=i-1;
            q=i+1;
            a=j-1;
            b=j+1;
            e1=abs(x(j,i)-x(j,p));
            e2=abs(x(j,i)-x(j,q));
            d1=abs(x(j,i)-x(a,i));
            d2=abs(x(j,i)-x(b,i));
            f1=abs(x(j,i)-x(a,p));
            f2=abs(x(j,i)-x(b,q));
            g1=abs(x(j,i)-x(a,q));
            g2=abs(x(j,i)-x(b,p));

            if (d1 >= mm) || (d2 >= mm)
            yyy(j,i)=255;
            elseif (e1 >= mm) || (e2 >= mm)
            yyy(j,i)=255;
            elseif (f1 >= mm) || (f2 >= mm)
            yyy(j,i)=255;
            elseif (g1 >= mm) || (g2 >= mm)
            yyy(j,i)=255;

            else
                yyy(j,i)=0;  

            end
        end
    end
    x=yyy;
%

    se = strel('line',10, 0);
    x=imdilate(x, se);
    x=x(2:end-1, 2:end-1);
    x=imresize(x, [r c]);
%     Y=bwconncomp(x);
%     
%     
%     bL=length(Y.PixelIdxList);
%     for i=1:bL
%         if length(Y.PixelIdxList{1, i})<width.*2
% level = graythresh(x)
% % x = imbinarize(x,level.*2);
%     x=im2bw(G1t, level./2);
%          disp ('XXX')
%          break
%         end
%     end
%         

% 
% level = graythresh(x)
%     x=im2bw(G1t, level./2);


       figure(11110011), imshow(x);

    %% Smoothing
 
    params.bz_thresh = -2; %thresholding for distant image
    params.as_scale = 1; %1/4 default
    params.debug = 0; % show all figures for 1


    B_our = bwsmooth_2( x, params );

    B_our = imcomplement(B_our);  %%% white to black and black to white
    figure(99), imshow(B_our);
    

    %%    
    C1=bwconncomp(B_our);
   
    cL=length(C1.PixelIdxList);
 
    if C1.NumObjects==0
        continue
     end
 
    cid=0;
    cmax=0;
    if ii==tmfr(1) %for the first frame
%%
        varmax=height/2;
        for ci=1:cL
            cL2=length(C1.PixelIdxList{1,ci});

            [jlist,ilist]=ind2sub([r,c],C1.PixelIdxList{1,ci});

            imean=round(mean(ilist));
            jmean=round(mean(jlist));
            varpos=1./abs(imean-height/2)

            if varpos<varmax && cL2>minarea %checking for every white curve
               varmax=varpos
               cid=ci;
               cmax=length(C1.PixelIdxList{1,ci});
            end

        end
        %%
        C_prev=C1;
        cid_prev=cid;
        cmax_prev=cmax;
        
       
     
    else %for other images

        diffrnmax=r*c; % maximum difference possible
        [jprev,iprev]=ind2sub([r,c],C_prev.PixelIdxList{1,cid_prev});
        
        C2=zeros(size(B_our));
        for iii=1:length(jprev)
            C2(jprev(iii),iprev(iii))=1;
        end
        
        divisor=2; % sometimes, curve was not found for minarea, so minarea

        loopcount = 1; % how many iteration should be considered is initialized
        while cid==0
            for ci=1:cL
                cL2=length(C1.PixelIdxList{1,ci}); % as loop is continuing,

                C3=zeros(size(B_our));

                [jnow,inow]=ind2sub([r,c],C1.PixelIdxList{1,ci});

                for iii=1:length(jnow) %this image curve is taken and 
                    C3(jnow(iii),inow(iii))=1;
                end            

                        
                    diffrn=sqrt((mean(jnow)-mean(jprev)).^2+ (mean(inow)-mean(iprev)).^2);
%                 diffrn=sum(sum(abs(C2-C3)));%first sum is for column total,

                if diffrn<diffrnmax  && cL2>minarea/divisor % here divisor 
                    %is used for easier findings of curve 
    
                   diffrnmax=diffrn;
                   cid=ci;
                   cmax=length(C1.PixelIdxList{1,ci});
                end
               
            end

            divisor=divisor*2; %%%%%%%%%%%%%
            
            loopcount = loopcount+1;
            if loopcount > 1000 %if a curve is not found, the iteration can
                % run maximum 1000 times
                break
            end
        end

        divisor;
        C_prev=C1;
        cid_prev=cid;
        cmax_prev=cmax;    

    end


    C2=zeros(size(B_our));
    [jnow,inow]=ind2sub([r,c],C1.PixelIdxList{1,cid});
    for iii=1:length(jnow)
        C2(jnow(iii),inow(iii))=1;
    end
    count2g=sum(sum((I1(:, :, 2))))/(r*c)
     C3=C2;
    D2=imcomplement(C2);
       
     bound=10;
    for iii=1:width
        for jjj=1:height
            if iii<=bound || iii>=width-bound || jjj<=bound || jjj>=height-bound
                D2(jjj,iii)=0;
            end
        end
    end
      
    E1=bwconncomp(D2);
    eL=length(E1.PixelIdxList);
    if E1.NumObjects==0
        continue
     end
    eid=0;

    
    if ii==tmfr(1)
        varmax=height/2;
        for ei=1:eL
            eL2=length(E1.PixelIdxList{1,ei});

            [jlist,ilist]=ind2sub([r,c],E1.PixelIdxList{1,ei});
            imean=round(mean(ilist));
            jmean=round(mean(jlist));
            varpos=1./abs(imean-height/2);
%             varpos=sum(((jlist-jmean).^2+(ilist-imean).^2).^0.5)/length(jlist);
            % comparison parameter
            if varpos<varmax %area comparison is not needed
               varmax=varpos;
               eid=ei;

            end

        end
        E_prev=E1;
        eid_prev=eid;
    else
        diffrnmax=r*c;
        [jprev,iprev]=ind2sub([r,c],E_prev.PixelIdxList{1,eid_prev});
        E2=zeros(size(B_our));
        for iii=1:length(jprev)
            E2(jprev(iii),iprev(iii))=1;
        end
        divisor=2;
        
        loopcount = 1;
        while eid==0
            
            for ei=1:eL
                eL2=length(E1.PixelIdxList{1,ei});

                E3=zeros(size(B_our));

                [jnow,inow]=ind2sub([r,c],E1.PixelIdxList{1,ei});
                
                for iii=1:length(jnow)
                    E3(jnow(iii),inow(iii))=1;
                end            
%                 imean=mean(inow);
                diffrn=sqrt((mean(jnow)-mean(jprev)).^2+ (mean(inow)-mean(iprev)).^2);
%                 diffrn=sum(sum(abs(E2-E3)));

                if diffrn<diffrnmax && eL2>minarea/divisor

                   diffrnmax=diffrn;
                   eid=ei;

                end
            end

            divisor=divisor*2; %%%%%%%%%%%
            
            loopcount = loopcount+1;
            if loopcount > 1000
                break
            end
        end
        
%         if loopcount > 1000
%             continue
%         end
%         
        divisor;
        E_prev=E1;
        eid_prev=eid;
   end
    
        C4=zeros(size(C2));
    for ei=1:eL
        if ei==eid
            [ejlist,eilist]=ind2sub([r,c],E1.PixelIdxList{1,ei}); 
            for jjj=1:length(ejlist)
                C4(ejlist(jjj),eilist(jjj))=1;
            end
        end
    end   
    
    figure(104)
    imshow(C4)  

    countblack=0;

    for j=1:r % r is no of rows 
        for i=1:c  % c is no of columns 
            if C4(j,i)==0
               countblack=countblack+1; 
            end
        end
    end
    count2=countblack;
    if ii==tmfrinit
       count1=countblack;
count1g=sum(sum((I1(:, :, 2))))/(r*c)
 
    else
        sg(ii-tmfrinit+1)=(count2-count1);
        sgg(ii-tmfrinit+1)=(count2g-count1g); %for intensity (number of values of pixels(greyscale))

    end        

end

% for k=1:Fs:tmfrmax
%      if k>tmfrmax-Fs
%         break 
%     end
% sgn(k:k+Fs)=detrend(movmean((sg(k:k+Fs))/max(abs(sg)), 5));
% sggn(k:k+Fs)=detrend(movmean(sgg(k:k+Fs)/max(abs(sgg)), 5)); %Intensity Curve
% end

sgn=movmean((sg)/max(abs(sg)), 5);
sggn=movmean(sgg/max(abs(sgg)), 5); %Intensity Curve

grid on
figure
plot(time,(sgn))
xlabel('time (sec)','FontName','Times New Roman','FontSize', 20)
ylabel('normialized signal','FontName','Times New Roman','FontSize', 20)
legend('For Fingerprint')

% figure
% subplot(2,1,2)
hold on
plot(time,(-sggn))
legend('From Fingertip Movement','From Average Intensity')
% 
% p1=90
% p2=100
% 
% 
% [pks,locs1] = findpeaks(sgn(p1*30:p2*30),'MinPeakDistance',20);
% RRI1=(diff(locs1))/30;
% % M1=1./median(RRI1)*60
% M1=1./median(RRI1)*60;
% sprintf('From FIngertip image, Heartrate= %2d ', M1)
% 
% [pks,locs2] = findpeaks(sggn(p1*30:p2*30),'MinPeakDistance',20);
% RRI2=(diff(locs2))/30;
% M2=1./median(RRI2)*60;
% sprintf('From Avarage intensity, Heartrate= %2d ', M2)

% save ('signals.mat', 'sgn' , 'sggn');