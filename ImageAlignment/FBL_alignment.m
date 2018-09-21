%This code is adopted from:
%https://ww2.mathworks.cn/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation

f=im2double(imread('template_image.tif'));
g=im2double(imread('align_image.tif'));
usfac = 100;
merr=zeros(1,360);
si=size(f);
mGreg=zeros(si(1),si(2),360);
%Rotate the aligning image 360 degree by one degree a step
for i=0:359
    A=imrotate(g,i,'bilinear','crop');
    [output, Greg] = dftregistration(fft2(f),fft2(A),usfac);
    merr(i+1)=output(1);
    mGreg(:,:,i+1)=abs(ifft2(Greg));
end
[mi,idx]=min(merr);
total=(f+mGreg(:,:,idx))/2;
imwrite(mGreg(:,:,idx),'registerision_image.tif');
imwrite(total,'new_template.tif');
%{
for n=3:length(file)
    f=total;
    g=im2double(imread(strcat('try_result/',file(n).name)));
    usfac = 100;
    merr=zeros(1,360);
    si=size(f);
    mGreg=zeros(si(1),si(2),360);
    for i=0:359
        A=imrotate(g,i,'bilinear','crop');
        [output, Greg] = dftregistration(fft2(f),fft2(A),usfac);
        merr(i+1)=output(1);
        mGreg(:,:,i+1)=abs(ifft2(Greg));
    end
    [mi,idx]=min(merr);
    total=(total*(n-1)+mGreg(:,:,idx))/n;
    imwrite(mGreg(:,:,idx),strcat('FBL_try_register/reg_',num2str(n-1),'.tif'));
    imwrite(total,strcat('FBL_try_ref/ref_',num2str(n-1),'.tif'));
end
%}    
%total=(f*17+mGreg(:,:,idx))/18;
%imshow(total);
%imwrite(mGreg(:,:,idx),'FBL_register/reg_17.tif');
%imwrite(total,'FBL_ref/ref_17.tif');

%[ma,idx2]=max(merr);
%tr2=(f+mGreg(:,:,idx2))/2;

%% 
% The pixel shift error (difference between the true and obtained shifts)
% is 0.0016 and 0.0043 in the x and y directions respectively. Well within
% the expected accuracy of 0.01. Notice that using the conventional zero-padded 
% FFT approach with the same accuracy, would
% require computation of a 25,600x25,600 FFT, which would require more than
% 19 Gbytes of RAM and a very comfortable chair.
%
% The following plot shows the reference image and the registered image.
