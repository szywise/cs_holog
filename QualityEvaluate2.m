close all; clear all; clc;

for i = 1:2:10
    
    img1 = imread(  sprintf('D:/BaiduYunDownload/SPmag-CS-Code_new/img_580_p.png')  );
    img1 = rgb2gray(img1);
    
    img1r = imread(  sprintf('D:/BaiduYunDownload/SPmag-CS-Code_new/img_%d.bmp',i)  );
    
    [m, n] = size(img1);
    
    disp( sprintf('## %d ##', i) )
    disp('## psnr ##')
    psnr_r = cpsnr(img1, img1r)
    
end