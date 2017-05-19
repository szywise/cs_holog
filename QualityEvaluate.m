close all; clear all; clc;
%% z = 500
z = 500;

img1 = imread(  sprintf('../data/img_%d_R.bmp',z)  );
img2 = imread(  sprintf('../data/img_%d_G.bmp',z)  );
img3 = imread(  sprintf('../data/img_%d_B.bmp',z)  );

img1r = imread(  sprintf('../data/img_%d_R_decompress.bmp',z)  );
img2r = imread(  sprintf('../data/img_%d_G_decompress.bmp',z)  );
img3r = imread(  sprintf('../data/img_%d_B_decompress.bmp',z)  );

[m, n] = size(img1);

disp('## z = 500 ##')
disp('## psnr ##')
psnr_r = cpsnr(img1, img1r);
psnr_g = cpsnr(img2, img2r);
psnr_b = cpsnr(img3, img3r);
psnr1 = (psnr_r + psnr_g + psnr_b)/3

disp('## ncc ##')
ncc_r = cncc(img1, img1r);
ncc_g = cncc(img2, img2r);
ncc_b = cncc(img3, img3r);
ncc1 = (ncc_r + ncc_g + ncc_b)/3

disp('## k ##')
k_r = ck(img1, img1r);
k_g = ck(img2, img2r);
k_b = ck(img3, img3r);
k1 = (k_r + k_g + k_b)/3

disp('## q ##')
q_r = cq(img1, img1r);
q_g = cq(img2, img2r);
q_b = cq(img3, img3r);
q1 = (q_r + q_g + q_b)/3

%% z = 530
z = 580;

img1 = imread(  sprintf('../data/img_%d_R.bmp',z)  );
img2 = imread(  sprintf('../data/img_%d_G.bmp',z)  );
img3 = imread(  sprintf('../data/img_%d_B.bmp',z)  );

img1r = imread(  sprintf('../data/img_%d_R_decompress.bmp',z)  );
img2r = imread(  sprintf('../data/img_%d_G_decompress.bmp',z)  );
img3r = imread(  sprintf('../data/img_%d_B_decompress.bmp',z)  );

[m, n] = size(img1);

disp('## z = 530 ##')
disp('## psnr ##')
psnr_r = cpsnr(img1, img1r); 
psnr_g = cpsnr(img2, img2r);
psnr_b = cpsnr(img3, img3r);
psnr2 = (psnr_r + psnr_g + psnr_b)/3

disp('## ncc ##')
ncc_r = cncc(img1, img1r);
ncc_g = cncc(img2, img2r);
ncc_b = cncc(img3, img3r);
ncc2 = (ncc_r + ncc_g + ncc_b)/3

disp('## k ##')
k_r = ck(img1, img1r);
k_g = ck(img2, img2r);
k_b = ck(img3, img3r);
k2 = (k_r + k_g + k_b)/3

disp('## q ##')
q_r = cq(img1, img1r);
q_g = cq(img2, img2r);
q_b = cq(img3, img3r);
q2 = (q_r + q_g + q_b)/3

eva1 = ( ((psnr1-30)/50+0.9) + ncc1 + k1 + q1 )/4;
eva2 = ( ((psnr2-30)/50+0.9) + ncc2 + k2 + q2 )/4;

%% image show
img1  = imread( sprintf('../data/img_%d.bmp',580) );
img1r = imread( sprintf('../data/img_%d_decompress.bmp',580) );

img2  = imread( sprintf('../data/img_%d.bmp',500) );
img2r = imread( sprintf('../data/img_%d_decompress.bmp',500) );

figure;imshow(img1);title('未压缩再现图_z530');
figure;imshow(img1r);title('压缩后再现图_z530');
figure;imshow(img2);title('未压缩再现图_z500');
figure;imshow(img2r);title('压缩后再现图_z500');