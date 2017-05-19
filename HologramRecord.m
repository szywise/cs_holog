clear all;close all;clc;
%% read images
img1 = imread('../data/1.png');
img2 = imread('../data/2.png');

%% Global constraints
z1 = 500;
z2 = 580;

%% record hologram Red
hologR = myRecord2(img1, img2, z1, z2, 1);
imwrite(hologR, '../data/hologram_R.bmp');

%% record hologram Green
hologG = myRecord2(img1, img2, z1, z2, 2);
imwrite(hologG, '../data/hologram_G.bmp');

%% record hologram Blue
hologB = myRecord2(img1, img2, z1, z2, 3);
imwrite(hologB, '../data/hologram_B.bmp');
