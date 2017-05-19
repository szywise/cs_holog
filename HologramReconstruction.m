%% ÕºœÒ∂¡»°
clear all; close all; clc;

%% Global constraints
z1 = 500;
z2 = 580;%540;

%% reconstruction hologram Red
hologR = imread('../data/hologram_R.bmp');
img1r = myReconstruction(hologR, z1);
img2r = myReconstruction(hologR, z2);
imwrite(img1r, sprintf('../data/img_%d_R.bmp', z1));
imwrite(img2r, sprintf('../data/img_%d_R.bmp', z2));

%% reconstruction hologram Green
hologG = imread('../data/hologram_G.bmp');
img1g = myReconstruction(hologG, z1);
img2g = myReconstruction(hologG, z2);
imwrite(img1g, sprintf('../data/img_%d_G.bmp', z1));
imwrite(img2g, sprintf('../data/img_%d_G.bmp', z2));

%% reconstruction hologram Blue
hologB = imread('../data/hologram_B.bmp');
img1b = myReconstruction(hologB, z1);
img2b = myReconstruction(hologB, z2);
imwrite(img1b, sprintf('../data/img_%d_B.bmp', z1));
imwrite(img2b, sprintf('../data/img_%d_B.bmp', z2));

%% write color image
[a,b] = size(img1r);
img = zeros(a, b, 3);
img(:, :, 1) = img1r; img(:, :, 2) = img1g; img(:, :, 3) = img1b;
imwrite((img),  sprintf('../data/img_%d.bmp',z1)  );
figure;imshow(img);title('Œ¥—πÀı‘Ÿœ÷Õºz500');

img(:, :, 1) = img2r; img(:, :, 2) = img2g; img(:, :, 3) = img2b;
imwrite((img),  sprintf('../data/img_%d.bmp',z2)  );
figure;imshow(img);title('Œ¥—πÀı‘Ÿœ÷Õºz580');