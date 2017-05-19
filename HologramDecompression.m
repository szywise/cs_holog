clear all; close all; clc;

%%
load ../data/compressData.mat

%% decompress hologram Red
hologR = myDecompress(sampleR);
imwrite(hologR, '../data/hologram_R_decompress.bmp')

%% decompress hologram Green
hologG = myDecompress(sampleG);
imwrite(hologG, '../data/hologram_G_decompress.bmp')

%% decompress hologram Blue
hologB = myDecompress(sampleB);
imwrite(hologB, '../data/hologram_B_decompress.bmp')