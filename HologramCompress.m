clc; clear all ;close all;

%% read images
hologR = imread('../data/hologram_R.bmp'); fileR = dir('../data/hologram_R.bmp');
hologG = imread('../data/hologram_G.bmp'); fileG = dir('../data/hologram_G.bmp');
hologB = imread('../data/hologram_B.bmp'); fileB = dir('../data/hologram_B.bmp');

%% Global constraints
ratio = 0.05;

%% compress hologram Red
[sampleR, cfileR, ratio_t] = myCompress(hologR, fileR, ratio);
disp( sprintf('compress size / file size = %d / %d bytes, compress ratio = %2.2f %%', cfileR.bytes, fileR.bytes, ratio_t*100));

%% compress hologram Green
[sampleG, cfileG, ratio_t] = myCompress(hologG, fileG, ratio);
disp( sprintf('compress size / file size = %d / %d bytes, compress ratio = %2.2f %%', cfileG.bytes, fileG.bytes, ratio_t*100));

%% compress hologram Blue
[sampleB, cfileB, ratio_t] = myCompress(hologB, fileB, ratio);
disp( sprintf('compress size / file size = %d / %d bytes, compress ratio = %2.2f %%', cfileB.bytes, fileB.bytes, ratio_t*100));

%% write compress data
save ../data/compressData.mat sampleR sampleG sampleB
file_c = dir('../data/compressData.mat');
ratio_c = file_c.bytes / (fileR.bytes+fileG.bytes+fileB.bytes);
disp( sprintf('final compress size / file size = %d / %d bytes, compress ratio = %2.2f %%', file_c.bytes, (fileR.bytes+fileG.bytes+fileB.bytes), ratio_c*100') );