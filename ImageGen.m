clc; close all; clear;

%% read images
img1 = imread('img1.bmp');
img2 = imread('img2.bmp');

%% padding
image1 = zeros(512, 512);
imgae2 = zeros(512, 512);