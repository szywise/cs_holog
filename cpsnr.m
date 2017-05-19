function [result] = cpsnr(img1, img2)
%figure;imshow(img1);
%figure;imshow(img2);
img1 = double(img1);
img2 = double(img2);

[m, n] = size(img1);

%img1 = img1 - mean(img1(:));
%img2 = img2 - mean(img2(:));

d2sum = sum ( sum ( (img1 - img2).^2 ) );

result = d2sum/m/n;

result = 10 * log10(255^2/result);