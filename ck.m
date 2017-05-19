function [result] = ck(img1, img2)

[m, n] = size(img1);
img1 = img1(m/2:m, n/2:n);
img2 = img2(m/2:m, n/2:n);

img1 = double(img1);
img2 = double(img2);

w = [-1 , 0 , 1;
      0 , 0.1 , 0;
      1 , 0 ,-1];

y1 = imfilter(img1, w);
y2 = imfilter(img2, w);

sum1 = sum( sum ( abs(y1) ) );
sum2 = sum( sum ( abs(y2) ) );

result = sum2/sum1;