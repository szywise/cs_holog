function [result] = cq(img1, img2)
img1 = double(img1);
img2 = double(img2);

[m, n] = size(img1);

sum1 = sum ( sum( img1 ) );
sum2 = sum ( sum( img2 ) );

avg1 = sum1 / m / n;
avg2 = sum2 / m / n;

cov1 = sqrt ( sum ( sum ( ( (img1 - avg1).^2 ) ) ) /m/n);
cov2 = sqrt ( sum ( sum ( ( (img2 - avg2).^2 ) ) ) /m/n);

covxy = sqrt ( sum ( sum ( ( (img1 - avg1).*(img2 - avg2) ) ) ) /m/n);

result = 4 * avg1 * avg2 * covxy^2 / (cov1^2 + cov2^2) / (avg1^2 + avg2^2);