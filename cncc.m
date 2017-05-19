function [result] = cncc(img1, img2)
img1 = double(img1);
img2 = double(img2);

mul = img1.*img2;
squ = img1.^2;

sum_mul = sum( sum (mul) );
sum_squ = sum( sum (squ) );

result = sum_mul/sum_squ;
if(result>1)
    result = 1/result;
end
