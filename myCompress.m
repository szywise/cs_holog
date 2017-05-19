function [result, file_t, ratio_t] = myCompress(img, file_o, ratio)
M_image = double(img);
[a, b] = size(M_image);

mult = 2; 
ratio_t = 1;
while(ratio_t > ratio)
    mult = mult+0.1;
    mask_matrix = zeros(a, b);
    mask_matrix( floor(a/2-a/mult+0.01):floor(a/2+a/mult+0.01), floor(b/2-b/mult+0.01):floor(b/2+b/mult+0.01)) = 1;
    M_measure = mask_matrix.*M_image; 
    result = uint8( [mult*10; M_measure(mask_matrix~=0)] );
    save ./tmp/tmp.mat result mult
    file_t = dir('./tmp/tmp.mat');
    ratio_t = double(file_t.bytes) / double(file_o.bytes);
end
ss = sum(sum(mask_matrix));
disp( sprintf('sample num = %d, sample ratio = %2.2f %%', ss, ss/a/b*100) );