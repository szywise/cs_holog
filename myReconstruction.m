function [img] = myReconstruction(holog, z0)

I1 = holog;
lv = 1.0;
ll = 10;

[m,n,p] = size(I1);
if p == 3
    I1 = rgb2gray(I1);
end

f0 = double(I1);

%% 滤波
f0=f0-mean(mean(f0));                                  %中值滤波
%figure;imshow(f0,[]);

%fm=filter2(fspecial('average',6),f0);                   %局部平均值消除法，3*3格取平均值                  
%f0=f0-fm;

%R=abs(refer).^2;                             %HRO相减法，有问题，结果不对
%O=am.^2;
%f0=f0-R-O;

%gap = 50;
%[m,n] = size(f0);
%B = ones(m, n);
%B(m/2-gap:m/2+gap, n/2-gap:n/2+gap) = 0;
%figure;imshow(B);
%f0 = fftshift(fft2( (f0) ));
%afimg = abs(f0);
%ma = max( max( afimg ) );
%mi = min( min( afimg ) );
%figure;imshow( afimg, [mi ma/1] )
%f0 = f0.*B;
%f0 = (ifft2(fftshift(f0)));
%figure; imshow(f0, []);


%% 再现
[m,n] = size(f0);
A = f0;
lamd=0.532e-003;

%L0 = sqrt(lamd*z0*n);           %FFT计算时同时满足振幅及相位取样条件的物光场宽度
L0 = ll;
n1=1:n;
x = -L0/2 + L0/n*(n1-1);	 					
y = x;
[Y,X] = meshgrid(y,x);

%L = lamd*z0*n/L0;               %FFT计算后观测屏的物理宽度
L = ll;
x = -L/2 + L/n*(n1-1);	 					
y = x;
[Y1,X1] = meshgrid(y,x); 

%pp = 1.95;
pp = 0;
pr=2*pi*X1*sin(pp*pi/180)/lamd+2*pi*Y1*sin(pp*pi/180)/lamd;
refer = exp(-1i*pr);%再现

%xr = (L0 * -1/4);                                                         %第一象限的中点作为再现像原像的中点  
%u1 = exp(1i*pi*((X1+xr).^2+(Y1+xr).^2)/(lamd*z0));                 %参考光，参考光的X1坐标加上了一个量Xr,球面波中心坐标为（X1+xr，Y1+xr）

u2 = exp(1i*pi*(X.^2+Y.^2)/(lamd*z0));     %菲涅尔全息的第二个因子，表示物体发出的光波经过一段距离后的分布情况
u = A.*u2.*refer;
%u = A.*u2;

w = fftshift(fft2( (u) ));
phase = exp(1i*2*pi*z0/lamd)/(1i*lamd*z0)*exp(1i*pi/lamd/z0*(X1.^2+Y1.^2));    %菲涅尔衍射积分号前的相位因子
w = w.*phase;
T = L0/n;               %空域取样间隔
U0 = w*T*T;              %二维离散变换量值补偿

If = U0*conj(U0);
Gmax = max(max(abs(U0)));
Gmin = min(min(abs(U0)));
%figure,imshow(abs(U0),[Gmin Gmax/lv]);
%colormap(gray);
U1=abs(U0);
U1=(U1-Gmin)./(Gmax/lv-Gmin);
%imwrite(abs(U1), outname);
img = abs(U1);