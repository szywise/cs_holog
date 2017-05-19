function [holog] = myRecord(img1, img2, z1, z2, color)

ima = img2;
%% 图像扩展 
[M,N,p] = size(ima);
c = color;
ll = 10;

if(1)
    m = 2 * M;
    n = 2 * N;
    posx = (m-M) / 2;
    posy = (n-N) / 2;

    A = zeros(m,n);
    A(posx:posx+M-1, posy:posy+N-1) = double(ima(:, :, c))/255;
else
    A = double(ima(:, :, c))/255;
    m = M;
    n = N;
end
%figure;imshow(A, []);

%for i = 300:20:300+20
for i = 256:30:256+30
    
    ima = zeros(1024, 1024);
    for j = 256:768
        ima(i:i+29, j) = A(i:i+29, j);
    end
    
    %figure;imshow(uint8(ima));
    
    %% 全息图制作
    lamd=0.532e-003;
    z0 = z1 + (z2-z1)/20;

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

    %r = 0.1.*rand(m,n);                             %随机位相,为物光波加上一个随机相位，平滑fft变换之后的频率普
    %u1 = exp(1i*2*pi*r);                       %物体自身的随机相位因子
    u2 = exp(1i*pi*(X.^2+Y.^2)/(lamd*z0));     %菲涅尔全息的第二个因子，表示物体发出的光波经过一段距离后的分布情况
    %u = A.*u1.*u2;
    u = ima.*u2;

    w = fftshift(fft2( (u) ));
    phase = exp(1i*2*pi*z0/lamd)/(1i*lamd*z0)*exp(1i*pi/lamd/z0*(X1.^2+Y1.^2));    %菲涅尔衍射积分号前的相位因子
    w = w.*phase;
    T = L0/n;               %空域取样间隔
    w2 = w*T*T;              %二维离散变换量值补偿
    wt = w2/7.5;
end
cnt = 1;
%for i = 400:20:768
for i = 286:30:556
    
    ima = zeros(1024, 1024);
    for j = 256:768
        ima(i:i+29, j) = A(i:i+29, j);
    end
    
    %figure;imshow(uint8(ima));
    
    %% 全息图制作
    lamd=0.532e-003;
    z0 = z1 + (z2-z1)/9*cnt;
    cnt = cnt+1;

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

    %r = 0.1.*rand(m,n);                             %随机位相,为物光波加上一个随机相位，平滑fft变换之后的频率普
    %u1 = exp(1i*2*pi*r);                       %物体自身的随机相位因子
    u2 = exp(1i*pi*(X.^2+Y.^2)/(lamd*z0));     %菲涅尔全息的第二个因子，表示物体发出的光波经过一段距离后的分布情况
    %u = A.*u1.*u2;
    u = ima.*u2;

    w = fftshift(fft2( (u) ));
    phase = exp(1i*2*pi*z0/lamd)/(1i*lamd*z0)*exp(1i*pi/lamd/z0*(X1.^2+Y1.^2));    %菲涅尔衍射积分号前的相位因子
    w = w.*phase;
    T = L0/n;               %空域取样间隔
    w2 = w*T*T;              %二维离散变换量值补偿
    wt = (wt + w2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ima = imread('../data/3.png');
[M,N,p] = size(ima);
c = color;
ll = 10;

if(1)
    m = 2 * M;
    n = 2 * N;
    posx = (m-M) / 2;
    posy = (n-N) / 2;

    A = zeros(m,n);
    A(posx:posx+M-1, posy:posy+N-1) = double(ima(:, :, c))/255;
else
    A = double(ima(:, :, c))/255;
    m = M;
    n = N;
end
%figure;imshow(A, []);

for i = 300:20:300+20
    
    ima = zeros(1024, 1024);
    for j = 256:768
        ima(i:i+19, j) = A(i:i+19, j);
    end
    
    %figure;imshow(uint8(ima));
    
    %% 全息图制作
    lamd=0.532e-003;
    z0 = z1 + (z2-z1)/20;

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

    %r = 0.1.*rand(m,n);                             %随机位相,为物光波加上一个随机相位，平滑fft变换之后的频率普
    %u1 = exp(1i*2*pi*r);                       %物体自身的随机相位因子
    u2 = exp(1i*pi*(X.^2+Y.^2)/(lamd*z0));     %菲涅尔全息的第二个因子，表示物体发出的光波经过一段距离后的分布情况
    %u = A.*u1.*u2;
    u = ima.*u2;

    w = fftshift(fft2( (u) ));
    phase = exp(1i*2*pi*z0/lamd)/(1i*lamd*z0)*exp(1i*pi/lamd/z0*(X1.^2+Y1.^2));    %菲涅尔衍射积分号前的相位因子
    w = w.*phase;
    T = L0/n;               %空域取样间隔
    w2 = w*T*T;              %二维离散变换量值补偿
    wt2 = w2;
end
cnt = 1;
for i = 400:20:768
    
    ima = zeros(1024, 1024);
    for j = 256:768
        ima(i:i+19, j) = A(i:i+19, j);
    end
    
    %figure;imshow(uint8(ima));
    
    %% 全息图制作
    lamd=0.532e-003;
    z0 = z2 - (z2-z1)/20*cnt;
    cnt = cnt+1;

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

    %r = 0.1.*rand(m,n);                             %随机位相,为物光波加上一个随机相位，平滑fft变换之后的频率普
    %u1 = exp(1i*2*pi*r);                       %物体自身的随机相位因子
    u2 = exp(1i*pi*(X.^2+Y.^2)/(lamd*z0));     %菲涅尔全息的第二个因子，表示物体发出的光波经过一段距离后的分布情况
    %u = A.*u1.*u2;
    u = ima.*u2;

    w = fftshift(fft2( (u) ));
    phase = exp(1i*2*pi*z0/lamd)/(1i*lamd*z0)*exp(1i*pi/lamd/z0*(X1.^2+Y1.^2));    %菲涅尔衍射积分号前的相位因子
    w = w.*phase;
    T = L0/n;               %空域取样间隔
    w2 = w*T*T;              %二维离散变换量值补偿
    wt2 = (wt2 + w2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ima = img1;
%% 图像扩展 
[M,N,p] = size(ima);
if (1)
    m = 2 * M;
    n = 2 * N;
    posx = (m-M) / 2;
    posy = (n-N) / 2;

    A = zeros(m,n);
    A(posx:posx+M-1, posy:posy+N-1) = double(ima(:, :, c))/255;
else
    A = double(ima(:, :, c))/255;
    m = M;
    n = N;
end
%figure;imshow(A, []);

%% 全息图制作
%lamd=0.532e-003;
z0=z2;

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
  
%r = 0.1.*rand(m,n);                             %随机位相,为物光波加上一个随机相位，平滑fft变换之后的频率普
%u1 = exp(1i*2*pi*r);                       %物体自身的随机相位因子
u2 = exp(1i*pi*(X.^2+Y.^2)/(lamd*z0));     %菲涅尔全息的第二个因子，表示物体发出的光波经过一段距离后的分布情况
%u = A.*u1.*u2;
u = A.*u2;

w = fftshift(fft2( (u) ));
phase = exp(1i*2*pi*z0/lamd)/(1i*lamd*z0)*exp(1i*pi/lamd/z0*(X1.^2+Y1.^2));    %菲涅尔衍射积分号前的相位因子
w = w.*phase;
T = L0/n;               %空域取样间隔
w = w*T*T;              %二维离散变换量值补偿
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pp = 1.95;
pp = .75;
pr=2*pi*X1*sin(pp*pi/180)/lamd+2*pi*Y1*sin(pp*pi/180)/lamd;
refer = exp(-1i*pr);%再现

%xr = (L0 * 1/2);                                                         %第一象限的中点作为再现像原像的中点  
%refer = exp(1i*pi*((X1+xr).^2+(Y1+xr).^2)/(lamd*z0));                 %参考光，参考光的X1坐标加上了一个量Xr,球面波中心坐标为（X1+xr，Y1+xr）
%refer = exp(1i*pi*(X1+Y1)) / (lamd*z0);

%w1 = w2 + 5*refer;
w1 = w + wt + wt2 + 10*refer;
w1 = abs(w1).^2;                                                   %求出全息平面上的强度分布
holog = mat2gray(w1);                                                 %进行归一化处理，得到我们想要的干涉条纹强度分布


%% 全息图储存
%figure;
%imshow(w1,[]);
%colormap(gray)
%imwrite(w1,outname)                                          %将所得到的干涉图（全息图保存下来）

