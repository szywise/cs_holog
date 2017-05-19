function [holog] = myRecord(img1, img2, z1, z2, color)

ima = img2;
%% ͼ����չ 
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
    
    %% ȫϢͼ����
    lamd=0.532e-003;
    z0 = z1 + (z2-z1)/20;

    %L0 = sqrt(lamd*z0*n);           %FFT����ʱͬʱ�����������λȡ����������ⳡ���
    L0 = ll;
    n1=1:n;
    x = -L0/2 + L0/n*(n1-1);	 					
    y = x;
    [Y,X] = meshgrid(y,x);

    %L = lamd*z0*n/L0;               %FFT�����۲�����������
    L = ll;
    x = -L/2 + L/n*(n1-1);	 					
    y = x;
    [Y1,X1] = meshgrid(y,x); 

    %r = 0.1.*rand(m,n);                             %���λ��,Ϊ��Ⲩ����һ�������λ��ƽ��fft�任֮���Ƶ����
    %u1 = exp(1i*2*pi*r);                       %��������������λ����
    u2 = exp(1i*pi*(X.^2+Y.^2)/(lamd*z0));     %������ȫϢ�ĵڶ������ӣ���ʾ���巢���ĹⲨ����һ�ξ����ķֲ����
    %u = A.*u1.*u2;
    u = ima.*u2;

    w = fftshift(fft2( (u) ));
    phase = exp(1i*2*pi*z0/lamd)/(1i*lamd*z0)*exp(1i*pi/lamd/z0*(X1.^2+Y1.^2));    %������������ֺ�ǰ����λ����
    w = w.*phase;
    T = L0/n;               %����ȡ�����
    w2 = w*T*T;              %��ά��ɢ�任��ֵ����
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
    
    %% ȫϢͼ����
    lamd=0.532e-003;
    z0 = z1 + (z2-z1)/9*cnt;
    cnt = cnt+1;

    %L0 = sqrt(lamd*z0*n);           %FFT����ʱͬʱ�����������λȡ����������ⳡ���
    L0 = ll;
    n1=1:n;
    x = -L0/2 + L0/n*(n1-1);	 					
    y = x;
    [Y,X] = meshgrid(y,x);

    %L = lamd*z0*n/L0;               %FFT�����۲�����������
    L = ll;
    x = -L/2 + L/n*(n1-1);	 					
    y = x;
    [Y1,X1] = meshgrid(y,x); 

    %r = 0.1.*rand(m,n);                             %���λ��,Ϊ��Ⲩ����һ�������λ��ƽ��fft�任֮���Ƶ����
    %u1 = exp(1i*2*pi*r);                       %��������������λ����
    u2 = exp(1i*pi*(X.^2+Y.^2)/(lamd*z0));     %������ȫϢ�ĵڶ������ӣ���ʾ���巢���ĹⲨ����һ�ξ����ķֲ����
    %u = A.*u1.*u2;
    u = ima.*u2;

    w = fftshift(fft2( (u) ));
    phase = exp(1i*2*pi*z0/lamd)/(1i*lamd*z0)*exp(1i*pi/lamd/z0*(X1.^2+Y1.^2));    %������������ֺ�ǰ����λ����
    w = w.*phase;
    T = L0/n;               %����ȡ�����
    w2 = w*T*T;              %��ά��ɢ�任��ֵ����
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
    
    %% ȫϢͼ����
    lamd=0.532e-003;
    z0 = z1 + (z2-z1)/20;

    %L0 = sqrt(lamd*z0*n);           %FFT����ʱͬʱ�����������λȡ����������ⳡ���
    L0 = ll;
    n1=1:n;
    x = -L0/2 + L0/n*(n1-1);	 					
    y = x;
    [Y,X] = meshgrid(y,x);

    %L = lamd*z0*n/L0;               %FFT�����۲�����������
    L = ll;
    x = -L/2 + L/n*(n1-1);	 					
    y = x;
    [Y1,X1] = meshgrid(y,x); 

    %r = 0.1.*rand(m,n);                             %���λ��,Ϊ��Ⲩ����һ�������λ��ƽ��fft�任֮���Ƶ����
    %u1 = exp(1i*2*pi*r);                       %��������������λ����
    u2 = exp(1i*pi*(X.^2+Y.^2)/(lamd*z0));     %������ȫϢ�ĵڶ������ӣ���ʾ���巢���ĹⲨ����һ�ξ����ķֲ����
    %u = A.*u1.*u2;
    u = ima.*u2;

    w = fftshift(fft2( (u) ));
    phase = exp(1i*2*pi*z0/lamd)/(1i*lamd*z0)*exp(1i*pi/lamd/z0*(X1.^2+Y1.^2));    %������������ֺ�ǰ����λ����
    w = w.*phase;
    T = L0/n;               %����ȡ�����
    w2 = w*T*T;              %��ά��ɢ�任��ֵ����
    wt2 = w2;
end
cnt = 1;
for i = 400:20:768
    
    ima = zeros(1024, 1024);
    for j = 256:768
        ima(i:i+19, j) = A(i:i+19, j);
    end
    
    %figure;imshow(uint8(ima));
    
    %% ȫϢͼ����
    lamd=0.532e-003;
    z0 = z2 - (z2-z1)/20*cnt;
    cnt = cnt+1;

    %L0 = sqrt(lamd*z0*n);           %FFT����ʱͬʱ�����������λȡ����������ⳡ���
    L0 = ll;
    n1=1:n;
    x = -L0/2 + L0/n*(n1-1);	 					
    y = x;
    [Y,X] = meshgrid(y,x);

    %L = lamd*z0*n/L0;               %FFT�����۲�����������
    L = ll;
    x = -L/2 + L/n*(n1-1);	 					
    y = x;
    [Y1,X1] = meshgrid(y,x); 

    %r = 0.1.*rand(m,n);                             %���λ��,Ϊ��Ⲩ����һ�������λ��ƽ��fft�任֮���Ƶ����
    %u1 = exp(1i*2*pi*r);                       %��������������λ����
    u2 = exp(1i*pi*(X.^2+Y.^2)/(lamd*z0));     %������ȫϢ�ĵڶ������ӣ���ʾ���巢���ĹⲨ����һ�ξ����ķֲ����
    %u = A.*u1.*u2;
    u = ima.*u2;

    w = fftshift(fft2( (u) ));
    phase = exp(1i*2*pi*z0/lamd)/(1i*lamd*z0)*exp(1i*pi/lamd/z0*(X1.^2+Y1.^2));    %������������ֺ�ǰ����λ����
    w = w.*phase;
    T = L0/n;               %����ȡ�����
    w2 = w*T*T;              %��ά��ɢ�任��ֵ����
    wt2 = (wt2 + w2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ima = img1;
%% ͼ����չ 
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

%% ȫϢͼ����
%lamd=0.532e-003;
z0=z2;

%L0 = sqrt(lamd*z0*n);           %FFT����ʱͬʱ�����������λȡ����������ⳡ���
L0 = ll;
n1=1:n;
x = -L0/2 + L0/n*(n1-1);	 					
y = x;
[Y,X] = meshgrid(y,x);

%L = lamd*z0*n/L0;               %FFT�����۲�����������
L = ll;
x = -L/2 + L/n*(n1-1);	 					
y = x;
[Y1,X1] = meshgrid(y,x); 
  
%r = 0.1.*rand(m,n);                             %���λ��,Ϊ��Ⲩ����һ�������λ��ƽ��fft�任֮���Ƶ����
%u1 = exp(1i*2*pi*r);                       %��������������λ����
u2 = exp(1i*pi*(X.^2+Y.^2)/(lamd*z0));     %������ȫϢ�ĵڶ������ӣ���ʾ���巢���ĹⲨ����һ�ξ����ķֲ����
%u = A.*u1.*u2;
u = A.*u2;

w = fftshift(fft2( (u) ));
phase = exp(1i*2*pi*z0/lamd)/(1i*lamd*z0)*exp(1i*pi/lamd/z0*(X1.^2+Y1.^2));    %������������ֺ�ǰ����λ����
w = w.*phase;
T = L0/n;               %����ȡ�����
w = w*T*T;              %��ά��ɢ�任��ֵ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pp = 1.95;
pp = .75;
pr=2*pi*X1*sin(pp*pi/180)/lamd+2*pi*Y1*sin(pp*pi/180)/lamd;
refer = exp(-1i*pr);%����

%xr = (L0 * 1/2);                                                         %��һ���޵��е���Ϊ������ԭ����е�  
%refer = exp(1i*pi*((X1+xr).^2+(Y1+xr).^2)/(lamd*z0));                 %�ο��⣬�ο����X1���������һ����Xr,���沨��������Ϊ��X1+xr��Y1+xr��
%refer = exp(1i*pi*(X1+Y1)) / (lamd*z0);

%w1 = w2 + 5*refer;
w1 = w + wt + wt2 + 10*refer;
w1 = abs(w1).^2;                                                   %���ȫϢƽ���ϵ�ǿ�ȷֲ�
holog = mat2gray(w1);                                                 %���й�һ�������õ�������Ҫ�ĸ�������ǿ�ȷֲ�


%% ȫϢͼ����
%figure;
%imshow(w1,[]);
%colormap(gray)
%imwrite(w1,outname)                                          %�����õ��ĸ���ͼ��ȫϢͼ����������

