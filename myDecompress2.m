function [result] = myDecompress2(sample)

a = 1024; b = 1024;
sample = double(sample);
mult = double(sample(1, 1)/10);
len = size(sample, 1);
meanS = mean(sample);

M_measure = zeros(a, b);
mask_matrix = zeros(a, b);
mask_matrix( floor(a/2-a/mult+0.01):floor(a/2+a/mult+0.01), floor(b/2-b/mult+0.01):floor(b/2+b/mult+0.01)) = 1;

id = 2;
for i = 1:b
    for j = 1:a
        if(mask_matrix(j, i) > 0)
            M_measure(j, i) = sample(id, 1);
            id = id+1;
        else
            M_measure(j, i) = meanS;
        end
    end
end

MI = mean(M_measure(:));
Scale = max(max(M_measure));
M_measure = M_measure/Scale;

shift_v = -0.45;
mask_matrix=round(rand(a,b)-shift_v);
mask_matrix( floor(a/2-a/mult+0.01):floor(a/2+a/mult+0.01), floor(b/2-b/mult+0.01):floor(b/2+b/mult+0.01)) = 1;
M_measure = mask_matrix.*M_measure; 

%mask_matrix(mask_matrix~=1) = 1;
%  ͼ�������ֵ
M_0 = zeros(a,b);

%  2�������ܱ���Ȩ����ֵ
lambda=0.1;
%  �ݶ�����
grad=grad_2norm(mask_matrix,M_0,M_measure)+lambda*grad_1norm_tv(M_0);
gradx=grad;
dire=-grad;

%  ��������Line Search������
alpha=0.01;
beta=0.6;
tau0=1;
index_1=0;

%  �ݶ�������׼
epsi=1e-4;
%  ����������
max_iter=50;
%  ��ǰ��������
k=0;
%  ���ջָ���ͼ��
M_recover=M_0;

%  ��������
while(norm(grad,'fro')>epsi && k<max_iter)
    
    %  ��ֵ
    tau=tau0;
    num=0;
    
    %  ��������Line Search��
    while ((f_2norm(mask_matrix,M_recover+tau*dire,M_measure)+...
            lambda*f_1norm_tv(M_recover+tau*dire))>...
           (f_2norm(mask_matrix,M_recover,M_measure)+...
            lambda*f_1norm_tv(M_recover)+alpha*tau*real(conj(grad).*dire)))
        tau=beta*tau;
        num=num+1;
    end
    
    %  ����Ӧ��ֵ
	if num>2
		tau0 = tau0*beta;
	end 
	if num<1
		tau0 = tau0/beta;
    end
    
    %  �ָ�ͼ������
    M_recover=M_recover+tau*dire;
    grad_0=grad;  
    
    %  �ݶ���ʾ
    grad_show=norm(grad,'fro');
    %disp('�ݶ���')
    %disp(grad_show)
    
    %  �ݶ�����
    grad=grad_2norm(mask_matrix,M_recover,M_measure)+lambda*grad_1norm_tv(M_recover);
    gamma=norm(grad,'fro')^2/norm(grad_0,'fro')^2;
    dire=-grad+gamma*dire;  
    
    %  ������������
    k=k+1;
    %disp('����������')
    %disp(k)
    
    %  2�������ܱ���Ȩ����ֵ�������������������lambda�ɹ̶���
    lambda=lambda*0.90;

end

%  �����ʾ
result = uint8(Scale*M_recover);
%result = M_measure;

%  2����
function TT=f_2norm(mask_matrix,T,S)  
TT=norm(mask_matrix.*T-S,'fro')^2;

%  �ܱ��
function TT=f_1norm_tv(Solution)
Solution=[Solution(:,1) Solution Solution(:,end)];
Solution=[Solution(1,:);Solution;Solution(end,:)];
df_x=(Solution(2:end-1,3:end)-Solution(2:end-1,1:end-2))/2;
df_y=(Solution(3:end,2:end-1)-Solution(1:end-2,2:end-1))/2;
TT=sum(sum(sqrt(df_x.^2+df_y.^2)));

%  2�������ݶ�
function TT=grad_2norm(mask_matrix,T,S)   
TT=2*mask_matrix.*(mask_matrix.*T-S);

%  �ܱ����ݶ�
function TT=grad_1norm_tv(Solution)

epsx=1e-14;  %  ��ֹ�ݶ����޴�

Solution=[Solution(:,1) Solution Solution(:,end)];
Solution=[Solution(1,:);Solution;Solution(end,:)];

%  ������
xx_1=Solution(2:end-1,2:end-1)-Solution(2:end-1,3:end);
yy_1=Solution(2:end-1,2:end-1)-Solution(3:end,2:end-1);
%  ��ߵ���
xx_2=Solution(2:end-1,1:end-2)-Solution(2:end-1,2:end-1);
yy_2=Solution(2:end-1,1:end-2)-Solution(3:end,1:end-2);
%  �ϱߵ���
xx_3=Solution(1:end-2,2:end-1)-Solution(1:end-2,3:end);
yy_3=Solution(1:end-2,2:end-1)-Solution(2:end-1,2:end-1);

%  �ݶ�
grad_1=sqrt(xx_1.^2+yy_1.^2+epsx);
grad_2=sqrt(xx_2.^2+yy_2.^2+epsx);
grad_3=sqrt(xx_3.^2+yy_3.^2+epsx);
            
%  �ܱ����ݶ�
TT=(xx_1./grad_1+yy_1./grad_1-xx_2./grad_2-yy_3./grad_3);



