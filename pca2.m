function [y,D] = pca2(sig)%mixedsig的行表示信号，列表示采样点

%程序说明：y = pca(mixedsig)，程序中mixedsig为 n*T 阶混合数据矩阵，n为信号个数，T为采样点数
% y为 m*T 阶主分量矩阵。
% n是维数，T是样本数。

if nargin == 0
    error('You must supply the mixed data as input argument.');
end
if length(size(sig))>2
    error('Input data can not have more than two dimensions. ');
end
if any(any(isnan(sig)))
    error('Input data contains NaN''s.');
end
%——————————————去均值———————————— %，消除量纲、量级不统一而造成结果不精确
[m,n]=size(sig);
sigmean=mean(sig');

for i =  1:m           %数据矩阵的标准化 
    for j = 1:n
        sig(i,j) = (sig(i,j)-sigmean(i));
    end
end

% [Dim,NumofSampl] = size(mixedsig);%对初始化后的矩阵再次求其维数
% oldDimension = Dim;
%把新矩阵的行赋给oldDimension
% fprintf('Number of signals: %d\n',m);
% fprintf('Number of samples: %d\n',n);
% fprintf('Calculate PCA...');
%covarianceMatrix = cov(mixedsig');    %计算协方差矩阵
covarianceMatrix =(sig*transpose(sig))./(n-1); %计算协方差矩阵

[E,D] = eig(covarianceMatrix); 
[l,q]=size(E);
%计算协方差矩阵的特征值和特征向量
%——————————降序排列特征值与特征向量——————————
neweigenvalues=diag(D);
eigenvalues = flipud(neweigenvalues);%对特征值进行降序排序,即从大到小
% eigenvalues = flipud(eigenvalues);%对特征值升序排列，即从小到大
Enew=fliplr(E);
c=sum(neweigenvalues);
a=neweigenvalues./c;
p0=1;
h=0;
for i=1:l
     h=h+eigenvalues(i)/c;
     if(h<0.8)  
      p0=p0+1;
     end
end 
%}
% save w
%fprintf('p is: %d',p);
%diag(D)是把特征值矩阵写成特征值向量的形式
 %sort(diag(D))表示把一个向量按从小到大的顺序排序
 %-------------------主分量的计算及重构-------------------%
 y1=Enew(1:l,1:p0).'*sig;
% y1=Enew(1:l,1:p0).'*sig;
% y=Enew(1:l,1:p0)*y1+sum(meanValue)/m*ones(m,n);
y=Enew(1:l,1:p0)*y1+sigmean'*ones(1,n);
%  y=Enew(1:l,1:p0)*y1+ones(m,1)*sigmean;

