function [y,D] = pca2(sig)%mixedsig���б�ʾ�źţ��б�ʾ������

%����˵����y = pca(mixedsig)��������mixedsigΪ n*T �׻�����ݾ���nΪ�źŸ�����TΪ��������
% yΪ m*T ������������
% n��ά����T����������

if nargin == 0
    error('You must supply the mixed data as input argument.');
end
if length(size(sig))>2
    error('Input data can not have more than two dimensions. ');
end
if any(any(isnan(sig)))
    error('Input data contains NaN''s.');
end
%����������������������������ȥ��ֵ������������������������ %���������١�������ͳһ����ɽ������ȷ
[m,n]=size(sig);
sigmean=mean(sig');

for i =  1:m           %���ݾ���ı�׼�� 
    for j = 1:n
        sig(i,j) = (sig(i,j)-sigmean(i));
    end
end

% [Dim,NumofSampl] = size(mixedsig);%�Գ�ʼ����ľ����ٴ�����ά��
% oldDimension = Dim;
%���¾�����и���oldDimension
% fprintf('Number of signals: %d\n',m);
% fprintf('Number of samples: %d\n',n);
% fprintf('Calculate PCA...');
%covarianceMatrix = cov(mixedsig');    %����Э�������
covarianceMatrix =(sig*transpose(sig))./(n-1); %����Э�������

[E,D] = eig(covarianceMatrix); 
[l,q]=size(E);
%����Э������������ֵ����������
%��������������������������������ֵ������������������������������
neweigenvalues=diag(D);
eigenvalues = flipud(neweigenvalues);%������ֵ���н�������,���Ӵ�С
% eigenvalues = flipud(eigenvalues);%������ֵ�������У�����С����
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
%diag(D)�ǰ�����ֵ����д������ֵ��������ʽ
 %sort(diag(D))��ʾ��һ����������С�����˳������
 %-------------------�������ļ��㼰�ع�-------------------%
 y1=Enew(1:l,1:p0).'*sig;
% y1=Enew(1:l,1:p0).'*sig;
% y=Enew(1:l,1:p0)*y1+sum(meanValue)/m*ones(m,n);
y=Enew(1:l,1:p0)*y1+sigmean'*ones(1,n);
%  y=Enew(1:l,1:p0)*y1+ones(m,1)*sigmean;

