%%------------------------------------------------------------------------%%
% This code is to detect the cirrus in remote sensing images by merging 
% saliency features in different scales based on PCA. 
% This is section 4.2 of Yuhan Liu's thesis. 
% Coder: Yuhan Liu
% Date: 2020.11.12
%%------------------------------------------------------------------------%%
clear all;
close all;
clc;
tic
I=imread('data278.tif');% input the data
GT=logical(imread('GT278.jpg'));
I1=double(I).^0.67;% gamma transform to adjust the image
% figure,imshow(I1,[],'border','tight')
%% acquire ground truth,when run the algorithm,comment this part
% figure,imshow(I1,[])
% h = drawfreehand;
% Mask1 = createMask(h);
% figure,imshow(I1,[])
% h = drawfreehand;
% Mask2 = createMask(h);
% figure,imshow(I1,[])
% h = drawfreehand;
% Mask3 = createMask(h);
% Mask=Mask1+Mask2+Mask3;
% GT = activecontour(I1,Mask ,1000);
% figure,imshow(GT)
% imwrite(GT,'GT278.jpg');
% figure,imshow(I1,[],'border','tight')
%% feature extraction in different scales
sigma=[0.5,1.5,2.5,3.5];
for i=1:length(sigma)
G=fspecial('gaussian',[7 7],sigma(i));% generate gaussian filters
IF=imfilter(I1,G,'replicate');
FT{i}=FrequencyTuned(IF);
% entr(i)=entropy(FT{i});
end
%% feature fusion based on PCA
[M N]=size(FT{1});
S=zeros(length(sigma),M*N);
S(1,:)=FT{1}(1:M*N);
S(2,:)=FT{2}(1:M*N);
S(3,:)=FT{3}(1:M*N);
S(4,:)=FT{4}(1:M*N);
[Y d]=pca2(S);
% P1=reshape(Y(1,:),M,N);
% P2=reshape(Y(2,:),M,N);
% P3=reshape(Y(3,:),M,N);
% P4=reshape(Y(4,:),M,N);
LF=FT{1}*(sum(d(1,:))/sum(d(:)))+FT{2}*(sum(d(2,:))/sum(d(:)))+FT{3}*(sum(d(3,:))/sum(d(:)))+FT{4}*(sum(d(4,:))/sum(d(:)));% low frequency fusion

for i=2:length(sigma)
    HFT{i-1}=FT{i-1}-FT{i};
end
S=zeros(length(sigma)-1,M*N);
S(1,:)=FT{1}(1:M*N);
S(2,:)=FT{2}(1:M*N);
S(3,:)=FT{3}(1:M*N);
% S(4,:)=FT{4}(1:M*N);
[YH dh]=pca2(S);
HF=HFT{1}*(sum(dh(1,:))/sum(dh(:)))+HFT{2}*(sum(dh(2,:))/sum(dh(:)))+HFT{3}*(sum(dh(3,:))/sum(dh(:)));% high frequency fusion
FF=HF+LF;% final fusion
% for i=1:length(sigma)
%     Diff{i}=FT{i}-FF;
% end

FF=(FF-min(FF(:)))./(max(FF(:))-min(FF(:)));
% figure,imshow(FF,[],'border','tight')
% figure,imshow(activecontour(I,im2bw(FF)),'border','tight')
 detect=activecontour(I,im2bw(FF));
 toc
 ju %for experiment stopping
%% cloud detection
T=0:0.01:0.99;
for i=1:length(T)
  
%         detect=activecontour(I,im2bw(FF,T(i)));
    
% figure,imshow(detect,'border','tight')
detect=im2bw(FF,T(i));
f2=detect;
[m n]=size(GT);
tp=sum(sum(GT.*f2));
fn=sum(sum(GT))-tp;
fp=sum(sum(f2))-tp;
tn=m*n-tp-fn-fp;
acc=(tp+tn)/(m*n);
recall=tp/(tp+fn);
pre=tp/(tp+fp);
f=((0.3+1)*pre*recall)/(0.3*pre+recall);
fpr=fp/(tn+fp);
po=acc;
pe=((tp+fn)*(tp+fp)+(fp+tn)*(fn+tn))/((m*n).^2);
k=(po-pe)/(1-pe);
FPR(i)=fpr;
TPR(i)=recall;
PRE(i)=pre;
REC(i)=recall;
K(i)=k;
F(i)=f;
end

%    F=(F-min(F))./(max(F)-min(F));

% % % 
figure(1),hold on, plot(FPR,TPR)
figure(2),hold on, plot(TPR,PRE)
%     FPR=(FPR-min(FPR))./(max(FPR)-min(FPR));
%   TPR=(TPR-min(TPR))./(max(TPR)-min(TPR));
%    PRE=(PRE-min(PRE))./(max(PRE)-min(PRE));

 AUC=trapz(FPR,TPR);
AUCPR=trapz(TPR,PRE);
mean(K)
mean(F)