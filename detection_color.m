% ========================================================================
% Detail-Preserving Multi-Exposure Fusion with Edge-Preserving Structural
% Patch Decomposition, IEEE TCSVT,2021
% algorithm Version 1.0
% Copyright(c) 2021, Hui Li, Tsz Nam Chan, Xianbiao Qi and Wuyuan Xie
% All Rights Reserved.
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is hereby
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.


function imgSeqcor= detection_color(imgSeqColor,r)

[h,w,~,n]=size(imgSeqColor);
N = boxfilter(ones(h, w), r);
i_mean=zeros(h,w,n);

C = 0.03 ^ 2 / 2;
structureThres=0.8;
consistencyThres=0.2;
exposureThres=0.02;
ref=ceil(n/2);

%   ref = selectRef(imgSeqColor);
se = strel('disk',2*r+1);
% se = strel('square',10);

R_img=imgSeqColor(:,:,:,ref);
R_mean=boxfilter(R_img(:,:,1), r)./ N+boxfilter(R_img(:,:,2), r)./ N+boxfilter(R_img(:,:,3), r)./ N;
R_mean= R_mean./3;

R_var= boxfilter(R_img(:,:,1).*R_img(:,:,1), r)./ N+ boxfilter(R_img(:,:,2).*R_img(:,:,2), r)./ N+...
    boxfilter(R_img(:,:,3).*R_img(:,:,3), r)./ N- R_mean.* R_mean*3;

R_var = R_var./3;
R_std=sqrt(max(R_var,0));

i_mean(:,:,ref)=R_mean;
muIdxMap = R_mean < exposureThres | R_mean > 1 - exposureThres;
% R=2*r+1;
% countWindowt = ones(R, R);
% countWindow=countWindowt./sum(countWindowt(:));

for i = 1 : n
    if i ~= ref
        img= imgSeqColor(:,:,:,i);
        
        
        i_mean(:,:,i)=boxfilter(img(:,:,1), r)./ N+boxfilter(img(:,:,2), r)./ N+boxfilter(img(:,:,3), r)./ N;
        i_mean(:,:,i)= i_mean(:,:,i)./3;
        
        i_var= boxfilter(img(:,:,1).*img(:,:,1), r)./ N+ boxfilter(img(:,:,2).*img(:,:,2), r)./ N+...
            boxfilter(img(:,:,3).*img(:,:,3), r)./ N- i_mean(:,:,i).* i_mean(:,:,i)*3;
        
        i_var = i_var./3;
        i_std(:,:,i)=sqrt(max(i_var,0));
        
        
        
        mean_iR=boxfilter(img(:,:,1).*R_img(:,:,1), r)./ N+boxfilter(img(:,:,2).*R_img(:,:,2), r)./ N+...
            boxfilter(img(:,:,3).*R_img(:,:,3), r)./ N;
        mean_iR=mean_iR/3;
        cov_iR = mean_iR -  i_mean(:,:,i).* R_mean;
        
        %         mean_iR=boxfilter(img(:,:,1).*R_img(:,:,1), r)./ N+boxfilter(img(:,:,2).*R_img(:,:,2), r)./ N+...
        %             boxfilter(img(:,:,3).*R_img(:,:,3), r)./ N-i_mean(:,:,i).* R_mean*3;
        %         cov_iR=mean_iR./3;
        %
        
        sMap(:,:,i) = (cov_iR + C) ./ (R_std.* i_std(:,:,i) + C);
        sRefMap(:,:,i)=(sMap(:,:,i) >structureThres);
        
        
        t= sRefMap(:,:,i) ;
        t(muIdxMap) = 1;
        sRefMap(:,:,i) =t;
        
        sRefMap(:,:,i) = imopen(sRefMap(:,:,i),se);
        
        
        
        
        %     sRefMap(:,:,i) =bwareaopen( sRefMap(:,:,i),150);
        %  sRefMap(:,:,i) = imopen(sRefMap(:,:,i),se);
        
        
        
        %   figure,imshow( sRefMap(:,:,i))
        
        %   sRefMap(:,:,i) = imopen(sRefMap(:,:,i),se);
        
        %
        cMu  = imhistmatch(i_mean(:,:,i), R_mean,256);
        diff = abs(cMu - R_mean);
        cMap = diff <= consistencyThres;
        sRefMap(:,:,i) =  sRefMap(:,:,i).*cMap;
        
        
        
        
        
        
        %     figure,imshow( sRefMap(:,:,i))
        
        %     sRefMap(:,:,i) =bwareaopen( sRefMap(:,:,i),150);
        
        %         figure,imshow( sRefMap(:,:,i))
        
        %         sRefMap(:,:,i) = imopen(sRefMap(:,:,i),se);
        
        
        
        temp = imhistmatch(imgSeqColor(:,:,:,ref),...
            imgSeqColor(:,:,:,i),256);
        %         temp(temp<0) = 0;
        %         temp(temp>1) = 1;
        
        %         figure,imshow(temp)
        
        imgSeqcor(:,:,:,i)=imgSeqColor(:,:,:,i).*repmat(sRefMap(:,:,i),[1 1 3])+...
            temp.*repmat(1-sRefMap(:,:,i),[1 1 3]);
    end
    imgSeqcor(:,:,:,ref)=imgSeqColor(:,:,:,ref);
    
    
end
