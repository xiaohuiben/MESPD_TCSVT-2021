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

function [fI3, i_meant,aa,N1] = scale_coarse(imgSeqColor,r,lambda)

[h,w,n]=size(imgSeqColor);
N = boxfilter(ones(h, w), r);

tem=ones(h, w);
tem(:,2:2:w)=0;
tem(2:2:h,:)=0;
N1= boxfilter(tem, r);

% R=2*r+1;
% countWindowt = ones(R, R);
% countWindow=countWindowt./sum(countWindowt(:));

p=6;

[WD, Cmax,i_mean]= weight_spd_detail(imgSeqColor,p,r,lambda);
WD=WD.*repmat(Cmax,[1 1 n]);


% WD= weight_cal_gray(imgSeqColor,[0 0 0 0 1 0 0]);

% brightness='direct';
% brightness='truncated';
% brightness='gussian_lg';
% brightness='gussian_l';
% brightness='nontrucated';
brightness='gamma';

WB= weight_cal_base(i_mean,brightness);

% WB_mask= weight_cal_gray(imgSeqColor,[0 0 0 0 0 0 1],brightness);
% WB=WB.*WB_mask;


F_temp2_detail=zeros(h,w,n);
F_temp2_base=zeros(h,w,n);

i_meant=zeros(ceil(h/2),ceil(w/2),2);

%% approximate aggregation through averaging(mean filter) the weight map


tic
for i = 1:n

    
%     aa=i_mean(:,:,i).*tem;
%     i_meant(:,:,i)=aa(1:2:h,1:2:w);
% 
%     W_D2=boxfilter(WD(:,:,i).*tem, r)./ N1;
%     W_B2=boxfilter(WB(:,:,i).*tem, r)./ N1;
%    
%    F_temp2_detail(:,:,i)=t*(W_D2.*(imgSeqColor(:,:,i)-i_mean(:,:,i)));
%    
%     F_temp2_base(:,:,i)= W_B2.* i_mean(:, :, i);
    
    
        aa=i_mean(:,:,i).*tem;
    i_meant(:,:,i)=aa(1:2:h,1:2:w);

%     W_D1=boxfilter(aa.*WD(:,:,i), r)./ N1;
%     W_D2=boxfilter(WD(:,:,i).*tem, r)./ N1;
%     W_B2=boxfilter(aa.*WB(:,:,i), r)./ N1;

    W_D1=boxfilter(i_mean(:,:,i).*WD(:,:,i), r)./ N;
    W_D2=boxfilter(WD(:,:,i), r)./ N;
    W_B2=boxfilter(i_mean(:,:,i).*WB(:,:,i), r)./ N;
    F_temp2_detail(:,:,i)=W_D2.*(imgSeqColor(:,:,i))-W_D1;
    F_temp2_base(:,:,i)= W_B2;  
    
    
%          W_D2=boxfilter(WD(:,:,i), r)./ N;
%     W_B2=boxfilter(WB(:,:,i), r)./ N;
%     F_temp2_detail(:,:,i)=W_D2.*(imgSeqColor(:,:,i)-i_mean(:,:,i));
%     F_temp2_base(:,:,i)= W_B2.* i_mean(:, :, i);
%     
    
    
end

 F_temp2=F_temp2_detail+F_temp2_base;

fI3=sum(F_temp2,3);
toc

% figure, imshow(fI3,[])


end



