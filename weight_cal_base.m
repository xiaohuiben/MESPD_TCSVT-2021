
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
function W= weight_cal_base(I,brightness)

r = size(I,1);
c = size(I,2);
N = size(I,3);


W = ones(r,c,N);

for i = 1:N
   

     img=I(:,:,i);
     
     
 M=ones(size(I,1),size(I,2))*mean(img(:));
 W(:,:,i)=meanfun(img,M,brightness);
end



%normalize weights: make sure that weights sum to one for each pixel
W = W + 1e-12; %avoids division by zero
W = W./repmat(sum(W,3),[1 1 N]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





