
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
%----------------------------------------------------------------------
% This is a static scene implementation of "Detail-Preserving Multi-Exposure Fusion with Edge-Preserving Structural Patch Decomposition"
% Please refer to the following paper:
% H. Li et al., "Detail-Preserving Multi-Exposure Fusion with Edge-Preserving Structural Patch Decomposition, 2021" In press
% IEEE Transactions on Circuits and Systems for Video Technology
% Please kindly report any suggestions or corrections to xiaohui102788@126.com
%----------------------------------------------------------------------

clear ;
close all;
addpath(genpath(pwd));


folderNames = {'Arno_Bartlomiej Okonek'};

for i = 1:size(folderNames,2)
%     for i =1
    
    folderName = folderNames{i};
    
    Dir =sprintf('%s/',folderName);   %static scenes
    %     Dir =sprintf('C:/Users/hui li/Downloads/Dataset_Part1/Dataset_Part1/%s/',folderName);  %static scenes
    
    %       Dir =sprintf(' E:/20150816????????????/TIP_MSPD/TIP_MSPDcode/image_sequence(downsampled)/%s/',folderName);  %static scenes
    imgSeqColor= loadImg(Dir); % [0,1]
    imgSeqColor = downSample(imgSeqColor, 1024);
    imgSeqColor2 = uint8(imgSeqColor*255); % use im2double
%     figure,imshow(imgSeqColor(:,:,:,1))
%     figure,imshow(imgSeqColor(:,:,:,2))
%     figure,imshow(imgSeqColor(:,:,:,3))
    
    %    tic
    
    
    r1=4;
   
%     
    lambda=0.25;
    T=0.5;
% scale=1000;
% T=1;
    %% single scale
    % C_out= SPD_fast_single3(imgSeqColor,r,t,eps);
    %% multi-scale scale
    [ D1,i_mean1,aa1,N1] = scale_fine(imgSeqColor,r1,lambda);
%     figure,imshow(i_mean1(:,:,1))
%     figure,imshow(i_mean1(:,:,2))
%     figure,imshow(i_mean1(:,:,3))
    
    [w,h,~,~]=size(imgSeqColor);
    nlev = floor(log(min(w,h)) / log(2))-5;
    
    D2 = cell(nlev,1);
    aa2= cell(nlev,1);
    N2= cell(nlev,1);
    r2=4;
   

    lambda=lambda*T;
    
    for ii=1:nlev
        [ D2{ii},i_mean2,aa2{ii},N2{ii}] = scale_interm(i_mean1,r2,lambda);
        i_mean1=i_mean2;
        lambda=lambda*T;
%         figure,imshow(i_mean2(:,:,1))
%         figure,imshow(i_mean2(:,:,2))
%         figure,imshow(i_mean2(:,:,3))
    end
    
    
    %% the coarsest  scale
    r3=3;
%     t=1;
%     scale=0.005;
    [fI3,i_mean3,aa3,N3] = scale_coarse(i_mean2,r3,lambda);
    
    %% reconstruct
    %% Intermediate layers
    for ii=nlev:-1:1
        temp=aa2{ii};
        fI=zeros(size(temp));
        fI(1:2:size(temp,1),1:2:size(temp,2))=fI3;
        B2=boxfilter(fI, r2)./ N2{ii}+D2{ii};
        
        fI3=B2;
    end
    %% finest layers
    fI=zeros(size(aa1));
    fI(1:2:size(aa1,1),1:2:size(aa1,2))=B2;
    B1=boxfilter(fI, r1)./ N1;
    % C_out=repmat(B1,[1 1 3])+2/pi*atan(2.3*D1);
    C_out=repmat(B1,[1 1 3])+D1;
    toc
    
    %     % scale=1;
    %     % for ii=1:scale
    %     % [ fI2_detail(:,:,ii),i_mean2] = scale2(i_mean1,r,t);
    %     % i_mean1=i_mean2;
    %     % end
    %     % fI2_detail_sum=sum(fI2_detail,3);
    %
    %
    %     r3=5;
    %     t=1;
    %     scale=0.08;
    %     [fI3,i_mean3,aa3,N3] = scale3(i_mean2,r3,t,scale);
    %
    %     %% reconstruct
    %
    %     fI=zeros(size(aa2));
    %     fI(1:2:size(aa2,1),1:2:size(aa2,2))=fI3;
    %     B2=boxfilter(fI, r2)./ N2+D2;
    %
    %     fI=zeros(size(aa1));
    %     fI(1:2:size(aa1,1),1:2:size(aa1,2))=B2;
    %     B1=boxfilter(fI, r1)./ N1;
    %     C_out=repmat(B1,[1 1 3])+D1;
    
    % C_out=repmat(fI2,[1 1 3])+repmat(fI2_detail_sum,[1 1 3])+fI1_detail;
    
    % C_out=repmat(fI2,[1 1 3])+repmat(fI2_detail_sum,[1 1 3])+fI1_detail;
    
    %  toc
    figure,imshow(C_out,[])
    
    
    
    %% estimating images
    %     filename = sprintf('%s',folderName);
    %     %  imwrite(C_out,strcat('./exr_result',filename,'.png'));
    %     filename = ['/Users/huili/Desktop/Multi_mean_ifv/comparison methods/ADF/res_adf_test/',filepath1(i).name(1:end-4),'_adf',filepath1(i).name(end-5:end)];
    %     imwrite(C_out,strcat(filename,'.png'));
    
    
    filename = [folderName,'_mespd','.png'];
    %  imwrite(C_out,strcat('./exr_result',filename,'.png'));
    imwrite(C_out,filename);
    
    %% calculating objective metrics
    [s1, s2, s3, s4] = size(imgSeqColor2);
    imgSeq = zeros(s1, s2, s4);
    for ii = 1:s4
        imgSeq(:, :, ii) =  rgb2gray( squeeze( imgSeqColor2(:,:,:,ii) ) ); % color to gray conversion
    end
    Dir =sprintf('/Users/huili/Desktop/TMSD/TIP_MSPDcode/fast_multi-scale/MESPD_TCSVT/%s',filename);
    
    
    fI1 = imread(Dir);
    fI1 = double(rgb2gray(fI1));
    [Q(1), Qs1, QMap1] = mef_ms_ssim(imgSeq, fI1);
    %     score=Q(1)
    %
    score(1,i)=roundn(Q(1),-3);
    
    
end
roundn(mean(score),-3)