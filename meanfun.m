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
function Y=meanfun(X,X_mean,index,width,xielv)


switch index
    
    case 'direct'
       Y=X;
Y(Y>0.5)=1-Y(Y>0.5); 
       Y2=X_mean;
Y2(Y2>0.5)=1-Y2(Y2>0.5); 
Y=Y+Y2;
    case 'truncated'
if 0==exist('width','var')
    width=0.9;%????????
end
if 0==exist('xielv','var')
   
    xielv=0.7; %????????
end     
% if 0==exist('gamma','var')
%    
%     gamma=2; %????????
% end     
Mexp=1-(max(abs(X-0.5),0.0)-0.0)/(0.5-0.0);
Mexp=min(Mexp,width)/width;
m=0.5;
Y=(atan((Mexp-m).*(tan(pi*m-xielv*m)-tan(-pi*m+xielv*m)))+pi*m-xielv*m)./(pi-xielv); 

 width=0.95;
Mexp=1-(max(abs(X_mean-0.5),0.0)-0.0)/(0.5-0.0);
Mexp=min(Mexp,width)/width;
 xielv=0.05;
Y_g=(atan((Mexp-m).*(tan(pi*m-xielv*m)-tan(-pi*m+xielv*m)))+pi*m-xielv*m)./(pi-xielv);
Y=Y+Y_g;


    case 'gussian_lg'
        lSig=0.5;
       gSig=0.2;
      Y =  exp( -.5 * ( (X_mean - .5).^2 /gSig.^2 +  (X - .5).^2 /lSig.^2 ) );     
     case 'gussian_l'
         lSig=0.2;
         Y=exp((-(X-0.5).^2)./(lSig.^2*2));
    case 'nontrucated'    %% arctan
%      Y5=2/pi*atan(100*(X.*(X<=0.5)))+2/pi*atan(100*(X_mean.*(X_mean<=0.7)));
% Y6=2/pi*atan(100*((1-X).*(X>0.5)))+2/pi*atan(100*((1-X_mean).*(X_mean>0.3)));
% Y5=2/pi*atan(30*X.*(X<=0.5))+2/pi*atan(10*X_mean.*(X_mean<=0.2));
% Y6=2/pi*atan(30*(1-X).*(X>0.5))+2/pi*atan(10*(1-X_mean).*(X_mean>0.2));
% Y7=Y5+Y6;
% Y=Y7;

     Y5=2/pi*atan(20*(X.*(X<=0.5)));
Y6=2/pi*atan(20*((1-X).*(X>0.5)));
Y7=Y5+Y6;
Y=Y7;
    case 'gamma'
gamma=5;
% lSig=0.5;
y1=X.^gamma*(0.25.^(1-gamma)).*(X<=0.25);
 y2=(0.5-(0.5-X).^gamma*(0.25.^(1-gamma))).*(X>0.25&X<0.5);
  y3=(0.5-(0.5-(1-X)).^gamma*(0.25.^(1-gamma))).*(X>=0.5&X<0.75);
  y4=(1-X).^gamma*(0.25.^(1-gamma)).*(X>=0.75&X<=1);
  
  y5=X_mean.^gamma*(0.25.^(1-gamma)).*(X_mean<=0.25);
 y6=(0.5-(0.5-X_mean).^gamma*(0.25.^(1-gamma))).*(X_mean>0.25&X_mean<0.5);
  y7=(0.5-(0.5-(1-X_mean)).^gamma*(0.25.^(1-gamma))).*(X_mean>=0.5&X_mean<0.75);
  y8=(1-X_mean).^gamma*(0.25.^(1-gamma)).*(X_mean>=0.75&X_mean<=1);
%   y5=exp((-(X_mean-0.5).^2)./(lSig.^2*2));
y=y1+y2+y3+y4+y5+y6+y7+y8;
% y=y1+y2+y3+y4;
Y=2*y;
% Y=Y7./repmat(max(Y7),size(Y7,1),1);
end

% Y=(atan((Mexp-m).*(tan(pi*m-xielv*m)-tan(-pi*m+xielv*m)))+pi*m-xielv*m)./(pi-xielv);
% Y=exp(-(X-0.5).^2)./(lSig.^2*2);
% Y =  exp( -.5 * ( (X_mean - .5).^2 /gSig.^2 +  (X - .5).^2 /lSig.^2 ) ); 



%% huatu
%  width=0.4;%????????
%  xielv=0.8; %????????
% X=0:0.01:1;
% Mexp=1-(max(abs(X-0.5),0.0)-0.0)/(0.5-0.0);
% Mexp=min(Mexp,width)/width;
% m=0.5;
% Y=(atan((Mexp-m).*(tan(pi*m-xielv*m)-tan(-pi*m+xielv*m)))+pi*m-xielv*m)./(pi-xielv);
% plot(X,Y)


% sig=0.2;
% Y=exp(-(X-0.5).^2)./(sig.^2*2);

% Y=X;
% Y(Y>0.5)=1-Y(Y>0.5);
% plot(X,Y)