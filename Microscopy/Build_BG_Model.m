function [bgInterp,shadeInterp]=Build_BG_Model(tImg)

% Generates a model to perform background correction on IF Images
% INPUT
% tImg - an IF slide on which background correction is desired, specified in the form of an instance of the TissueImage class
%               
% OUTPUT
% The background correction model is specified in terms of two set of parameters, 
% each of which is calculated at the different magnification levels and for all 
% the different markers (these are therefore 3x4 cell arrays for our data). 
% Each interpolant models is specified in terms of griddedInterpolant functions 
% enabling imputation of intermediate values for background subtraction. 
%
% bgInterp -  This is meant to correct for overall background variation, 
% and is therefore a 2D interpolant.
% shadeInterp    - This is meant to control for the observed striping pattern, 
% and thus is a 1D interpolant along the x-axis (columns) of the image

imgLowRes=tImg.LoadImage('MagLevel',tImg.numberOfMagLevels);
intLowRes=sum(imgLowRes.^2,3);

intLowRes=intLowRes/max(intLowRes(:));
maskLowRes=imbinarize(intLowRes,adaptthresh(intLowRes,0.25));

se=strel('disk',30);
%se1=strel('line',yres_small/3,90);
maskLowRes=imfill(imclose(maskLowRes,se),'holes');
props=regionprops(~maskLowRes,'Area','PixelIdxList');
holes={props([props.Area]<2E4).PixelIdxList};
maskLowRes(vertcat(holes{:}))=true;
%imagesc(maskLowRes)
%%
imgSizeLow=size(maskLowRes);
blockSize=[30,98];
nBlocks=ceil(imgSizeLow./blockSize);
postMul=repmat((1:imgSizeLow(2))'<=blockSize(2),1,nBlocks(2));
for i=1:nBlocks(2)
    postMul(:,i)=circshift(postMul(:,i),(i-1)*blockSize(2));
end
preMul=repmat((1:imgSizeLow(1))<=blockSize(1),nBlocks(1),1);
for i=1:nBlocks(1)
    preMul(i,:)=circshift(preMul(i,:),(i-1)*blockSize(1));
end
%%
mags=bsxfun(@rdivide,tImg.dimensions,tImg.dimensions(end,:));
bgInterp=cell(tImg.numberOfMagLevels,tImg.numberOfChannels);
shadeInterp=cell(tImg.numberOfMagLevels,tImg.numberOfChannels);
stripeWidth=97.5;
xVals=repmat((1:size(imgLowRes,2))-floor((1:size(imgLowRes,2))/...
    stripeWidth)*stripeWidth,size(imgLowRes,1),1);
[uniqueX,~,ic]=unique(xVals);
for channelNum=1:tImg.numberOfChannels
    bgInt=imgLowRes(:,:,channelNum); bgInt(maskLowRes)=0;
    avgBgInt=(preMul*bgInt*postMul)./(preMul*double(~maskLowRes)*postMul);
    isFgBlock=(preMul*double(maskLowRes)*postMul/prod(blockSize))>0.25;
    
    
    %
    
    maxiter       = 20;
    tol           = 1e-5;
    param.lambda  = 10^9;   % weight on data fidelity (should usually be large).
    param.alpha   = 1;  % regularisation parameters \alpha.
    param.gamma   = 0.5;    % regularisation parameters \gamma.
    param.epsilon = 0.1;    % accuracy of Ambrosio-Tortorelli approximation of the edge set.
    
    [imputedImg,~]=inpainting_mumford_shah(avgBgInt,isFgBlock,maxiter,tol,param);
    %
    [blockY,blockX]=meshgrid(((1:nBlocks(2))-0.5)*blockSize(2),((1:nBlocks(1))-0.5)*blockSize(1));
    [Y,X]=meshgrid(1:imgSizeLow(2),1:imgSizeLow(1));
    interpLowRes=interp2(blockY,blockX,imputedImg,Y,X);
    
    %
    
    img2Stripe=imgLowRes(:,:,channelNum)./interpLowRes;
    img2Stripe(maskLowRes)=NaN;
    
    
    %stripeVals=zeros(length(uniqueX),1);
    %   stripeErr=zeros(length(uniqueX),1);
    %idxNan=find(isnan(img2Stripe));
    ic1=ic;
    ic1(isnan(img2Stripe(:)))=0;
    [temp,name]=grpstats(img2Stripe(:),ic1,{'mean','gname'});
    stripeVals=temp(~strcmp(name,'0'));
%     for i=1:length(uniqueX)
%         stripeVals(i)=mean(img2Stripe(ic1==i));
%         %       stripeErr(i)=nanstd(img2Stripe(ic==i));
%     end
    shadingImg=zeros(size(img2Stripe));
    shadingImg(:)=stripeVals(ic);
    
    
    
    
    for magCounter=1:tImg.numberOfMagLevels
        bgInterp{magCounter,channelNum}=griddedInterpolant(...
            (mags(magCounter,1)*blockX)-(mags(magCounter,1)-1)/2,...
            (mags(magCounter,2)*blockY)-(mags(magCounter,2)-1)/2,imputedImg);
        shadeInterp{magCounter,channelNum}=griddedInterpolant(...
            (mags(magCounter,2)*(1:size(shadingImg,2)))-mags(magCounter,2)/2,...
            shadingImg(1,:));
    end
end
end
