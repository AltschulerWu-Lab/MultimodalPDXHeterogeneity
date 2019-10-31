function [tissueRegionList,subregionPos]=TissueRegionList_Generator(tImg,subImgDims,outMagLevel)

% A function that takes a tissue-image, identifies the tumor portion, and breaks it up into a grid of sub-images.
% INPUTS
% tImg	      -	A instance of class TissueImage	or TissueImageCorrected pointing to the image of a tissue slide
% subImgDims  - A vector of length 2 denoting the dimensions in pixels (rows x cols) of each sub-image
% outMagLevel - The magnification level (in the image pyramid) at which calculations are performed
%
% OUTPUT
% tissueRegionList -	A struct array, with each element corresponding to a different sub-image, with fields containing information enabling the sub-image to be read using the TissueImage interface.
% 			tImg - the TissueImage instance which can be sued to read the image
% 			pixelRegion - the location in pixels of the sub-image, specified in the pixelRegion format used by TissueImage
%			magLevel - the magnification level at which the sub-image should be loaded.
%
% subregionPos    -	A 2D array of size, number_of_subimagesx2 containing the grid position of each sub-image. This function identifies the area of the image corresponding to the tissue
% 			and breaks it up into a grid of non-overlapping sub-images of size specified by subImgDims. The values in this array are the row and column number in that grid for 
% 			a specific sub-image
subImgXRes=subImgDims(2);
subImgYRes=subImgDims(1);
%regionIdx=1;% For Tumor


outImgDim=tImg.dimensions(outMagLevel,:);
inImgDim=tImg.dimensions(tImg.numberOfMagLevels,:);
magFactor=mean(outImgDim./inImgDim);

%Marker Scales and Intensity Cutoff

% Identify foreground at low res

%[~,pixelRegion]=Identify_Region_Masks(tImg);
%rowRange=round(magFactor*pixelRegion{regionIdx}{1});rowRange(rowRange>outImgDim(1))=outImgDim(1);
%colRange=round(magFactor*pixelRegion{regionIdx}{2});colRange(colRange>outImgDim(2))=outImgDim(2);

[~,pixelRegion]=Identify_Tumor_Mask(tImg);
rowRange=round(magFactor*pixelRegion{1});rowRange(rowRange>outImgDim(1))=outImgDim(1);
colRange=round(magFactor*pixelRegion{2});colRange(colRange>outImgDim(2))=outImgDim(2);



nY=ceil(((rowRange(2)-rowRange(1))+1)/subImgYRes);
nX=ceil(((colRange(2)-colRange(1))+1)/subImgXRes);

nSubRegions=nX*nY;
subregionCounter=1;
tissueRegionList=struct;
subregionPos=zeros(nSubRegions,2);
for x=1:nX
    for y=1:nY
        subregionPos(subregionCounter,:)=[y,x];
        subRowRange=(y-1)*subImgYRes+[0,subImgYRes-1]+rowRange(1);
        subRowRange(subRowRange<1)=1; subRowRange(subRowRange>outImgDim(1))=outImgDim(1);
        subColRange=(x-1)*subImgXRes+[0,subImgXRes-1]+colRange(1);
        subColRange(subColRange<1)=1; subColRange(subColRange>outImgDim(2))=outImgDim(2);
        tissueRegionList(subregionCounter).tImg=tImg;
        tissueRegionList(subregionCounter).pixelRegion={subRowRange,subColRange};
        tissueRegionList(subregionCounter).magLevel=outMagLevel;
        subregionCounter=subregionCounter+1;
    end
end
