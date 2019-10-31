function [tumorMask,pixelRegion]=Identify_Tumor_Mask(tImg,closeRadius)

% Convenience function to identify parts of slide corresponding to
% the PDX tumor 
% Each slide contains 1 tissue section plus two control blocks.
% INPUT
% tImg - an IF slide specified as an instance of the TissueImage class
% closeRadius - the radius (on maglevel 3) used for image closing to avoid holes
%               messing up determination of tissue areas
% OUTPUT
% The background correction model is specified in terms o
%
% tumorMask -  a binary mask of size 1/16 the original image marking the
%              PDX tissue area
% pixelRegion    - a cell array  marking the bounding box of the 
%               tissue region specifed in the pixelRegion format used in 
%               the TissueRegion class.

if(nargin<2)
    closeRadius=25;
end
img=tImg.LoadImage('MagLevel',tImg.numberOfMagLevels);
img=sum(img.^2,3);
img=img/max(img(:));
thresh=multithresh(img,4);
fullMask=img>thresh(1);




%fullMask=imbinarize(img,adaptthresh(img,0.01));

% remove small empty space between cells etc
areaOfTumor=0;
%closeRadius=25;
minTumorArea=1E5;
while(areaOfTumor<minTumorArea)
se = strel('disk',closeRadius);
maskClose=imclose(fullMask,se);

props=regionprops(maskClose,{'Area','PixelIdxList','BoundingBox'});
[areaOfTumor,idx]=max([props.Area]);
closeRadius=closeRadius+10;
end
%maskClose(vertcat(props(idx((numberOfRegions+1):end)).PixelIdxList))=false;
props=props(idx);
tumorMask=false(size(fullMask));
tumorMask(props.PixelIdxList)=true;
rowRange=ceil(cumsum(props.BoundingBox([2,4]))); 
rowRange(rowRange>size(img,1))=size(img,1);
colRange=ceil(cumsum(props.BoundingBox([1,3])));
colRange(colRange>size(img,2))=size(img,2);
pixelRegion={rowRange,colRange};


end