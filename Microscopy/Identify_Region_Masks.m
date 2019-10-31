function [regionMask,pixelRegions]=Identify_Region_Masks(tImg)

% Convenience function to identify parts of slide corresponding to
% different cellular areas. 
% Each slide contains 1 tissue section plus two control blocks.
% INPUT
% tImg - an IF slide specified as an instance of the TissueImage class
%               
% OUTPUT
% The background correction model is specified in terms o
%
% regionMask -  a numeric array of size 1/16 the original slide 
% taking values 0-BG, 1-tissue, 2/3- cell blocks denoting the different
% areas
% pixelRegion    - a cell array with 3 elements marking the image location of
% the bounding box of the image areas described above. The bounding box is
% specifed in the pixelRegion format used in the TissueRegion class.


img=tImg.LoadImage('MagLevel',tImg.numberOfMagLevels);
img=sum(img.^2,3);
img=img/max(img(:));
thresh=multithresh(img,4);
fullMask=img>thresh(1);




%fullMask=imbinarize(img,adaptthresh(img,0.01));

% remove small empty space between cells etc
se = strel('disk',25);
maskClose=imclose(fullMask,se);

%keep the 3 biggest regions
numberOfRegions=3;
props=regionprops(maskClose,{'Area','PixelIdxList','BoundingBox','centroid'});
[~,idx]=sort([props.Area],'descend');
%maskClose(vertcat(props(idx((numberOfRegions+1):end)).PixelIdxList))=false;
props(idx(4:end))=[];

%props=regionprops(maskClose,{'BoundingBox','centroid'});
centroids=cell2mat({props(:).Centroid}');
% Tumor assumed to be the right most
[~,tumorIdx]=max(centroids(:,1));
% Remaing two are the pellets
cellPelletIdx=setdiff(1:3,tumorIdx);

[~,topPelletIdx]=min(centroids(cellPelletIdx,2));
topPelletIdx=cellPelletIdx(topPelletIdx);
bottomPelletIdx=setdiff(cellPelletIdx,topPelletIdx);
regionOrder=[tumorIdx,topPelletIdx,bottomPelletIdx];
pixelRegions=cell(numberOfRegions,1);
regionMask=zeros(size(fullMask));
for regionCounter=1:numberOfRegions
   regionNumber=regionOrder(regionCounter); 
   rowRange=ceil(cumsum(props(regionNumber).BoundingBox([2,4]))); 
   rowRange(rowRange>size(img,1))=size(img,1);
   colRange=ceil(cumsum(props(regionNumber).BoundingBox([1,3])));
   colRange(colRange>size(img,2))=size(img,2);
   pixelRegions{regionCounter}={rowRange,colRange};
   regionMask(props(regionNumber).PixelIdxList)=regionCounter;
end

%regionOrder=[0 regionOrder];
%regionMask=regionOrder(bwlabel(maskClose)+1);
    
end