function tImg=Afi2TImg(afiFileName)

% Convenience function to generate an instance of TissueImage given an
% afi file
% INPUT
% afiFileName - name of the afi file
% OUTPUT
% tImg - an instance of class TissueImage enabling us to interface with the
% image data corresponding to the afi file.

[afiPath,~,~]=fileparts(afiFileName);

fileList=cellfun(@(x) fullfile(afiPath,x),ReadAfiFilenames(afiFileName),...
    'Unif',false);
numberOfChannels=length(fileList);
infoList=cellfun(@(x) imfinfo(x),fileList,'Unif',false);
layerList=cellfun(@(x) find(strcmp({x.ColorType},'grayscale')),infoList,...
    'Unif',false);
layerWidths=cellfun(@(x) [x.Width],infoList,'Unif',false);
layerHeights=cellfun(@(x) [x.Height],infoList,'Unif',false);
if(~(isequal(layerList{:})&&isequal(layerWidths{:})&&isequal(layerHeights{:})))
   error('Mismatch in channel images'); 
end
layerList=layerList{1}(:);
numberOfLayers=length(layerList);

imgFileMat=repmat(fileList',numberOfLayers,1);
imgLayerMat=repmat(layerList,1,numberOfChannels);

tImg=TissueImage(imgFileMat,imgLayerMat);

end
