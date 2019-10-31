function img=read_and_scale_timg(tissueRegion,xres,yres,markerScales)
% Convenience function used to return a rescaled image for a tissueImage corresponding to the size and marker-scales provided.

img=zeros(xres,yres,tissueRegion.tImg.numberOfChannels);
tempImg=tissueRegion.tImg.LoadImage('MagLevel',tissueRegion.magLevel,...
    'pixelRegion',tissueRegion.pixelRegion);
nR=min(xres,size(tempImg,1));
nC=min(yres,size(tempImg,2));
img(1:nR,1:nC,:)=tempImg(1:nR,1:nC,:);
for c=1:size(img,3)
    temp=img(:,:,c);
    temp=100*(temp-markerScales(c,1))/(diff(markerScales(c,:)));
    temp(temp<0)=0;temp(temp>100)=100;
    img(:,:,c)=temp;
end

end
