function imagescNaN(dataMat,minMax,cmap,nanColor)
if(isempty(minMax))
   minMax=[nanmin(dataMat(:)),nanmax(dataMat(:))];
end
nColors=size(cmap,1);
scaledVals=(dataMat-minMax(1))/(minMax(2)-minMax(1));
scaledVals(scaledVals<0)=0;
scaledVals(scaledVals>1)=1;

colorIdx=round(scaledVals*(nColors-1))+1;
colorIdx(isnan(colorIdx))=nColors+1;
cmapPlus=[cmap;nanColor];
colorMat=zeros(size(dataMat,1),size(dataMat,2),3);
for i=1:3
    cmapChannel=cmapPlus(:,i);
    colorMat(:,:,i)=cmapChannel(colorIdx);
   
end
image(colorMat);

end