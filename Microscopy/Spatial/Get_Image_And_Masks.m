function [img,masks]=Get_Image_And_Masks(tImg)
	
% Convenience function to get the portion of a slide image corresponding to tumor and matching
% nuclear and cytoplasmic masks.
% INPUT
% tImg   - A slide from which we want the tumor image, passed as as instance of the TissueImage or TissueImageCorrected class.
%              
% OUTPUT
% img - A 3D matrix representing the portion of the slide corresponding to tumor tissue. i
% 	The 3rd dimension is the channel axis.
% masks - A struct with fields corresponding to nuclear/cellular/cytoplasmic/stromal masks matching img.
    masks=struct;
    [tumorMaskLow,pixelRegion]=Identify_Tumor_Mask(tImg);
    
    pixelRegionFull=cellfun(@(x) [16*x(1),16*x(2)],pixelRegion,'Unif',false);
    img=tImg.LoadImage('PixelRegion',pixelRegionFull);
    disp('Image Loaded');
    
    % Load Masks
    dapiImg=img(:,:,1)/max(max(img(:,:,1)));
    masks.dapiMask=imbinarize(dapiImg,adaptthresh(dapiImg,0.5));
    
    vimImg=img(:,:,4)/max(max(img(:,:,4)));
    vimImgSmooth=imgaussfilt(vimImg,10);
    masks.stromalMask=imbinarize(vimImgSmooth,adaptthresh(vimImgSmooth,0.01));
    
    tumorMaskLowCropped=tumorMaskLow(pixelRegion{1}(1):pixelRegion{1}(2),pixelRegion{2}(1):pixelRegion{2}(2));
    masks.tumorMask=imresize(tumorMaskLowCropped,size(dapiImg));
    
    allChannels=sum(img,3);
    allChannels=allChannels/max(max(allChannels));
    masks.cellularMask=imbinarize(allChannels,adaptthresh(allChannels,0.9));
  
    masks.dapiMask=masks.dapiMask&masks.tumorMask;
    se=strel('disk',25);
    
    masks.extraNuclearMask=imdilate(masks.dapiMask,se);
    
    masks.stromalMask=masks.stromalMask&masks.tumorMask;
    masks.cellularMask=masks.cellularMask&masks.tumorMask& ~masks.stromalMask;
    masks.cytoplasmicMask=masks.cellularMask&masks.extraNuclearMask&~masks.dapiMask & ~masks.stromalMask;
    masks.nuclearMask=masks.dapiMask&~masks.stromalMask;
    disp('Masks Calculated');
    

end
