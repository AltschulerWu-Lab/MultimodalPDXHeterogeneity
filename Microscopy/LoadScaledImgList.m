function scaledImgList=LoadScaledImgList()

% Convenience function to load contents of params.microscopy.bgSubtractedImgList
% The primary purpose of this function (As opposed to simply using load) is to
% force all paths to be relative to params.microscopy.ToDepositDir. 
% INPUT
% OUTPUT
% scaledImgList - contents of params.microscopy.bgSubtractedImgList with all tissue images 
% 		  with root directory reset to params.microscopy.ToDepositDir

params=GetParams('microscopy');
scaledImgList=load(params.microscopy.bgSubtractedImgList);
s=scaledImgList.scaledImageList;
for a=1:size(s,1)
    for b=1:size(s,2)
        for c=1:length(s{a,b})
            s{a,b}{c}=RebaseTissueImage(s{a,b}{c},params.microscopy.ToDepositDir,2);
        end
    end
    scaledImgList.scaledImageList=s;
    
end
