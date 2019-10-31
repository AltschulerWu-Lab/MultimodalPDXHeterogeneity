function [tImg]=RebaseTissueImage(tImg,newRoot,fixedLevels)
    
% Convenience function to change the filepath in an existing TissueImage
% The primary use is if we want to use saved tissueImages but the locations of the underlying 
% files has moved.
% INPUT
% tImg - an instance of class TissueImage or TissueImageCorrected
% newRoot - new root directory
% fixedLevels - an integer specifying low many levels of the path must be preserved. For example
% if old path was a/b/c/d.svs and new path is x/y/c/d.svs this would be 2 (since c&d are the same).
% OUTPUT
% tImg - a tissueImage with file paths pointing to a new root directory
    imgFileMat=tImg.imgFileMat;
    for r=1:size(imgFileMat,1)
        for c=1:size(imgFileMat,2)
            imgFileMat{r,c}=rebasePath(imgFileMat{r,c},newRoot,fixedLevels);
        end
    end
    tImg.imgFileMat=imgFileMat;
end

function [newPath]= rebasePath(origFile,newRoot,fixedLevels)
    fileParts=split(origFile,'/');
    preserved=fileParts(end-(fixedLevels-1):end);
    newPath=fullfile(newRoot,preserved{:});
end
