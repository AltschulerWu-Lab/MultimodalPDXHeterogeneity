function info =GetImagingInfo(imgFilename)
% Convenience function to pull various image parameters from image headers
allInfo=imfinfo(imgFilename);

infoSplit=regexp(strsplit(allInfo(1).ImageDescription,'|'),' = ','split');

info=struct;
for fieldCounter=1:length(infoSplit)
    if(length(infoSplit{fieldCounter})==2)
        
        fieldName=matlab.lang.makeValidName(infoSplit{fieldCounter}{1});
        fieldVal=str2double(infoSplit{fieldCounter}{2});
        if(~isnan(fieldVal))
            info.(fieldName)=fieldVal;
        else
            info.(fieldName)=infoSplit{fieldCounter}{2};
        end
    end
end



end