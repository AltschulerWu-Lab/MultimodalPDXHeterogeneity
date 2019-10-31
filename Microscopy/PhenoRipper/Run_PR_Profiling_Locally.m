% A convenience script demonstrating how to profile the whole set of samples using PhenoRipper without submitting jobs to a cluster
addpath('../');
addpath('../../');
addpath('../../Growth_Curves/');
addpath('../../Pathology/');
%%
numberOfSamples=36;
numberOfMarkerSets=4;
numberOfModelGroups=10; %this isn't really used
modelsPerModelGroup=10; % Be careful about changing this since models are stored in groups of 10
params=GetParams('microscopy');

for jobCounter=1:144
    
    
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    temp=JobIdNumToTasks(jobCounter,...
        [numberOfSamples,numberOfMarkerSets,numberOfModelGroups]);
    sampleNumber=temp(1);
    markerSet=temp(2);
    modelGroup=temp(3);
    %resFile=fullfile(params.microscopy.rawPRDir,['PR_Results_MS'...
    %    num2str(markerSet) '_Sample' num2str(sampleNumber) '_Group1.mat']);
    resFile=fullfile(params.microscopy.bgSubtractedPRDir,['PR_Results_MS'...
        num2str(markerSet) '_Sample' num2str(sampleNumber) '_Group1.mat']);
    isNotFinished=true;
    if(exist(resFile,'File'))
        try
            resData=load(resFile);
            isNotFinished=any(cellfun(@isempty,resData.superblockProfiles));
        end
    end
    
    if(isNotFinished)
       SerializePhenoRipper(jobCounter);
    else
        warning([resFile ' skipped!']);
    end
end
