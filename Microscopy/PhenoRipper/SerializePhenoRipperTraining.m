function SerializePhenoRipperTraining(runId)

% Cluster friendly function to generate multiple PhenoRipper models for
% each marker set with differing random seed. These models will serve as
% the basis for PhenoRipper profiling, and so having multiple models allows
% us to study the effect of randomness of PhenoRipper profiles.
% INPUT
% runId   -     An integer (supplied by the cluster), which serves as a
%               unique identifier for a single job. Each job produces multiple
%               PhenoRipper models for a single marker set, with the mapping from runID to 
%               markerSet and run within markerSet determined by the 
%               threadsPerMarkerSet & runsPerThread variables as in the formula below.
% OUTPUT
% Each run results in a single PhenoRipper model stored as a mat file in
% the directory specified in params.microscopy.bgSubtractedPRDir if
% cleanImage=true or params.microscopy.rawPRDir if cleanImage=false.

cleanImage=true;% set true if you want to perform background subtraction
threadsPerMarkerSet=10;
runsPerThread=10;
markerSet=ceil(runId/threadsPerMarkerSet);
threadNumber=rem(runId-1,threadsPerMarkerSet)+1;
disp(['Run ID:' num2str(runId) ' markerSet:' num2str(markerSet) ...
    ' threadNumber:' num2str(threadNumber)])
rng('shuffle');



addpath('../');
addpath('../../');
addpath('../../Growth_Curves/');
addpath('../../Pathology/');
addpath(genpath('../../Common/'));
addpath('engine/');
addpath('util/');
params=GetParams({'samples','microscopy'});
if(cleanImage)
    saveDir=params.microscopy.bgSubtractedPRDir;
    %temp=load(params.microscopy.bgSubtractedImgList);
    temp=LoadScaledImgList();
    tImgList=temp.scaledImageList;
    markerScales=params.microscopy.prBgSubMarkerScales;
else
    saveDir=params.microscopy.rawPRDir;
    temp=load(params.microscopy.rawImgList);
    tImgList=temp.rawImageList;
    markerScales=params.microscopy.prRawMarkerScales;
end
sectionNumberList=temp.sectionNumbers;
clear('temp');
%%
numberOfSamples=params.samples.numberOfSamples;
numberOfReplicates=3;

outMagLevel=1;
% blockSize=20;
% cutoffIntensity=10;
% subImgDims=[1300,1300];
% numberOfTrainingImages=100;
% if(cleanImage)
%     markerScales=repmat([0,30000],4,1);
% else
%     markerScales=repmat([5000,30000],4,1);
% end


blockSize=params.microscopy.prBlockSize;
cutoffIntensity=params.microscopy.prCutoffIntensity;
subImgDims=params.microscopy.prSubImgDims;
numberOfTrainingImages=params.samples.numberOfSamples;


tissueRegionList=[];
subRegionPos=[];
sampleNumbers=[];
replicateNumbers=[];
sectionNumbers=[];

for sampleNumber=1:numberOfSamples
    % tic;
    %[afiFiles,fileInfo]=GetImages('IF','sampleNumber',sampleNumber,...
    %    'markerSet',markerSet);
    
    for repCounter=1:numberOfReplicates
        tImg=tImgList{markerSet,sampleNumber}{repCounter};
        sectionNumber=sectionNumberList{markerSet,sampleNumber}(repCounter);
        [tissueRegionListTemp,subRegionPosTemp]=TissueRegionList_Generator(tImg,...
            subImgDims,outMagLevel);
        if(isempty(tissueRegionList))
            tissueRegionList=tissueRegionListTemp;
        else
            tissueRegionList=[tissueRegionList,tissueRegionListTemp];
        end
        subRegionPos=[subRegionPos;subRegionPosTemp];
        sampleNumbers=[sampleNumbers;sampleNumber*ones(length(tissueRegionListTemp),1)];
        replicateNumbers=[replicateNumbers;repCounter*ones(length(tissueRegionListTemp),1)];
        sectionNumbers=[sectionNumbers;sectionNumber*...
            ones(length(tissueRegionListTemp),1)];
    end
    
    %toc;
end



disp(['Training Marker Set ' num2str(markerSet) ' Started!']);
trainingRegionLists=cell(runsPerThread,1);
prModels=cell(runsPerThread,1);
for runNumber=1:runsPerThread
    trainingRegionLists{runNumber}=tissueRegionList(...
        randsample(length(tissueRegionList),numberOfTrainingImages));
end
parfor runNumber=1:runsPerThread
    
    prModels{runNumber}=TrainPR(trainingRegionLists{runNumber},...
        blockSize,markerScales,cutoffIntensity);
    disp(['Profiling Marker Set ' num2str(markerSet) ' Run ' num2str(runNumber) ' done!']);
    
    
end


saveFile=fullfile(saveDir, ['PR_Models-MS' num2str(markerSet) '-R' num2str(threadNumber) '.mat']);
save(saveFile,'prModels','cleanImage','subImgDims','outMagLevel',...
    '-v7.3');


end

