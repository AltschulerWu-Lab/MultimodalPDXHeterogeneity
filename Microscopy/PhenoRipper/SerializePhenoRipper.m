function SerializePhenoRipper(idNum)

% Cluster friendly function to generate PhenoRipper profiles for each of
% 36 samples imaged across 3 replicate sections with 4 marker sets 
% based on 10 PhenoRipper models (differing only by random seed) per marker set. 
% These profiles serve as the basis for comparing image similarity.
% INPUT
% idNum   -     An integer (supplied by the cluster), which serves as a
%               unique identifier for a single job. Usage for a 
%               local (i.e. non-cluster) is demonstrated in the script
%               Run_PR_Profiling_Locally.m. Essentially, the idNum needs 
%               to take values from 1 to 144 to generate the full set of 
%               profiles, with each idNum generating profiles for a single 
%               sample+marker_set combination.
%              
% OUTPUT
% Each run generates profiles for a single sample and marker set.Since each
% sample is represented by three replicate sections and 10 models are used
% for each marker set, 30 sets of profiles are generated. These will be stored as .mat
% files in the the directory specified by params.microscopy.bgSubtractedPRDir 
% (or rawPRDir if cleanImage = false) in GetParams.m.
% cleanImage=true or params.microscopy.rawPRDir if cleanImage=false.
cleanImage=true;

%% Setup
numberOfSamples=36;
numberOfMarkerSets=4;
numberOfModelGroups=10; %this isn't really used
modelsPerModelGroup=10; % Be careful about changing this since models are stored in groups of 10


temp=JobIdNumToTasks(idNum,...
    [numberOfSamples,numberOfMarkerSets,numberOfModelGroups]);
sampleNumber=temp(1);
markerSet=temp(2);
modelGroup=temp(3);

disp(['Run ID:' num2str(idNum) ' markerSet:' num2str(markerSet) ...
    ' sampleNumber:' num2str(sampleNumber) ' modelGroup:' num2str(modelGroup) ]);
%saveDir='/work/bioinformatics/srajar/Leidos/Data/PR_Results/';

addpath('../');
addpath('../../');
addpath('../../Growth_Curves/');
addpath('../../Pathology/');
addpath(genpath('../../Common/'));
addpath('engine/');
addpath('util/');

%temp=load('../scaledImages.mat');
%tImgList=temp.scaledImageList;
params=GetParams({'samples','microscopy'});
if(cleanImage)
    saveDir=params.microscopy.bgSubtractedPRDir;
    %temp=load(params.microscopy.bgSubtractedImgList);
    temp=LoadScaledImgList();
    tImgList=temp.scaledImageList;
   
else
    saveDir=params.microscopy.rawPRDir;
    temp=load(params.microscopy.rawImgList);
    tImgList=temp.rawImageList;
end
sectionNumberList=temp.sectionNumbers;
clear('temp');
%% Generate TissueRegionList for Sample


numberOfReplicates=3;
outMagLevel=1;

subImgDims=params.microscopy.prSubImgDims;

tissueRegionList=[];
subRegionPos=[];
sampleNumbers=[];
replicateNumbers=[];
sectionNumbers=[];
disp(['Generating TissueRegionList MS:' num2str(markerSet) ' Sample ' ...
    num2str(sampleNumber)]);
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

%% Profiling

superblockProfiles=cell(modelsPerModelGroup,1);
superblockImages=cell(modelsPerModelGroup,1);
nFGBlocks=cell(modelsPerModelGroup,1);



modelFile=fullfile(saveDir, ['PR_Models-MS' num2str(markerSet) '-R'...
    num2str(modelGroup) '.mat']);
modelList=load(modelFile);

for numberInGroup=1:modelsPerModelGroup
    tic;
    prModel=modelList.prModels{numberInGroup};
    
    [superblockProfiles{numberInGroup},~,nFGBlocks{numberInGroup},prResults]=...
        profilePR(tissueRegionList,prModel);
    
    superblockImages{numberInGroup}=cellfun(@(x) x.image_superblock_states,...
        prResults,'Unif',false);
    
    saveFile=['PR_Results_MS' num2str(markerSet) ...
        '_Sample' num2str(sampleNumber)...
        '_Group' num2str(modelGroup) '.mat'];
    save(fullfile(saveDir,saveFile),'superblockProfiles','nFGBlocks',...
        'sectionNumbers','replicateNumbers','sampleNumbers','subRegionPos',...
        'superblockImages');
    
    disp(['Finished profiling Marker Set ' num2str(markerSet) ' Sample ' ...
        num2str(sampleNumber) ' ModelGroup ' num2str(modelGroup) ...
        ' Model ' num2str(numberInGroup)]);
    toc;
end



disp(['Marker Set ' num2str(markerSet) ' Sample' num2str(sampleNumber) ' done!']);

end

