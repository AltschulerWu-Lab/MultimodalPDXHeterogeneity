%% A script that saves the (previously calculated) PhenoRipper profiling results for all 36 samples across 4 markers sets in a single results file 
%  While phenoripper results were calculated for multiple runs (differing by initial conditions) for each marker set, only one representative run is saved here.
addpath('../../');
addpath('../');
addpath('../../Growth_Curves/');
addpath('../../Pathology/');
addpath(genpath('../../Common/Plotting/'));
params=GetParams({'all'});

%%
[combinedProfiles,tRegionInfo,nFGBlocksByMarkerSet,replicateSectionNumbers,tRegionPos]=...
    Load_PR_Results();
%%

numberOfMarkerSets=length(params.microscopy.markerSets);

%% Load the background subtraction model
%temp=load(params.microscopy.bgSubtractedImgList);
temp=LoadScaledImgList();% BEWARE: this will force all paths to be relative to rootDir in getParams
tImgList=temp.scaledImageList;
sectionNumberList=temp.sectionNumbers;
%% Load all the tissue images for a given marker set


numberOfReplicates=3;
subImgDims=[1300,1300];
outMagLevel=1;
msTRegionLists=cell(numberOfMarkerSets,1);

for markerSet=1:numberOfMarkerSets
    tissueRegionList=[];
    subRegionPos=[];
    sampleNumbers=[];
    replicateNumbers=[];
    sectionNumbers=[];
    
    
    for sampleNumber=1:params.samples.numberOfSamples
        tic;
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
        toc;
    end
    msTRegionLists{markerSet}=tissueRegionList;
    
end



%% Load PR Models
modelGroupsPerMarkerSet=10;
modelsPerModelGroup=10;
prResultDir='/work/bioinformatics/srajar/Leidos/Data/PR_Results/';
prModels=cell(numberOfMarkerSets,1);
for markerSet=1:numberOfMarkerSets
    prModels{markerSet}=cell(modelGroupsPerMarkerSet*modelsPerModelGroup,1);
    for modelGroup=1:modelGroupsPerMarkerSet
        modelFile=fullfile(prResultDir, ['PR_Models-MS' num2str(markerSet) '-R' ...
            num2str(modelGroup) '.mat']);
        models=load(modelFile);
        for modelCounter=1:modelsPerModelGroup
            modelNumber=(modelGroup-1)*modelsPerModelGroup+modelCounter;
            prModels{markerSet}{modelNumber}=models.prModels{modelCounter};
        end
    end
end
%%
numberOfSamples=params.samples.numberOfSamples;
numberOfSuperBlocks=30;
runNumberToDeposit=3;

cleanResult=struct;

% Load pr models and store in cleanResult.prModel
cleanResult.prModel=cell(numberOfMarkerSets,1);

cleanResult.sampleProfiles=cell(numberOfMarkerSets,1);%done

% Section level profiling
cleanResult.sampleReplicateCombos=cell(numberOfMarkerSets,1); %done
cleanResult.sampleReplicateFilenames=cell(numberOfMarkerSets,1); %done
cleanResult.sampleReplicateProfiles=cell(numberOfMarkerSets,1); % done
% Region level profiling
cleanResult.sampleReplicateRegionNFG=cell(numberOfMarkerSets,1);%done
cleanResult.sampleReplicateRegionPixels=cell(numberOfMarkerSets,1);%done
cleanResult.sampleReplicateRegionProfiles=cell(numberOfMarkerSets,1);%done
for markerSet=1:numberOfMarkerSets
    cleanResult.prModel{markerSet}=prModels{markerSet}{runNumberToDeposit};
    
    profiles=combinedProfiles{markerSet}{runNumberToDeposit};
    cleanResult.sampleProfiles{markerSet}=zeros(numberOfSamples,numberOfSuperBlocks);
    
    cleanResult.sampleReplicateCombos{markerSet}=...
        zeros(numberOfSamples*numberOfReplicates,2);
    cleanResult.sampleReplicateFilenames{markerSet}=...
        cell(numberOfSamples,numberOfReplicates);
    
    cleanResult.sampleReplicateRegionNFG{markerSet}=cell(numberOfSamples,numberOfReplicates);
    cleanResult.sampleReplicateRegionProfiles{markerSet}=cell(numberOfSamples,numberOfReplicates);
    cleanResult.sampleReplicateRegionPixels{markerSet}=cell(numberOfSamples,numberOfReplicates);
    % Sample Level quantification
    for sampleNumber=1:numberOfSamples
        % load profiles belonging to sample
        isInSample=tRegionInfo(markerSet).sample==sampleNumber;
        cleanResult.sampleProfiles{markerSet}(sampleNumber,:)=...
            AverageProfile(profiles(isInSample,:),nFGBlocksByMarkerSet{markerSet}(isInSample));
        

        
        % Sample Rep level quantification
        
        
        for repNumber=1:numberOfReplicates
            counter=(sampleNumber-1)*numberOfReplicates+repNumber;
            cleanResult.sampleReplicateCombos{markerSet}(counter,1)=sampleNumber;
            cleanResult.sampleReplicateCombos{markerSet}(counter,2)=repNumber;
            
            % Find all profiles&nFG in sample and replicate, perform
            % averaging
            isInSection=isInSample&(replicateSectionNumbers{markerSet}==repNumber);
            cleanResult.sampleReplicateProfiles{markerSet}(counter,:)=...
                AverageProfile(profiles(isInSection,:),nFGBlocksByMarkerSet{markerSet}(isInSection));
            cleanResult.sampleReplicateFilenames{markerSet}{sampleNumber,repNumber}=...
                msTRegionLists{markerSet}(find(isInSection,1)).tImg.imgFileMat{1,1};
            
            
            
            %Region level profiling
            % Find number of regions in section image
            regionIdx=tRegionPos{markerSet}(isInSection,:);
            nrRegions= max(regionIdx,[],1);
            cleanResult.sampleReplicateRegionNFG{markerSet}{sampleNumber,repNumber}=zeros(nrRegions);
            cleanResult.sampleReplicateRegionProfiles{markerSet}...
                {sampleNumber,repNumber}=zeros(nrRegions(1),nrRegions(2),numberOfSuperBlocks);
            cleanResult.sampleReplicateRegionPixels{markerSet}{sampleNumber,repNumber}=cell(nrRegions);
            regionsInSection=find(isInSection);
            
            
            %loop over regions and fill in values for profiles, nFG,
            %pixelRegions
            for regionCounter=1:length(regionsInSection)
                r=regionsInSection(regionCounter);
                i=regionIdx(regionCounter,1);
                j=regionIdx(regionCounter,2);
                cleanResult.sampleReplicateRegionNFG{markerSet}{sampleNumber,repNumber}(i,j)=...
                    nFGBlocksByMarkerSet{markerSet}(r);
                cleanResult.sampleReplicateRegionProfiles...
                    {markerSet}{sampleNumber,repNumber}(i,j,:)=...
                    profiles(r,:);
                cleanResult.sampleReplicateRegionPixels{markerSet}{sampleNumber,repNumber}{i,j}=...
                    msTRegionLists{markerSet}(r).pixelRegion;
            end
            
            
            
        end
        
    end
    cleanResult.sampleReplicateFilenames{markerSet}=...
        regexprep(cleanResult.sampleReplicateFilenames{markerSet},params.microscopy.ToDepositDir,'');

end

%%
save(params.microscopy.processedResultsFile,'-struct','cleanResult');

%%
function profAvg=AverageProfile(profiles,weights)

weights=weights(:)/sum(weights);
profiles(weights==0,:)=[];
weights(weights==0)=[];
profAvg=sum(bsxfun(@times,profiles,weights),1);
end
