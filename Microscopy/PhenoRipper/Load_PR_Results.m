function [combinedProfiles,tRegionInfo,nFGBlocksByMarkerSet,...
    replicateSectionNumbers,tRegionPos]=Load_PR_Results(useBGSubtractedResults)

% Function to load up and organize saved PhenoRipper results. 
% The function assumes PhenoRipper results for all 36 samples across 4 marker sets
% and multiple "random" seeded PhenoRipper runs have been saved in params.microscopy.bgSubtractedPRDir
% and named as indicated in SerializePhenoRipper.m. This function is purely
% meant to load up and organize those saved results.  
% INPUT
% useBGSubtractedResults   -     logical variable indicating whether background subtracted results 
% 				should be used. As raw results were not deposited, this should be set to true.
% OUTPUT
% All results are returned as organized by IF Marker Set (i.e. each contains 4 elements corresponding
% to different marker-sets). Additionally, PhenoRipper was run by dividing each slide into sub-images, and thus
% one (30-dimensional) PhenoRipper profiles is generated for each sub-image from a section stained with 
% a specific IF marker set. Note that each of the 36 regions (4 models x 3 replicate tumors x 3 regions per tumor)
% is represented by 3 replicate sections for a marker set. Thus this function also returns the mapping between
% the sub-image profiles and these variables.
%
% combinedProfiles	    -	Cell array with 4 elements (1 per marker set) containing the phenoRipper profiles
% 				for all sub-images from that marker set, across multiple PhenoRipper runs. 
% 				Each element is a cell array containing the profiles for a given PhenoRipper run 
% 				(i.e. starting with a random initial condition). The results of all runs for a marker set are organized identically 
% 				(i.e. they are applied the same-subimages in the same order, and hence other variable below reference each sub-image only once) 
% 				and stored as a 2D matrix of size number_of_subimages x 30, where 30 is the length of a phenoRipper profile.
%
% tRegionInfo		    -	Struct array with 4 elements (1 per marker set) containing information on the different 
% 				sub-images for that marker set. Each element contains 4 fields each of which is a vector
% 				of length number_of_subimages for that marker set, and contains information on the corresponding 
% 				aspect of the sub-image. The fields are
% 				1) sample - Takes values from 1 to 36 indicating which sample the subimage is from
% 				2) model - Takes values from 1 to 4 indicating which PDX model the sub-image was derived from
% 				3) tumor - Takes values between 1 to 3 indicating the which of the replicate tumors for a model the sub-image was derived from
% 				4) region - Takes values between 1 to 3 indicating the which of the three tumor region dorscal/ventral/ventral the sub-image was derived from
% 				5) section - Takes values between 1 and 108 =(36 samples x 3 replicates) for each section for a marker set that the sub-image was derived from.
% 				Note: the variable replicateSectionNumber below takes values 1-3 to denote which of the 3 replicate section for a region the sub-image was derived from.
%
%
% nFGBlocksByMarkerSet	    -	Cell array with 4 elements (1 per marker set) denoting the number of foreground blocks in each sub-images 
% 				from that marker set. Each element is a vector of length number_of_subimages in that marker set, containing
% 				the number or blocks that were considered as foreground (this is used to weight the different subimages)
%
% replicateSectionNumbers   -	Cell array with 4 elements (1 per marker set), each a vector of length number_of_subimages, denoting which 
% 				of the 3 replicates for a sample the sub-image was derived from.
%
% tRegionPos 		    -	Cell array with 4 elements (1 per marker set), each element being a 2D matrix of size number_of_subimages x 2
% 				that denote the position of the sub-image within the slide. Each PDX tissue region is broken into a grid of
% 				sub-images, and the two elements for each sub-image denote the row and column position within this grid.
if(nargin>1)
   useBGSubtractedResults=true; 
end
params=GetParams('microscopy');
if(useBGSubtractedResults)
    prResultDir=params.microscopy.bgSubtractedPRDir;
else
    prResultDir=params.microscopy.rawPRDir
end
numberOfMarkerSets=4;
modelGroupsPerMarkerSet=10;
modelsPerModelGroup=10;
numberOfSamples=36;
params=GetParams('samples');
modelNumbers=[params.samples.info.modelNumber];
tumorNumbers=[params.samples.info.repTumorNum];
regionNumbers=[params.samples.info.regionNumber];
%% Load Models
% 
% prModels=cell(numberOfMarkerSets,1);
% for markerSet=1:numberOfMarkerSets
%     prModels{markerSet}=cell(modelGroupsPerMarkerSet*modelsPerModelGroup,1);
%     for modelGroup=1:modelGroupsPerMarkerSet
%         modelFile=fullfile(prResultDir, ['PR_Models-MS' num2str(markerSet) '-R' ...
%             num2str(modelGroup) '.mat']);
%         models=load(modelFile);
%         for modelCounter=1:modelsPerModelGroup
%             modelNumber=(modelGroup-1)*modelsPerModelGroup+modelCounter;
%             prModels{markerSet}{modelNumber}=models.prModels{modelCounter};
%         end
%     end
% end

%% Load Results

subRegionPos=cell(numberOfMarkerSets,1);
sampleNumbers=cell(numberOfMarkerSets,1);
replicateNumbers=cell(numberOfMarkerSets,1);
sectionNumbers=cell(numberOfMarkerSets,1);
nFGBlocks=cell(numberOfMarkerSets,1);

prProfiles=cell(numberOfMarkerSets,1);
parfor markerSet=1:numberOfMarkerSets
    %prModels{markerSet}=cell(modelGroupsPerMarkerSet*modelsPerModelGroup,1);
    
    
    subRegionPos{markerSet}=cell(numberOfSamples,1);
    sampleNumbers{markerSet}=cell(numberOfSamples,1);
    replicateNumbers{markerSet}=cell(numberOfSamples,1); % e.g. 1,2,3
    sectionNumbers{markerSet}=cell(numberOfSamples,1); % this contains the physical section number e.g. 8,18,28
    nFGBlocks{markerSet}=cell(numberOfSamples,1);
    
    prProfiles{markerSet}=cell(modelGroupsPerMarkerSet*modelsPerModelGroup,1);
    prProfiles{markerSet}(:)={cell(numberOfSamples,1)};
    
    for modelGroup=1:modelGroupsPerMarkerSet
        
        
        for sampleNumber=1:numberOfSamples
            
            
            profilesFile=fullfile(prResultDir,['PR_Results_MS' num2str(markerSet) ...
                '_Sample' num2str(sampleNumber)...
                '_Group' num2str(modelGroup) '.mat']);
            if(exist(profilesFile,'file'))
                profilesData=load(profilesFile);
                
                subRegionPos{markerSet}{sampleNumber}=profilesData.subRegionPos;
                sampleNumbers{markerSet}{sampleNumber}=profilesData.sampleNumbers;
                replicateNumbers{markerSet}{sampleNumber}=profilesData.replicateNumbers;
                sectionNumbers{markerSet}{sampleNumber}=profilesData.sectionNumbers;
                nFGBlocks{markerSet}{sampleNumber}=profilesData.nFGBlocks{1};
                
                for modelCounter=1:modelsPerModelGroup
                    modelNumber=(modelGroup-1)*modelsPerModelGroup+modelCounter;
                    
                    prProfiles{markerSet}{modelNumber}{sampleNumber}=...
                        profilesData.superblockProfiles{modelCounter};
                end
            end
        end
    end
end

%%

isProfileCalculated=zeros(numberOfMarkerSets,...
    modelGroupsPerMarkerSet*modelsPerModelGroup,numberOfSamples);
for markerSet=1:numberOfMarkerSets
     for modelNumber=1:modelGroupsPerMarkerSet*modelsPerModelGroup
         isProfileCalculated(markerSet,modelNumber,:)=...
             ~cellfun(@isempty,prProfiles{markerSet}{modelNumber});
     end
end
%nMonnz(all(sum(isProfileCalculated,3)==numberOfSamples,1))

%%

tRegionPos=cell(numberOfMarkerSets,1);
tRegionInfo=struct;
combinedProfiles=cell(numberOfMarkerSets,1);
nFGBlocksByMarkerSet=cell(numberOfMarkerSets,1);
replicateSectionNumbers=cell(numberOfMarkerSets,1);
for markerSet=1:numberOfMarkerSets
    nFGBlocksByMarkerSet{markerSet}=vertcat(nFGBlocks{markerSet}{:});
    replicateSectionNumbers{markerSet}=vertcat(replicateNumbers{markerSet}{:});
    tRegionInfo(markerSet).sample=vertcat(sampleNumbers{markerSet}{:});
    tRegionInfo(markerSet).model=modelNumbers(tRegionInfo(markerSet).sample)';
    tRegionInfo(markerSet).tumor=tumorNumbers(tRegionInfo(markerSet).sample)';
    tRegionInfo(markerSet).region=regionNumbers(tRegionInfo(markerSet).sample)';
    [~,~,tRegionInfo(markerSet).section]=unique(...
        [tRegionInfo(markerSet).sample,replicateSectionNumbers{markerSet}],'rows');
    tRegionPos{markerSet}=vertcat(subRegionPos{markerSet}{:});
    completedModels=find(sum(isProfileCalculated(markerSet,:,:),3)==36);
    combinedProfiles{markerSet}=cell(length(completedModels),1);
    for modelCounter=1:length(completedModels)
        combinedProfiles{markerSet}{modelCounter}=vertcat(...
              prProfiles{markerSet}{completedModels(modelCounter)}{:});
    end
end
