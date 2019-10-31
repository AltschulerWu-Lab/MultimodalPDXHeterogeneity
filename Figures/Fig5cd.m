% Code for Figs 5c/d: Effect of different sampling strategies

addpath('../');
addpath(genpath('../Common/Plotting/'));
addpath('../Microscopy/PhenoRipper/')

% Data set fixed info
numberOfModels=4;
numberOfSamples=36;
numberOfSuperblocks=30;
numberOfReplicates=3;
params=GetParams();
numberOfMarkerSets=4;

% A model is characterized by its PhenoRipper profile calculated by
% analyzing combining results from all sub-images from the model.
% The goal is to sample sub-images in different ways from a model, and study how close
% close we get to the PhenoRipper profile for the model as a function of
% number of sub-images used

% Important params
scoreCutoff=0.01; % Max difference in PhenoRipper profile relative to true mode
% to be considered a good sampling of the model


%% Load profiles and info across all PhenoRipper runs.
% Note the PhenoRipper profiles are pre-calculated across multiple random
% initial seeds, and these results are being loaded here.
useBGSubtractedResults=true; % Set to false to generate non-bg subtracted images for supplement
[combinedProfiles,tRegionInfo,nFGBlocksByMarkerSet,replicateSectionNumbers]=...
    Load_PR_Results(useBGSubtractedResults);
%tRegionInfo contains the model/region/sample/section number for each
%sub-image

%% Calculate difference in PhenoRipper profiles etc
% Here we randomly sample a subset of areas based on different sampling strategies
% and calculate how different the PhenoRipper profile of the subset is from
% the true model profile.

% Params for #Samples calculations
minNFG=100;% Min Number of Foreground Blocks to USe image
subsetSizeRange=[1,5,20,30,40,50,75,100,125,150,175,200]; % range of sub-images used
numberOfRandomizations=100; % Number of random runs used to calculate results
differenceMeasure='AvgDiff'; % Measure used to compare profiles uses average difference in each component between true and sample


diffVals=cell(numberOfModels,numberOfMarkerSets);
prRunNumber=1;
for modelNumber=1:numberOfModels
    for markerSet=1:numberOfMarkerSets
        
        tic;
        numberOfPRRuns=length(combinedProfiles{markerSet});
        
        % We calculate difference between the reference profile and
        % subsamples for each run
        diffValsCombos=cell(numberOfPRRuns,1);
        
        nFG=nFGBlocksByMarkerSet{markerSet};
        fgSamplesInModel=(tRegionInfo(markerSet).model==modelNumber)&nFG>minNFG;
        
        
        markerSetProfiles= combinedProfiles{markerSet};
        parfor prRunNumber=1:numberOfPRRuns
            
            % nFg weighted wighted ref profile: this is the target
            refProfile=sum(bsxfun(@times,...
                markerSetProfiles{prRunNumber}(fgSamplesInModel,:),...
                nFG(fgSamplesInModel)/sum(nFG(fgSamplesInModel))),1);
            
            % define sampling superset: we will subsample from these
            % profiles
            allFgProfiles=markerSetProfiles{prRunNumber}(fgSamplesInModel,:);
            
            % Profile difference as a function of sample size
            groupingNames=fieldnames(tRegionInfo(markerSet));
            numberOfGroupings=length(groupingNames);
            
            diffValsCombos{prRunNumber}=zeros(numberOfGroupings,...
                length(subsetSizeRange),numberOfRandomizations);
            
            temp=zeros(numberOfGroupings,...
                length(subsetSizeRange),numberOfRandomizations);
            
            % Outer loop is the variable (e.g. tumor) value that all
            % subsamples will share
            for groupingCounter=1:numberOfGroupings
                field=groupingNames{groupingCounter};
                for subsetSizeCounter=1:length(subsetSizeRange)
                    subsetSize=subsetSizeRange(subsetSizeCounter);
                    profilesSuperSet=allFgProfiles;
                    
                    % This will generate a set of subsetSizeCounter profiles 
                    % from sub-images sharing the same value of the groupingField
                    subsetProfiles=GenerateRandomSubsetProfiles(profilesSuperSet,nFG(fgSamplesInModel),...
                        subsetSize,numberOfRandomizations,...
                        tRegionInfo(markerSet).(field)(fgSamplesInModel));
                    
                    % Calculate the difference between the subset profiles
                    % and the reference profile
                    diffValsCombos{prRunNumber}(groupingCounter,subsetSizeCounter,:)...
                        =Profile_Difference(refProfile,subsetProfiles,differenceMeasure);
                    
                end
            end
        end
        diffVals{modelNumber,markerSet}=diffValsCombos;
        toc;
    end
end
%% Engine: This performs the confidence calculations
% The basic idea is to get the distribution of differences from the true
% model profile using various sampling strategies x number of samples
% The confidence is the fraction of instances in which we get close enough
% (as defined by scoreCutoff above) to the true model profile

groupingNames=fieldnames(tRegionInfo);

%scores refer to max (saturation) confidence 
combinedScores=cell(numberOfModels,numberOfMarkerSets); % stores max confidence for each run
meanConfidence=cell(numberOfModels,numberOfMarkerSets); %stores mean (across random runs) of confidence for a sampling strategy as function of # samples
stdConfidence=cell(numberOfModels,numberOfMarkerSets);
meanScores=zeros(numberOfModels,numberOfMarkerSets,length(groupingNames));%average (across runs) of max confidence
stdScores=zeros(numberOfModels,numberOfMarkerSets,length(groupingNames));

%counter=1;
for modelNumber=1:numberOfModels
    for markerSet=1:numberOfMarkerSets
        
        numberOfPRRuns=length(diffVals{modelNumber,markerSet});
        combinedScores{modelNumber,markerSet}=zeros(length(groupingNames),numberOfPRRuns);
        allConfidences=zeros(length(groupingNames),length(subsetSizeRange),numberOfPRRuns);
        
        for prRunNumber=1:numberOfPRRuns
            confidence=100*sum(diffVals{modelNumber,markerSet}{prRunNumber}<scoreCutoff,3)/numberOfRandomizations;
            confidence(sum(~isnan(diffVals{modelNumber,markerSet}{prRunNumber}),3)~=numberOfRandomizations)=NaN; %only preserve runs with full data
            allConfidences(:,:,prRunNumber)=confidence;
            
            combinedScores{modelNumber,markerSet}(:,prRunNumber)=nanmax(confidence,[],2);

        end
        
        % Calculate averages/std across all runs
        meanScores(modelNumber,markerSet,:)=mean(combinedScores{modelNumber,markerSet},2);
        stdScores(modelNumber,markerSet,:)=std(combinedScores{modelNumber,markerSet},[],2);
        meanConfidence{modelNumber,markerSet}=mean(allConfidences,3);
        stdConfidence{modelNumber,markerSet}=std(allConfidences,[],3);
            
            
        
    end
end



%% Fig5c: Example curves for a couple of models/marker set
cmap=[0.5,0.5,0.5;0.8*colormap(brewermap(length(groupingNames)-1,'RdYlBu'))];
counter=1;
groupingOrder=[2,4,3,1,5];

markerSet=1;
for modelCounter=1:2
    modelNumber=params.figs.modelOrder(modelCounter);

    figure;
    for gCounter=1:length(groupingNames)
        gN=groupingOrder(gCounter);
        plot(subsetSizeRange,meanConfidence{modelNumber,markerSet}(gN,:),'LineWidth',3,...
            'Color',cmap(gCounter,:));
        hold on;
    end
    set(gca,'YLim',[0,100],'YTick',[0,50,100],'XTick',[0,100,200]);
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    
end

%% Fig5d: Barplot Summaries averaging across all models and marker sets
groupConfidenceScores=[];
groupNum=[];
colorOrder=[2,5,4,3,1];
for i=1:length(groupingNames)
    gN=groupingOrder(i);
    temp=meanScores(params.figs.modelOrder,params.microscopy.markerSetsToUse,gN);
    %temp(temp==0)=[];
    groupConfidenceScores=[groupConfidenceScores; temp(:)];
    groupNum=[groupNum;i*ones(numel(temp),1)];
end

figure;
boxplot(groupConfidenceScores,groupNum);
h=boxplot(groupConfidenceScores,groupNum,...
    'PlotStyle','traditional','Colors',cmap);
set(h,{'linew'},{2})
set(gca,'YTickLabel',[],'XTickLabel',[],'YTick',[0,50,100]);

