function params=GetParams(dataType)

% IMPORTANT:
% CHANGE THE FOLLOWING LINES TO CORRESPOND TO YOUR DATA LOCATIONS
rootDataDir='/work/bioinformatics/srajar/Leidos/Data/Deposited/';
rootInfoDir='/work/bioinformatics/srajar/Leidos/Data/Info_Final/';
figSaveDir='/home2/srajar/Work/Leidos/Paper/Figs/';
% INPUTS
% dataType      -	String indictaing which parameters we want. It can take one of the following values: 'all','samples','RPPA','growth','BI','pathology','RNA',
% 			 'microscopy','DNA','figs','pathway'. 
% OUTPUT	
% params        -	Is a struct with fields corresponding to different types of parameters, with the parameter values appearing as subfields within a field. 
% 			If the dataType is 'all' then all possible fields and their corresponding parameters are returned. Otherwise only the field specified (and its corresponding parameters) are returned.


validDataTypes={'all','samples','RPPA','growth','BI','pathology','RNA',...
    'microscopy','DNA','figs','pathway'};
if(nargin<1)
    dataType='all';
else
    
    if(~all(ismember(dataType,validDataTypes)))
        error(strjoin(['Invalid dataType. Must be one of:', validDataTypes],' '));
    end
end


params=struct;
warning('off','MATLAB:table:ModifiedVarnames');
%% Sample Info
if(any(ismember(dataType,{'all','samples'})))
    sampleInfoFile=fullfile(rootInfoDir,'Sample_Info.csv');
    info=readtable(sampleInfoFile);
    
    numberOfSamples=height(info);
    [modelNumbers,modelNames]=grp2idx(info.Model_);
    [regionNumbers,regionLetters]=grp2idx(info.Section);
    
    regionLetter2Name=containers.Map({'D','C','V'},{'Dorsal','Central','Ventral'});
    
    repTumorNumber=zeros(numberOfSamples,1);
    for modelCounter=1:length(modelNames)
        idx=modelNumbers==modelCounter;
        repNum=grp2idx(info.Mouse_(idx));
        repTumorNumber(idx)=repNum;
    end
    
    sampleInfo=struct();
    for sampleCounter=1:numberOfSamples
        sampleNum=info.SampleID(sampleCounter);
        sampleInfo(sampleNum).modelName=info.Model_{sampleCounter};
        sampleInfo(sampleNum).modelNumber=modelNumbers(sampleCounter);
        sampleInfo(sampleNum).ptcID=info.Study_{sampleCounter};
        sampleInfo(sampleNum).passageNum=info.Passage_(sampleCounter);
        sampleInfo(sampleNum).mouseNum=info.Mouse_(sampleCounter);
        sampleInfo(sampleNum).repTumorNum=repTumorNumber(sampleCounter);
        sampleInfo(sampleNum).flank=info.flank{sampleCounter};
        sampleInfo(sampleNum).regionName=regionLetter2Name(info.Section{sampleCounter});
        sampleInfo(sampleNum).regionNumber=regionNumbers(sampleCounter);
    end
    
    
    params.samples.modelNames=modelNames;
    params.samples.modelColonSide={'L','R','R','L'};
    params.samples.regionNames=cellfun(@(x) regionLetter2Name(x),...
        regionLetters,'Unif',false);
    params.samples.info=sampleInfo;
    params.samples.numberOfSamples=numberOfSamples;
    
    numberOfModels=length(params.samples.modelNames);
    numberOfRegions=length(params.samples.regionNames);
    params.samples.ReplicateTumorsPerModel=3;
    sampleMap=zeros(numberOfRegions,params.samples.ReplicateTumorsPerModel,numberOfModels);
    for regionCounter=1:numberOfRegions
        for replicateCounter=1:params.samples.ReplicateTumorsPerModel
            for modelCounter=1:numberOfModels
                sn=find([sampleInfo.regionNumber]==regionCounter & ...
                    [sampleInfo.repTumorNum]==replicateCounter & ...
                    [sampleInfo.modelNumber]==modelCounter);
                if(length(sn)~=1)
                    error('Oops something went wrong');
                else
                    sampleMap(regionCounter,replicateCounter,modelCounter)=sn;
                end
            end
        end
    end
    params.samples.sampleMap=sampleMap;
    
    
    params.samples.pairs.all=nchoosek(1:numberOfSamples,2);
    
    pairsPerModel= nchoosek(numberOfSamples/numberOfModels,2);
    params.samples.pairs.sameModel=zeros(numberOfModels*...
        pairsPerModel,2);
    for modelCounter=1:numberOfModels
        pairRange=(1:pairsPerModel)+(modelCounter-1)*pairsPerModel;
        params.samples.pairs.sameModel(pairRange,:)=nchoosek(...
            find([sampleInfo.modelNumber]==modelCounter),2);
    end
    
    pairsPerRegion= nchoosek(numberOfSamples/numberOfRegions,2);
    params.samples.pairs.sameRegion=zeros(numberOfRegions*...
        pairsPerRegion,2);
    for regionCounter=1:numberOfRegions
        pairRange=(1:pairsPerRegion)+(regionCounter-1)*pairsPerRegion;
        params.samples.pairs.sameRegion(pairRange,:)=nchoosek(...
            find([sampleInfo.regionNumber]==regionCounter),2);
    end
    
    sideNames={'L','R'};
    numberOfSides=2;%Left+Right
    pairsPerSide= nchoosek(numberOfSamples/numberOfSides,2);
    params.samples.pairs.sameSide=zeros(numberOfSides*...
        pairsPerSide,2);
    for sideCounter=1:numberOfSides
        pairRange=(1:pairsPerSide)+(sideCounter-1)*pairsPerSide;
        params.samples.pairs.sameSide(pairRange,:)=nchoosek(...
            find(strcmp(params.samples.modelColonSide([sampleInfo.modelNumber]),...
            sideNames{sideCounter})),2);
    end
    
    [tumorNumbers,uniqueMouseNumbers]=grp2idx([sampleInfo.mouseNum]);
    numberOfTumors=length(uniqueMouseNumbers);
    pairsPerTumor= nchoosek(numberOfSamples/numberOfTumors,2);
    params.samples.pairs.sameTumor=zeros(numberOfTumors*...
        pairsPerTumor,2);
    for tumorCounter=1:numberOfTumors
        pairRange=(1:pairsPerTumor)+(tumorCounter-1)*pairsPerTumor;
        params.samples.pairs.sameTumor(pairRange,:)=nchoosek(...
            find(tumorNumbers==tumorCounter),2);
    end
    
    [mrNumbers,uniqueMRs]=grp2idx(10*[sampleInfo.modelNumber]+...
        [sampleInfo.regionNumber]);
    numberOfMRs=length(uniqueMRs);
    pairsPerMR= nchoosek(numberOfSamples/numberOfMRs,2);
    params.samples.pairs.sameModelAndRegion=zeros(numberOfMRs*...
        pairsPerMR,2);
    for mrCounter=1:numberOfMRs
        pairRange=(1:pairsPerMR)+(mrCounter-1)*pairsPerMR;
        params.samples.pairs.sameModelAndRegion(pairRange,:)=nchoosek(...
            find(mrNumbers==mrCounter),2);
    end
    %params.samples.pairs.sameModelAndRegion
    
    
end
%% RPPA Info
if(any(ismember(dataType,{'all','RPPA'})))
    params.rppa.dataFile=fullfile(rootDataDir,...
        'RPPA/RPPA_Processed_Results.txt');
    params.rppa.antibodyFile=fullfile(rootDataDir,...
        'RPPA/RPPA_Antibody_Info.txt');
end

%% RNA Info
if(any(ismember(dataType,{'all','RNA'})))
    params.rna.combinedRPMFile=fullfile(rootInfoDir,'normalized_reads.xls');
    %params.rna.bedFile=fullfile(rootDataDir,'Reference',...
    %    'hg19_AmpliSeq_Transcriptome_21K_v1_Primer_to_visualize_in_UCSC_Browser.bed');
    params.rna.bedFile=fullfile(rootInfoDir,'hg19_AmpliSeq_Transcriptome_21K_v1.bed');
    params.rna.numberOfGenes=20812;
    params.rna.numberOfChips=4;
    params.rna.rootDir=rootDataDir;
    params.rna.badSamples=[7,8,10,12,24,37]; 
    
    % Note, this file corresponds to page 7 of deposited file 'RNA_Combined/AmpliSeq Transcriptome Performance Summary.pdf' 
    params.rna.qualityFile=fullfile(rootInfoDir,'RNA_Quality.txt'); 
end

%% Growth Curves
if(any(ismember(dataType,{'all','growth'})))
    params.growth.dataFolder=fullfile(rootDataDir,'Growth_Curves');
    %params.growth.dataFolder=fullfile(maikesAnnoDir,'GrowthCurves');
    
    %params.growth.model2File=containers.Map({'CN1571','CN1574','CN1572','CR-0104-O'},...
    %    fullfile(params.growth.dataFolder,{'PTC1943_71416.xlsx',...
    %    'PTC1925_72116.xlsx','PTC1966_71416.xlsx','PTC1951_80416.xlsx'}));

    params.growth.model2File=containers.Map({'CN1571','CN1574','CN1572','CR-0104-O'},...
        fullfile(params.growth.dataFolder,{'PTC1943_Leidos.xlsx',...
        'PTC1925_Leidos.xlsx','PTC1966_Leidos.xlsx','PTC1951_Leidos.xlsx'}));
    
end
%% Bioinformatics
if(any(ismember(dataType,{'all','BI'})))
    params.BI.pathwayDataDir=rootInfoDir;
    params.BI.antiBodyFile=fullfile(rootDataDir,'Reference',...
        'Antibody_Info.txt');
end
%% Pathology
if(any(ismember(dataType,{'all','pathology'})))
    %params.pathology.dataDir=maikesAnnoDir;
    %params.pathology.maikeHEScoringFile=fullfile(params.pathology.dataDir,'H&E_annotation_necrosis.xlsx');
    params.pathology.scottHEScoringFile=fullfile(rootDataDir,'Pathology.txt');
end
%% Microscopy
if(any(ismember(dataType,{'all','microscopy'})))
    params.microscopy.annoDir=rootInfoDir;
    %params.microscopy.aperioImgDir='/home/aperio/';
    params.microscopy.ToDepositDir=rootDataDir;
    params.microscopy.heLookupFile=fullfile(params.microscopy.annoDir,'HE_lookup.xlsx');
    params.microscopy.ifLookupFile=fullfile(params.microscopy.annoDir,'IF_Imaging_lookup.xlsx');
    params.microscopy.markerSets={{'DAPI','Beta-Catenin','pS6','Vimentin'},...
        {'DAPI','Ki67','pErk','Vimentin'},...
        {'DAPI','E-Cadherin','pStat3','Vimentin'},...
        {'DAPI','pAkt','pVEGFR','Vimentin'}};
    params.microscopy.markerSetsAlt={{'dapi_if','betacatenin_if',...
        'p_rps6_if','vimentin_if'},{'dapi_if','ki67_if','p_erk_if','vimentin_if'},...
        {'dapi_if','ecadherin_if','p_stat3_if','vimentin_if'},...
        {'dapi_if','p_akt_if','p_vegfr_if','vimentin_if'}};
    params.microscopy.isaTabAssayFile=fullfile(params.microscopy.ToDepositDir,...
        'a_casix_Microscopy.txt');
    params.microscopy.markerSetsToUse=1:3;
    params.microscopy.processedResultsFileRaw=fullfile(rootInfoDir,...
        'PR_Processed_Results_FullRes.mat');
     params.microscopy.processedResultsFile=fullfile(rootInfoDir,...
        'PR_Result_Clean.mat');
    params.microscopy.avgValueFile=fullfile(rootInfoDir,...
        'markerAvgIntensitiesBGSub.mat');
    params.microscopy.avgValueRawFile=fullfile(rootInfoDir,...
        'markerAvgIntensitiesRaw.mat');
    params.microscopy.bgSubtractedImgList=fullfile(rootInfoDir,'Images',...
        'scaledImages.mat');
    params.microscopy.rawImgList=fullfile(rootInfoDir,'Images',...
        'rawImages.mat');
    params.microscopy.referenceDapiExposure=125; % All DAPI intensities are rescaled to a level corresponding to this exposure time
    params.microscopy.scalingPctiles=[25,75]; % The intensity percentile distributions used to normalize across replicate sections
    params.microscopy.bgSubtractedPRDir=fullfile(rootInfoDir,'PR_Results/');
    params.microscopy.rawPRDir='/work/bioinformatics/srajar/Leidos/Data/PR_Results_Raw/';
    params.microscopy.downsamplingResultsDir=fullfile(rootInfoDir,'Spatial_Downsampling/');
    %params.microscopy.prResultDir='/work/bioinformatics/srajar/Leidos/Data/PR_Results/';
    %params.microscopy.prResultDirRaw='/work/bioinformatics/srajar/Leidos/Data/PR_Results_Raw/';
    % PhenoRipper analysis
    params.microscopy.prBlockSize=20;
    params.microscopy.prCutoffIntensity=10;
    params.microscopy.prSubImgDims=[1300,1300];
    params.microscopy.prNumberOfTrainingImages=100;
    params.microscopy.prBgSubMarkerScales=repmat([0,30000],4,1);
    params.microscopy.prRawMarkerScales=repmat([5000,30000],4,1);
    
end
%warning('on','MATLAB:table:ModifiedVarnames');
%% DNA Info
if(any(ismember(dataType,{'all','DNA'})))
    params.dna.processedResultsFile=fullfile(rootInfoDir,'DNA_Processed_Results_Final.mat');
    params.dna.badSamples=[15,17,25,33];
    params.dna.annovarResultDir=fullfile(rootInfoDir,'annovar_results');

end

%% Plotting options
if(any(ismember(dataType,{'all','figs'})))
    params.figs.modelOrder=[1,4,3,2]; %'CN1571','CR-0104-O', 'CN1574','CN1572'. First left sided then right.
    %params.figs.modelColors=brewermap(numberOfModels,'Set2');
    params.figs.modelColors=brewermap(numberOfModels,'PrGn');
    params.figs.saveDir=figSaveDir;
    params.figs.cmsColors=[232,158,51;0 115 172;208,121,164; 0 158 118]/255;
    params.figs.repTumorOrder=[1,3,2]; % for the hamburger plots swap the second and 3rd column

end

%% Pathways
if(any(ismember(dataType,{'all','pathway'})))
   params.pathway.gsvaResultsDir=fullfile(rootInfoDir,'gsva_results');
   params.pathway.gsvaFile=fullfile(params.pathway.gsvaResultsDir,...
       'gsvaHallmark_rnaseq.txt');
end



end

