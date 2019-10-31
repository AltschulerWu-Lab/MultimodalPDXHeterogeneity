function data=Load_Data(dataSet,samplesToGet,varargin)
% Main function to load up different types of data, and arrange them under a uniform data organization scheme.
% INPUTS
% Required:
% dataSet     	      -	Cell array listing the the types of data we want returned. Allowed values are 'RPPA','RNA','DNA','microscopy','pathway','all'. 
% 		    	The value 'all' is equivalent to listing all of the other entries.
% samplesToGet        -	Can be either 'tumor' if we want results for the 36 samples, or 'all' if we want to get the spike-in and other samples additionally.
% Optional:
% dropBadReadouts     - Boolean with default value true. This is only applied to the RPPA data, but if true will lead to antibodies failing QC to be dropped.
% 			when reporting expression for all samples.
% normalizeExpression - Boolean with default value true.  If true normalization for each feature is applied (for example gene expression is z-score normalized
% 			for each gene across samples).
%
% OUTPUT
% data                - Is a struct with fields for different kinds of information/assays (e.g. dna, rna etc to be described below). Each field is itself a struct
% 			containing with a fixed set of sub-fields as follows. 
% 			1) data: this sub-field is a 2D matrix with dimensions number_of_samples x number_of_readouts, i.e. rows describe the different samples,
% 			while columns represent different measurements/properties of the samples (the number/type of measurements is assay specific)
% 			2) columnLabels: A cell array of length number_of_readouts, describing the columns of the data matrix
% 			3) flaggedSamples: This optional field is a boolean vector of length number_of_samples, with true values indicating quality control issues 
% 			were detected in the measurement of the corresponding sample (i.e. row of data)
% 		      - The fields of data (structured as described above with data/columnLabels/flaggedSamples subfields) are
% 		        1) model: This field is always present (regardless of the dataSet parameter) and tells us which model each sample came from. 
% 		        Each of the 4 columns in the data sub-field represents a specific model, and the data matrix is true if the sample corresponding to that row, came from the 
% 		        corresponding model. The model field does not contain a flaggedSamples sub-field.
% 		        2) region: This field is always present (regardless of the dataSet parameter) and tells us which tumor region (dorsal/central/ventral)
% 		        each sample came from. Each of the 3 columns in the data sub-field represents a specific region, and the data matrix is true if the sample corresponding to that row, came from the 
% 		        corresponding region. The region field does not contain a flaggedSamples sub-field.
% 		        3) tumor: This field is always present (regardless of the dataSet parameter) and tells us which of the 3 replicate tumors for a PDX model
% 		        each sample came from. Each of the 3 columns in the data sub-field represents a specific replicate tumor, and the data matrix is true if the sample 
% 		        corresponding to that row, came from the corresponding tumor. The tumor field does not contain a flaggedSamples sub-field.
% 		        4) growth: This field is always present (regardless of the dataSet parameter) and tells us about the growth of the tumor from which the sample was
% 		        extracted. Here, the data-subfield has two columns, the first denoting the doubling time in days, and the second the final tumor volume in cubic mm.
% 		        Note: that if the normalizeExpression field is set to true, these values will be z-score normalized. The growth field does not contain a flaggedSamples sub-field.
% 		        5) pathology: This field is always present (regardless of the dataSet parameter) and tells us about the histopathological compostion of the sample as
% 		        assessed by a pathologist. Here, the data-subfield has three columns, corresponding to the percentage of the slide that was tumor, stroma and necrotic respectively.
% 		        Note: that if the normalizeExpression field is set to true, these values will be z-score normalized. The pathology field does not contain a flaggedSamples sub-field.
% 		        5) microscopy: This field is only present dataSet parameter contains 'microscopy' or 'all' and tells us about PhenoRipper composition of the sample across 4 
% 		        immunofluorescence data sets.The data-subfield has 120 columns, corresponding to the 30 components corresponding to the PhenoRipper profile for a single marker set
% 		        , which are then concatenated across the 4 marker sets.Note: that if the normalizeExpression field is set to true, these values will be z-score normalized. 
% 		        The microscopy field contains a flaggedSamples sub-field.
% 		        6) rppa: This field is only present dataSet parameter contains 'rppa' or 'all' and tells us about protein-level measurement for the sample across different protein
% 		        antibodies.The data-subfield has 45 columns (2 more ifdropBadReadout is set to false), corresponding to the 45 antibodies measured for each sample.
% 		        Note: that if the normalizeExpression field is set to true, these values will be z-score normalized. The rppa field contains a flaggedSamples sub-field.
% 		        7) dna: This field is only present dataSet parameter contains 'dna' or 'all' and tells us about the presence or absence of various SNV events. 
% 		        The data-subfield is a binary matrix with 44 columns corresponding to the 44 events that were present at least one sample (true value indicates presence of event).
% 		        The identifty of the event can be determined from the columnLabels sub-field. Note: that if the normalizeExpression field is set to true, these values will be z-score normalized. 
% 		        The dna field contains a flaggedSamples sub-field.
% 		        8) rna: This field is only present dataSet parameter contains 'rna' or 'all' and reports on the FPKM normalized transcriptomic profiling of the samples. 
% 		        The data-subfield is contains 20812 columns corresponding to the genes profiles, and the count data is FPKM normalized. Note: that if the normalizeExpression field is set to true, 
% 		        these values will be z-score normalized: without z-score normalization the rows should add up to 1E6 and with normalization the columns should have mean zero. 
% 		        The rna field contains a flaggedSamples sub-field.
% 			
	
validDataSets={'RPPA','RNA','DNA','microscopy','pathway','all'};
validSamplesToGet={'all','tumor'};

p=inputParser;
addRequired(p,'dataSet',@(x) (ischar(x) && ismember(x,validDataSets))||...
    (iscell(x) && all(ismember(x,validDataSets))));
addRequired(p,'samplesToGet',@(x) ischar(x) & ismember(x,validSamplesToGet));
addOptional(p,'dropBadReadouts',true,@(x) islogical(x));
addOptional(p,'normalizeExpression',true,@(x) islogical(x));
parse(p,dataSet,samplesToGet,varargin{:});


data =struct;
params=GetParams();
sampleInfo=params.samples.info;
modelNames=params.samples.modelNames;
regionNames=params.samples.regionNames;
numberOfSamples=params.samples.numberOfSamples;
numberOfReplicateTumors=length(unique([sampleInfo.repTumorNum]));


%% Load Tumor Growth Data
[growthCurves,doublingTimes,~]=GetSampleGrowthData();
finalTumorVolumes=cellfun(@(x) x(end,2),growthCurves);
[pctTumor,pctStroma,pctNecrosis]=GetPathologyInfo();

tumorGrowthMat=[doublingTimes,finalTumorVolumes];
pathMat=[pctTumor,pctStroma,pctNecrosis];

isModelMat=false(numberOfSamples,length(modelNames));
isModelMat(sub2ind(size(isModelMat),1:numberOfSamples,...
    [sampleInfo.modelNumber]))=true;


isRegionMat=false(numberOfSamples,length(regionNames));
isRegionMat(sub2ind(size(isRegionMat),1:numberOfSamples,...
    [sampleInfo.regionNumber]))=true;


isRepMat=false(numberOfSamples,numberOfReplicateTumors);
isRepMat(sub2ind(size(isRepMat),1:numberOfSamples,...
    [sampleInfo.repTumorNum]))=true;



normalize=@(x) bsxfun(@rdivide,bsxfun(@minus,x,mean(x,1)),...
    std(x,[],1));
%%
%switch(dataSet)
if(any(ismember(dataSet,{'all','RPPA'})))
    %case 'RPPA'
    
    % Load RPPA Data
    %[nums,txt]=xlsread(params.rppa.dataFile,'Satwik final');
    temp=readtable(params.rppa.dataFile);
    temp.CD133=nan*ones(height(temp),1);
    sampleNames=temp.Y;
    nums=table2array(temp(:,2:end));
    tumorData=nums(1:numberOfSamples,:);
    
    %antibodyNames=txt(1,2:end);
    fileID=fopen(params.rppa.dataFile,'r');
    temp=textscan(fileID,'%s',1,'delimiter','\n');
    fclose(fileID);
    antibodyNames=regexp(temp{1}{1},'\t','split');
    antibodyNames=strtrim(antibodyNames(2:end));
    abInfo=readtable(params.rppa.antibodyFile);
    [~,idx]=ismember(antibodyNames,abInfo.Combined_Result_Column);
    
    %data.rppa.columnLabels=antibodyNames;
    data.rppa.columnLabels=abInfo.Antibody_ID(idx);
    
    switch(samplesToGet)
        case 'tumor'
            data.rppa.data=tumorData;
        case 'all'
            
            wholeTissueData=nums((1:numberOfSamples)+numberOfSamples,:);
            residualData=nums((1:numberOfSamples)+2*numberOfSamples,:);
            pureMouseData=nums((1:5)+3*numberOfSamples,:);
            data.rppa.data=[tumorData;wholeTissueData;residualData;pureMouseData];
            
            mouseSampleNames=sampleNames((1:5)+3*numberOfSamples,1);
            mouseSampleNumbers=cellfun(@str2double,regexp(mouseSampleNames,' \d* ','match'));
            
            isModelMatMouse=false(length(mouseSampleNumbers),length(modelNames));
            isModelMatMouse(sub2ind(size(isModelMatMouse),1:length(mouseSampleNumbers),...
                [sampleInfo(mouseSampleNumbers).modelNumber]))=true;
            
            isRegionMatMouse=false(length(mouseSampleNumbers),length(regionNames));
            isRegionMatMouse(sub2ind(size(isRegionMatMouse),1:length(mouseSampleNumbers),...
                [sampleInfo(mouseSampleNumbers).regionNumber]))=true;
            
            isRepMatMouse=false(length(mouseSampleNumbers),numberOfReplicateTumors);
            isRepMatMouse(sub2ind(size(isRepMatMouse),1:length(mouseSampleNumbers),...
                [sampleInfo(mouseSampleNumbers).repTumorNum]))=true;
            
            
            isModelMat=[repmat(isModelMat,3,1);isModelMatMouse];
            isRegionMat=[repmat(isRegionMat,3,1);isRegionMatMouse];
            isRepMat=[repmat(isRepMat,3,1);isRepMatMouse];
            pathMat=[repmat(pathMat,3,1);pathMat(mouseSampleNumbers,:)];
            tumorGrowthMat=[repmat(tumorGrowthMat,3,1);...
                tumorGrowthMat(mouseSampleNumbers,:)];
    end
    if(p.Results.dropBadReadouts)
        isBadReadout=any(isnan(data.rppa.data),1);
        data.rppa.data=data.rppa.data(:,~isBadReadout);
        data.rppa.columnLabels=data.rppa.columnLabels(~isBadReadout);
    end
    if(p.Results.normalizeExpression)
        data.rppa.data=normalize(data.rppa.data);
    end
    
    data.rppa.flaggedSamples=false(size(data.rppa.data,1),1);
end
%case 'RNA'
if(any(ismember(dataSet,{'all','RNA'})))
    showBed=false;
    rnaData=readtable(params.rna.combinedRPMFile);
    
    rnaVals=table2array(rnaData(:,3:end));
    geneNames=rnaData.Gene;
    sampleNames=rnaData.Properties.VariableNames(3:end);
    sampleNumbers=cellfun(@(x) str2double(x{1}),regexp(sampleNames,'\d\d\d$','match'));
    
    %sampleNumbers(sampleNumbers==37)=8;
    data.rna.data=zeros(size(rnaVals'));
    data.rna.data(sampleNumbers,:)=rnaVals';
    data.rna.flaggedSamples=false(size(rnaVals,2),1);
    badSampleNumbers=params.rna.badSamples;
    data.rna.flaggedSamples(badSampleNumbers)=true;
    data.rna.columnLabels=geneNames;
    
    if(showBed)
        data.rna.ampliconNames=rnaData.Target;
        
        fileId=fopen(params.rna.bedFile);
        textscan(fileId, '%s',1,'Delimiter','\n');
        bedData=textscan(fileId, '%s %d %d %s %d %s %s %s','Delimiter','\t');
        fclose(fileId);
        numberOfAmplicons=length(bedData{1});
        temp=regexp(bedData{8},';','split');
        bedInfo=struct;
        for ampliconCounter=1:numberOfAmplicons
            bedInfo(ampliconCounter).NM=bedData{1}{ampliconCounter};
            bedInfo(ampliconCounter).start=bedData{2}(ampliconCounter);
            bedInfo(ampliconCounter).stop=bedData{3}(ampliconCounter);
            bedInfo(ampliconCounter).Amplicon=bedData{4}{ampliconCounter};
            temp1=regexp(temp{ampliconCounter},'=','split');
            for i=1:length(temp1)
                bedInfo(ampliconCounter).(temp1{i}{1})=temp1{i}{2};
            end
        end
        [isInBedInfo,idx]=ismember(rnaData.Target,bedData{4});
        if(all(isInBedInfo))
            data.rna.ampliconInfo=bedInfo(idx);
        else
            error('Amplicon Mismatch');
        end
    end
    switch(samplesToGet)
        case 'tumor'
            data.rna.data=data.rna.data(1:36,:);
            data.rna.flaggedSamples=data.rna.flaggedSamples(1:36);
        case 'all'
            isModelMat=[isModelMat;isModelMat(8,:)];
            isRegionMat=[isRegionMat;isRegionMat(8,:)];
            isRepMat=[isRepMat;isRepMat(8,:)];
            tumorGrowthMat=[tumorGrowthMat;tumorGrowthMat(8,:)];
            pathMat=[pathMat;pathMat(8,:)];
    end
    if(p.Results.normalizeExpression)
        data.rna.data=normalize(data.rna.data);
    end
end
%case 'DNA'
if(any(ismember(dataSet,{'all','DNA'})))
    dnaData=load(params.dna.processedResultsFile,...
        'uniqueIDs','isMutated','eventInfo','altAlleleFraction','eventNames');
    [dnaSamplesToUse,dnaMutationsToUse]=FilterDnaMutations(dnaData,...
        'mutationTypesToDrop',{'none'},'keepOnlyCosmic',false,...
        'filterMouse',true);
    switch(samplesToGet)
        case 'tumor'
            data.dna.data=dnaData.isMutated(1:36,dnaMutationsToUse);
            %data.dna.flaggedSamples=data.rna.flaggedSamples(1:36);
        case 'all'
            data.dna.data=dnaData.isMutated(:,dnaMutationsToUse);
            %data.dna.flaggedSamples=data.rna.flaggedSamples(1:36);
    end
    data.dna.columnLabels=dnaData.eventNames(dnaMutationsToUse);%dnaData.uniqueIDs;
    badSamples=params.dna.badSamples;
    data.dna.flaggedSamples=false(size(data.dna.data,1),1);
    data.dna.flaggedSamples(badSamples)=true;
    if(p.Results.normalizeExpression)
        data.dna.data=normalize(double(data.dna.data));
    end
end
%case 'microscopy'
if(any(ismember(dataSet,{'all','microscopy'})))
    prData=load(params.microscopy.processedResultsFile,'sampleProfiles');
    profiles=prData.sampleProfiles';
    columnLabels={};
    for ms=1:length(profiles)
        prefix=['MS' num2str(ms) '_SB'];
        columnLabels=[columnLabels arrayfun(@(x) strcat(prefix,num2str(x)),1:30,'Unif',false)];
    end
    
    %numberSbs=sum(cellfun(@(x) size(x,2), prData.sampleProfiles)
    
    switch(samplesToGet)
        case 'tumor'
            data.microscopy.data=horzcat(profiles{:});
            %data.dna.flaggedSamples=data.rna.flaggedSamples(1:36);
        case 'all'
            data.microscopy.data=horzcat(profiles{:});
            %data.dna.flaggedSamples=data.rna.flaggedSamples(1:36);
    end
    if(p.Results.normalizeExpression)
        data.microscopy.data=normalize(data.microscopy.data);
    end
    data.microscopy.columnLabels=columnLabels';
    data.microscopy.flaggedSamples=false(size(data.microscopy.data,1),1);
end

if(any(ismember(dataSet,{'all','pathway'})))
    pathwayRes=readtable(params.pathway.gsvaFile);
    pathwayNames=pathwayRes.Var1;
    gsvaSampleNumbers=cellfun(@(x) str2double(x{1}),regexp(...
        pathwayRes.Properties.VariableNames(2:end),'\d\d\d$','match'));
    gsvaData=table2array(pathwayRes(:,2:end))';
    data.pathway.data=zeros(size(gsvaData));
    data.pathway.data(gsvaSampleNumbers,:)=gsvaData;
    switch(samplesToGet)
        case 'tumor'
            data.pathway.data=data.pathway.data(1:params.samples.numberOfSamples,:);
            data.pathway.flaggedSamples=false(params.samples.numberOfSamples,1);
            temp=params.rna.badSamples;
            temp=temp(temp<=params.samples.numberOfSamples);
            data.pathway.flaggedSamples(temp)=true;
        case 'all'
            data.pathway.data=data.pathway.data;
            data.pathway.flaggedSamples=false(size(data.pathway.data,1),1);
            data.pathway.flaggedSamples(params.rna.badSamples)=true;
    end
    if(p.Results.normalizeExpression)
        data.pathway.data=normalize(data.pathway.data);
    end
    data.pathway.columnLabels=pathwayNames;
  
    
end

%     otherwise
%         error('Invalid Data Set');
%
% end

data.model.data=isModelMat;
data.model.columnLabels=modelNames;

data.region.data=isRegionMat;
data.region.columnLabels=regionNames;

data.tumor.data=isRepMat;
data.tumor.columnLabels={'1','2','3'};


if(p.Results.normalizeExpression)
    data.growth.data=normalize(tumorGrowthMat);
else
    data.growth.data=tumorGrowthMat;
end
data.growth.columnLabels={'2X Time (days)','Tumor Volume (cubic mm)'};


if(p.Results.normalizeExpression)
    data.pathology.data=normalize(pathMat);
else
    data.pathology.data=pathMat;
end
data.pathology.columnLabels={'% Tumor','% Stroma', '% Necrosis'};

% Reorder fields so that RNA/RPPA are at the end
preferredOrder={'model','region','tumor','growth','pathology','microscopy',...
    'rppa','rna','pathway','dna'};
data=orderfields(data,preferredOrder(ismember(preferredOrder,fieldnames(data))));

end
