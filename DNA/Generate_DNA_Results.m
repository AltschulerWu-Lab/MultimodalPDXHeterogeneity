function Generate_DNA_Results(saveFile)

% Convenience function to pull data from individual vcf files and combine into a single data struct	
% INPUT
% saveFilename   -     This is the name of the mat file where the data struct will be saved
%              
% OUTPUT
%   The output is struct saved as a mat file with fields containing information on the dna calls. 
%   Important fields are: 
%    uniqueIDs   - a cell array containing a unique identifier for mutational events allowing them to be compared across vcf files.
%   isMutated   - A logical array with rows corresponding to the mutational events specified in uniqueIDs and columns corresponding to the different samples. An element is true if that event is present in that sample.
%   eventInfo      - a struct array describing the events in uniqueIDs
%   altAlleleFraction - a numeric matrix with same structure as isMutated, denoting fraction of alt reads
    params=GetParams();
    %%
    annoDir=params.dna.annovarResultDir;

    fileNames=dir(fullfile(annoDir,'*multianno.vcf'));
    fileNames={fileNames.name};
    sampleNumbers=str2double(cellfun(@(x) x{1},regexp(fileNames,'\d\d\d','match'),'Unif',false));
    fileNames=fullfile(annoDir,fileNames);
    numberOfSamples=length(sampleNumbers);
    %% Load Data
    vcfData=cell(length(sampleNumbers),1);
    eventIDs=cell(length(sampleNumbers),1);
    for sampleNumber=1:length(sampleNumbers)
         txtFile=fileNames{sampleNumbers==sampleNumber};
         tic;
         [vcfData{sampleNumber},eventIDs{sampleNumber}]=Parse_VCF(txtFile);
         toc;
    end
    %% Organize Data into Table
    sN=arrayfun(@(x) x*ones(1,length(eventIDs{x})), 1:length(eventIDs),'Unif',false);
    sN=horzcat(sN{:});
    rN=arrayfun(@(x) 1:length(eventIDs{x}), 1:length(eventIDs),'Unif',false);
    rN=horzcat(rN{:});
    allEventIDs=horzcat(eventIDs{:});
    nEvents=cellfun(@(x) [x.nEvents],vcfData,'Unif',false);
    nEvents=horzcat(nEvents{:});
    [idNum,uniqueIDs]=grp2idx(allEventIDs);
    numberOfUniqueEvents=length(uniqueIDs);

    % Mutation/Alt allele fraction
    isMutated=false(numberOfSamples,numberOfUniqueEvents);
    altAlleleFraction=zeros(numberOfSamples,numberOfUniqueEvents);
    %vcfRes=struct;

    qt=false(size(allEventIDs));
    gt=cell(size(allEventIDs));
    for callCounter=1:length(allEventIDs)
        temp=vcfData{sN(callCounter)}(rN(callCounter));
        isMutated(sN(callCounter),idNum(callCounter))=(temp.QT==1)&strcmp(temp.filter,'PASS');
        altAlleleFraction(sN(callCounter),idNum(callCounter))=temp.AF;
        qt(callCounter)=temp.QT==1;
        gt{callCounter}=temp.GT;
    end

    % Mutation Information
    eventInfo=struct;
    firstOccStruct=struct;
    nCountsList=zeros(numberOfUniqueEvents,1);
    for eventCounter=1:numberOfUniqueEvents
       firstOcc=find(idNum==eventCounter,1);
       nCountsList(eventCounter)=max(nEvents(idNum==eventCounter));
       info=vcfData{sN(firstOcc)}(rN(firstOcc));
       %firstOccStruct(eventCounter)=info;
       eventInfo(eventCounter).func=info.Func_refGene;
       eventInfo(eventCounter).geneName=info.Gene_refGene;
       eventInfo(eventCounter).mutType=info.ExonicFunc_refGene;
       temp=regexp(info.AAChange_refGene,',','split');
       temp=regexp(temp{1},':','split');
       eventInfo(eventCounter).aaChange=strcat(temp{1},':',temp{end});
       eventInfo(eventCounter).avSnp=info.avsnp147;
       eventInfo(eventCounter).cosmic=info.id;
       eventInfo(eventCounter).pos=info.pos;
       eventInfo(eventCounter).ref=info.ref;
       eventInfo(eventCounter).alt=info.alt;
       eventInfo(eventCounter).filter=info.filter;
       if(nCountsList(eventCounter)>1)
           eventInfo(eventCounter).asterisk='*';
       else
           eventInfo(eventCounter).asterisk='';
       end
    end
    eventNames=strcat(uniqueIDs,{eventInfo.asterisk}');
    isSynonymous=strcmp({eventInfo.mutType},'synonymous_SNV');
    isGood=strcmp({eventInfo.filter},'PASS');
    isNonSynSnv=strcmp({eventInfo.mutType},'nonsynonymous_SNV');
    isExonic=~strcmp({eventInfo.mutType},'exonic');
    isSNP=~cellfun(@isempty,regexp({eventInfo.avSnp},'rs','match'));
    isCOSMIC=~cellfun(@isempty,regexp({eventInfo.cosmic},'COSM','match'));
    
    save(saveFile);
end

