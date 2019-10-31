function [vcfResults,eventID,format,info,contigs]=Parse_VCF(vcfFilename)

% Convenience function to parse annoVar generated vcf files	
% INPUT
% vcfFilename   -     Name of the vcfFile to open as a character array.
%              
% OUTPUT
%   vcfResult   - a struct array containing the vcf results.i Each element of the array corresponds
%   		to a single line (i.e. event) of the vcf and the fields of the struct are different
%   		"columns" of the vcf
%   eventId     - A cell array containing an identifier of the event, allowing events to be compared
%   		across vcf files
%   format      - a struct array describing the format and contents of different fields in vcfResult
%   
%   info       - a struct array describing the format and contents of the "Info" fields of the vcf file
%   
%   contigs    - a struct array describing the length etc of the different chromosomes 

fileId=fopen(vcfFilename,'r');
vcfData=textscan(fileId,'%s','delimiter','\n');
fclose(fileId);
vcfData=vcfData{1};
%%
numberOfLines=length(vcfData);

info=struct;
contigs=struct;
format=struct;
vcfResults=struct;
dataCounter=1;
infoCounter=1;
formatCounter=1;
for lineCounter=1:numberOfLines
    
    if(~isempty(regexp(vcfData{lineCounter},'^##','once'))) % Description Line
        if(~isempty(regexp(vcfData{lineCounter},'^##INFO=<ID=.*>$','once'))) % Variable Definition
            tokens=regexp(vcfData{lineCounter},'ID=(\S*),Number=(\S*),Type=(\S*),Description=(.*)>','tokens','once');
            name=tokens{1};
            number=tokens{2};
            varType=tokens{3};
            description=tokens{4};
            if(isfield(info,name))
                warning(['Duplicate info fields with name ' name]);
            end
            info(infoCounter).name=name;
            info(infoCounter).number=number;
            info(infoCounter).type=varType;
            info(infoCounter).description=description;
            infoCounter=infoCounter+1;
        elseif(~isempty(regexp(vcfData{lineCounter},'^##contig=<ID=.*>$','once'))) % Variable Definition
            tokens=regexp(vcfData{lineCounter},'ID=(\w*),length=(\d*),assembly=(\w*)>','tokens','once');
            name=tokens{1};
            chrLength=tokens{2};
            assembly=tokens{3};
            if(isfield(contigs,name))
                warning(['Duplicate contig fields with name ' name]);
            end
            contigs.(name).length=chrLength;
            contigs.(name).assembly=assembly;
        elseif(~isempty(regexp(vcfData{lineCounter},'^##FORMAT=<ID=.*>$','once'))) % Variable Definition
            tokens=regexp(vcfData{lineCounter},'ID=(\S*),Number=(\S*),Type=(\S*),Description=(.*)>','tokens','once');
            name=tokens{1};
            number=tokens{2};
            varType=tokens{3};
            description=tokens{4};
            if(isfield(format,name))
                warning(['Duplicate format fields with name ' name]);
            end
            format(formatCounter).name=name;
            format(formatCounter).number=number;
            format(formatCounter).type=varType;
            format(formatCounter).description=description;
            formatCounter=formatCounter+1;
        end
        
        
    elseif(~isempty(regexp(vcfData{lineCounter},'^#','once'))) % Header Line
        headerFields=strsplit(regexprep(vcfData{lineCounter},'#',''),'\t');
    else % Data
        %dataTokens=strsplit(vcfData{lineCounter},'\t');
        dataTokens=regexp(vcfData{lineCounter},'\t','split');
        chrom=dataTokens{strcmp(headerFields,'CHROM')};
        pos=dataTokens{strcmp(headerFields,'POS')};
        ref=dataTokens{strcmp(headerFields,'REF')};
        qual=dataTokens{strcmp(headerFields,'QUAL')};
        filter=dataTokens{strcmp(headerFields,'FILTER')};
        
        
        idListTemp=regexp(dataTokens{strcmp(headerFields,'ID')},';','split');
        altList=regexp(dataTokens{strcmp(headerFields,'ALT')},',','split');
        
        numberOfEvents=length(altList);
        if(length(idListTemp)~=numberOfEvents)
            %warning('Mismatch in lengths of ID and ALT events');
            idList=cell(size(altList));
            if(length(idListTemp)==1)
              
               idList(:)={idListTemp{1}};
            else
              
                idList(:)={dataTokens{strcmp(headerFields,'ID')}};
          
            end
        else
            idList=idListTemp;
        end
        
        
        infoString=dataTokens{strcmp(headerFields,'INFO')};
        infoOrig=regexp(infoString,'(.*)ANNOVAR_DATE','tokens','once');
        infoAnnoVar=regexp(infoString,'ANNOVAR_DATE=\d\d\d\d-\d\d-\d\d;(\S*?);ALLELE_END','tokens');
        infoAnnoVar=cellfun(@(x) x{1},infoAnnoVar,'Unif',false);
        
        if(length(infoAnnoVar)~=numberOfEvents)
            error('Error in AnnoVar info extraction');
        end
        
        
        
        %         infoByAlle= regexp(dataTokens{strcmp(headerFields,'INFO')},'ALLELE_END','split');
        %         infoTokens=cellfun(@(x) regexp(x,'=','split'), ...
        %             regexp(dataTokens{strcmp(headerFields,'INFO')},';','split'),'Unif',false);
        %         infoNames=cellfun(@(x) x{1},infoTokens,'Unif',false);
        %         infoVals=cellfun(@(x) x{numel(x)},infoTokens,'Unif',false);
        
        if(dataCounter==1)
            formatOrder=regexp(dataTokens{strcmp(headerFields,'FORMAT')},':','split');
            nDataFields=length(formatOrder);
            fData=struct;
            for fieldCounter=1:nDataFields % This needs to be run just once!
                fName=formatOrder{fieldCounter};
                fData(fieldCounter).name=fName;
                fData(fieldCounter).niceName=matlab.lang.makeValidName(fName);
                fIdx=find(strcmp({format.name},fName));
                fData(fieldCounter).length=str2double(format(fIdx).number);
                fData(fieldCounter).type=format(fIdx).type;
            end
            
            
            tokens=cellfun(@(x) regexp(x,'=','split'), ...
                regexp( infoAnnoVar{1},';','split'),'Unif',false);
            annoFNames=cellfun(@(x) x{1},tokens,'Unif',false);
            
            nAnnoFields=length(annoFNames);
            annoData=struct;
            for fieldCounter=1:nAnnoFields % This needs to be run just once!
                fName=annoFNames{fieldCounter};
                annoData(fieldCounter).name=fName;
                annoData(fieldCounter).niceName=matlab.lang.makeValidName(fName);
                fIdx=find(strcmp({info.name},fName));
                %annoData(fieldCounter).length=str2double(format(fIdx).number);
                annoData(fieldCounter).type=info(fIdx).type;
            end
            
            
        end
        
        dataOrder=regexp(dataTokens{end},':','split'); % Can add support for multiple samples
        
        if(length(dataOrder)~=nDataFields)
            error('Format String Does not Match Sample Data');
        end
        
        
        
        
        dataCell=cell(nDataFields,1);
        for fieldCounter=1:nDataFields
            fLength=fData(fieldCounter).length;
            if(isnan(fLength))
                fLength=numberOfEvents;
            end
            fType=fData(fieldCounter).type;
            tokens=regexp(dataOrder{fieldCounter},',','split');
            if(length(tokens)~=fLength)
                error('Problem with length of data');
                
            else
                switch(fType) % switch type then generate dataCell as appropriate
                    case 'Integer'
                        try
                            dataCell{fieldCounter}=cellfun(@(x) sscanf(x,'%d'),tokens);
                        catch
                            dataCell{fieldCounter}=NaN*ones(size(tokens));
                        end
                    case 'Float'
                        try
                            dataCell{fieldCounter}=cellfun(@(x) sscanf(x,'%f'),tokens);
                        catch
                            dataCell{fieldCounter}=NaN*ones(size(tokens));
                        end
                    case 'String'
                        dataCell{fieldCounter}=cellfun(@(x) sscanf(x,'%s'),tokens,'Unif',false);
                    otherwise
                        error('Unsopported Type');
                end
            end
        end
        
        
        
        for eventCounter=1:numberOfEvents
            %vcfResults(dataCounter).data=struct;
            %vcfResults(dataCounter).annoVar=struct;
            for fieldCounter=1:nDataFields
                fLength=fData(fieldCounter).length;
                if(isnan(fLength))
                    fLength=numberOfEvents;
                end
                if(fLength>1)
                    evIdx=eventCounter;
                else
                    evIdx=1;
                end
                
                if(iscell(dataCell{fieldCounter}))
                    %vcfResults(dataCounter).data.(fData(fieldCounter).niceName)=...
                    %    dataCell{fieldCounter}{evIdx};
                    vcfResults(dataCounter).(fData(fieldCounter).niceName)=...
                        dataCell{fieldCounter}{evIdx};
                else
                    %vcfResults(dataCounter).data.(fData(fieldCounter).niceName)=...
                    %    dataCell{fieldCounter}(evIdx);
                    vcfResults(dataCounter).(fData(fieldCounter).niceName)=...
                        dataCell{fieldCounter}(evIdx);
                end
                
                
            end
            
            vcfResults(dataCounter).nEvents=numberOfEvents;
            vcfResults(dataCounter).chrom=chrom;
            vcfResults(dataCounter).pos=pos;
            vcfResults(dataCounter).ref=ref;
            vcfResults(dataCounter).alt=altList{eventCounter};
            vcfResults(dataCounter).id=idList{eventCounter};
            vcfResults(dataCounter).qual=qual;
            vcfResults(dataCounter).filter=filter;
            
            
            tokens=cellfun(@(x) regexp(x,'=','split'), ...
                regexp( infoAnnoVar{eventCounter},';','split'),'Unif',false);
            annoVals=cellfun(@(x) x{2},tokens,'Unif',false);
            
            
            for fieldCounter=1:nAnnoFields
                
                
                %vcfResults(dataCounter).annoVar.(annoData(fieldCounter).niceName)=...
                %        annoVals{fieldCounter};
                vcfResults(dataCounter).(annoData(fieldCounter).niceName)=...
                    annoVals{fieldCounter};
                
            end
            
            
            
            dataCounter=dataCounter+1;
        end
        
        
        
    end
    
    
end
%%

eventID=strcat({vcfResults.chrom},'_',{vcfResults.pos},'_', {vcfResults.ref},'_',{vcfResults.alt});
[num,name]=grp2idx(eventID);
%repEvents=find(diff(num)==0);
%areRepsEqual=false(size(repEvents));
repList=false(size(eventID));
for i=1:length(name)
    if(nnz(num==i)>1)
        eventIdx=find(num==i);
        %areRepsEqual(i)=isequal(vcfResults(repEvents(i)),vcfResults(repEvents(i)+1));
        
        scores=[vcfResults(eventIdx).QT]+0.5*[vcfResults(eventIdx).AF];
        [~,idx]=max(scores);
        eventIdx(idx)=[];
        repList(eventIdx)=true;
    end
    
end

vcfResults(repList)=[];
eventID=strcat({vcfResults.chrom},'_',{vcfResults.pos},'_', {vcfResults.ref},'_',{vcfResults.alt});
end
