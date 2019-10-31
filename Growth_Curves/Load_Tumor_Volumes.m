function  tumorGrowth=Load_Tumor_Volumes(onlyUseDepositedTumors)

% Returns a structure array (with one element per model) quantifying
% various parameters of tumor growth.
% onlyUseDepositedTumors - if this is set to true (default), only the 3
% tumors deposited are used, whereas if false, all available tumors for
% each model are used.
%%
if (nargin<1)
    onlyUseDepositedTumors=true;
end
params=GetParams({'growth','samples'});

numberOfModels=length(params.samples.modelNames);
sampleInfo=params.samples.info;
mouseFlankCombos=strcat(cellfun(@num2str,{sampleInfo.mouseNum},'Unif',false)',...
    '_',{sampleInfo.flank}');
%%
model2MouseFlank=cell(numberOfModels,1);
for modelNumber=1:numberOfModels
    samplesInModel=[sampleInfo.modelNumber]==modelNumber;
    temp=regexp(unique(mouseFlankCombos(samplesInModel,:)),'_','split');
    temp=vertcat(temp{:});
    model2MouseFlank{modelNumber}=cell(size(temp));
    model2MouseFlank{modelNumber}(:,1)=cellfun(@str2double,temp(:,1),'Unif',false);
    model2MouseFlank{modelNumber}(:,2)=temp(:,2);
end

%%

xlsxDataRaw=cell(numberOfModels,1);
volumeInfo=cell(numberOfModels,1);

tumorGrowth=struct;
for modelCounter=1:numberOfModels
    [~,sheetNames]=xlsfinfo(params.growth.model2File(params.samples.modelNames{modelCounter}));
    tumorVolumeSheets=find(~cellfun(@isempty,regexp(sheetNames,'IM_Pool_I[1,2]I','match')));
    xlsxDataRaw{modelCounter}=cell(length(tumorVolumeSheets),1);
    volumeInfo{modelCounter}=struct;
    if(length(tumorVolumeSheets)==1)
        flankList={'R'};
    elseif(length(tumorVolumeSheets)==2)
        flankList={'L','R'};
    else
        error('Something went wrong in reading the growth curve sheets');
    end
    for sheetCounter=1:length(tumorVolumeSheets)
        xlsxDataRaw{modelCounter}{sheetCounter}=readtable(...
            params.growth.model2File(params.samples.modelNames{modelCounter}),...
            'Sheet',sheetNames{tumorVolumeSheets(sheetCounter)});
        daysRow=find(strcmp(xlsxDataRaw{modelCounter}{sheetCounter}{:,8},'Day'),1);
        lastDataRow=find(cellfun(@isempty,xlsxDataRaw{modelCounter}{sheetCounter}{(daysRow+1):end,9}),1)-1;
        dataRows=daysRow+(1:lastDataRow);
        
        volumeInfo{modelCounter}(sheetCounter).volumeData=str2double(table2array(xlsxDataRaw{modelCounter}{sheetCounter}(dataRows,9:end)));
        volumeInfo{modelCounter}(sheetCounter).days=str2double(table2array(xlsxDataRaw{modelCounter}{sheetCounter}(daysRow,9:end)));
        volumeInfo{modelCounter}(sheetCounter).animalIds=str2double(table2array(xlsxDataRaw{modelCounter}{sheetCounter}(dataRows,5)));
        
        volumeInfo{modelCounter}(sheetCounter).flank=cell(length(dataRows),1);
        volumeInfo{modelCounter}(sheetCounter).flank(:)=flankList(sheetCounter);
    end
    
    
    tumorGrowth(modelCounter).animalIds=vertcat(volumeInfo{modelCounter}.animalIds);
    tumorGrowth(modelCounter).flank=vertcat(volumeInfo{modelCounter}.flank);
    
    
    tumorGrowth(modelCounter).volume=vertcat(volumeInfo{modelCounter}.volumeData);
    tumorGrowth(modelCounter).days=vertcat(volumeInfo{modelCounter}.days);
    if(size(tumorGrowth(modelCounter).days,1)>1)
        if(isequal(tumorGrowth(modelCounter).days(1,:),tumorGrowth(modelCounter).days(2,:)))
            tumorGrowth(modelCounter).days=tumorGrowth(modelCounter).days(1,:);
        else
            error('Mismatch in days');
        end
    end
    
    
    if(onlyUseDepositedTumors)
        tumorsToUse=zeros(size(model2MouseFlank{modelCounter},1),1);
        for tumorCounter=1:size(model2MouseFlank{modelCounter})
            tumorsToUse(tumorCounter)=find(...
                tumorGrowth(modelCounter).animalIds==model2MouseFlank{modelCounter}{tumorCounter,1}&...
                strcmp(tumorGrowth(modelCounter).flank,model2MouseFlank{modelCounter}{tumorCounter,2}));
        end
    else
        tumorsToUse=1:length( tumorGrowth(modelCounter).animalIds);
    end
    
    tumorGrowth(modelCounter).animalIds=tumorGrowth(modelCounter).animalIds(tumorsToUse);
    tumorGrowth(modelCounter).flank=tumorGrowth(modelCounter).flank(tumorsToUse);
    tumorGrowth(modelCounter).volume=tumorGrowth(modelCounter).volume(tumorsToUse,:);
    
    tumorGrowth(modelCounter).doublingTime=zeros(size(tumorGrowth(modelCounter).animalIds));
    tumorGrowth(modelCounter).initVol=zeros(size(tumorGrowth(modelCounter).animalIds));
    tumorGrowth(modelCounter).finalVol=zeros(size(tumorGrowth(modelCounter).animalIds));
    for tumorCounter=1:length(tumorGrowth(modelCounter).animalIds)
        growthCurve=[tumorGrowth(modelCounter).days',...
            tumorGrowth(modelCounter).volume(tumorCounter,:)'];
        [tumorGrowth(modelCounter).doublingTime(tumorCounter),...
            tumorGrowth(modelCounter).initVol(tumorCounter),...
            tumorGrowth(modelCounter).finalVol(tumorCounter)]=...
            CalculateDoublingTime(growthCurve);
    end
    
    
    
end


end


function [doublingTime,initVol,finalVol]=CalculateDoublingTime(growthCurve)
%% Calculate the doubling time etc given the growth curve
t=growthCurve(:,1);
logVol=log(growthCurve(:,2));
t=t(isfinite(logVol));
logVol=logVol(isfinite(logVol));
regRes=[ones(size(t)),t]\logVol;
slope=regRes(2);
intercept=regRes(1);
doublingTime=log(2)/slope;
initVol=exp(intercept);
finalVol=exp(logVol(end));
if(isnan(doublingTime))
    error('Something Went Wrong In Calculating Doubling Times');
end
end