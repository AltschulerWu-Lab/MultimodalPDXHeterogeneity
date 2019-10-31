function [sampleGrowthCurves,doublingTimes,fitInitVol]=GetSampleGrowthData()

params=GetParams({'growth','samples'});
% OUTPUT is the growth curves as a cell array, with cells ordered
% by sample number, each containing the growth curves for that sample.
% Each growth curve is a 2 column matrix with the rows being different
% time points, the first column is day number and second column is
% volume in mm^3
sampleGrowthCurves=cell(params.samples.numberOfSamples,1);


% The code below assumes a specific order of models. Check to make sure
% upstream code does not break this.
if(~isequal(params.samples.modelNames,{'CN1571','CN1574','CN1572','CR-0104-O'}'))
    error('Model Names Dont Match Assumptions. Results Cannot be trusted');
end

tumorGrowth=Load_Tumor_Volumes(true);

fitInitVol=zeros(params.samples.numberOfSamples,1);
doublingTimes=zeros(params.samples.numberOfSamples,1);
for sampleNumber=1:params.samples.numberOfSamples
    info=params.samples.info(sampleNumber);
    model=info.modelNumber;
    animalID=info.mouseNum;
    mouseFlank=info.flank;
    growth=tumorGrowth(model);
    isSample= (growth.animalIds==animalID)&(strcmp(growth.flank,mouseFlank));
    fitInitVol(sampleNumber)=growth.initVol(isSample);
    doublingTimes(sampleNumber)=growth.doublingTime(isSample);
    days=growth.days;
    volume=growth.volume(isSample,:);
    sampleGrowthCurves{sampleNumber}=[days(~isnan(volume))',volume(~isnan(volume))'];
    
end

% OLD CODE
%     % Sheets for L/R flanks (side of mouse on which tumor is embedded).
%     % For data structure consistency, models where only one
%     % flank was used are denoted as having the same sheet for both flanks.
%     sheetsFromFlanks={[5,5],[5,5],[5,8],[5,8]};
%     % Range in excel sheet (specified above) containing tumor volume data
%     modelDataRange={'B9:S60','B9:AA60','B9:N15','B9:P15'};
%     % Rows and columns within specified range that specify date row and
%     % start of volume data. The file format used means these are preserved across models
%     dayRow=1;
%     volStartCol=8;
%
%
%     for modelNumber=1:length(params.samples.modelNames)
%
%         modelName=params.samples.modelNames{modelNumber};
%         dataFile=params.growth.model2File(modelName);
%         samplesFromModel=find(strcmp({params.samples.info.modelName},modelName));
%         mouseNumbersFromModel=[params.samples.info(samplesFromModel).mouseNum]';
%         mouseFlanksFromModel=strcmp({params.samples.info(samplesFromModel).flank},'R')'+1;
%
%         [sheetIdx,sheetNumbers]=grp2idx(sheetsFromFlanks{modelNumber}(mouseFlanksFromModel));
%
%         for sheetCounter=1:length(sheetNumbers)
%            sheetNumber=str2double(sheetNumbers(sheetCounter));
%            growthData=readtable(dataFile,'Sheet',sheetNumber,'Range',modelDataRange{modelNumber});
%
%            animalIds=str2double(growthData.AnimalID);
%            isInSheet=(sheetIdx==sheetCounter);
%            [~,rowNumbers]=ismember(mouseNumbersFromModel(isInSheet),animalIds);
%            volumeData=table2array(growthData(rowNumbers,volStartCol:end));
%            dayMeasured=table2array(growthData(dayRow,volStartCol:end));
%
%            samplesInSheet=samplesFromModel(isInSheet);
%            for sampleCounter=1:length(samplesInSheet)
%                sampleNum=samplesInSheet(sampleCounter);
%                growthCurves=[dayMeasured;volumeData(sampleCounter,:)]';
%                sampleGrowthCurves{sampleNum}=growthCurves(all(~isnan(growthCurves),2),:);
%
%            end
%         end
%
%
%     end
%     [doublingTimes,fitInitVol] =cellfun(@CalculateDoublingTime,sampleGrowthCurves);

end

function [doublingTime,initVol]=CalculateDoublingTime(growthCurve)
t=growthCurve(:,1);
logVol=log(growthCurve(:,2));
t=t(isfinite(logVol));
logVol=logVol(isfinite(logVol));
regRes=[ones(size(t)),t]\logVol;
slope=regRes(2);
intercept=regRes(1);
doublingTime=log(2)/slope;
initVol=exp(intercept);
if(isnan(doublingTime))
    error('Something Went Wrong In Calculating Doubling Times');
end
end