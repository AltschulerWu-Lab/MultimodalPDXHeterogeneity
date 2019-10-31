%% Code for Fig3d: Mouse contribution to RNA

addpath('../');
addpath('../Growth_Curves/');
addpath('../Microscopy/');
addpath('../Pathology/');
addpath(genpath('../Common/Plotting/'));
addpath('../RNA/');
params=GetParams();
numberOfSamples=params.samples.numberOfSamples;
%% Load RNA Data for original experiment

[combinedVals,singleChipVals,geneInfo,combinedData]=Load_Raw_RNA();
numberOfChips=size(singleChipVals,2);
%% Load RNA Data for mouse spike-in experiment

% Data for sample 3 in the spike-in experiment
spike3File=fullfile(params.rna.rootDir,'RNA_SpikeIn_Sample3',...
    'IonXpress_021_R_2017_06_28_14_31_35_LCM3_RNA.txt');
% Data for sample 3+mouse in the spike-in experiment
spike3PlusMouseFile=fullfile(params.rna.rootDir,'RNA_SpikeIn_Sample3_plus_mouse',...
    'IonXpress_022_R_2017_06_28_14_31_35_LCM3_+_mouse.cov.txt');

spike3Data=sortrows(readtable(spike3File),'attributes');
spike3Vals=spike3Data.overlaps;
spike3Vals=1E6*spike3Vals/sum(spike3Vals);

spike3PlusMouseData=sortrows(readtable(spike3PlusMouseFile),'attributes');
spike3PlusMouseVals=spike3PlusMouseData.overlaps;
spike3PlusMouseVals=1E6*spike3PlusMouseVals/sum(spike3PlusMouseVals);

% Values for Sample 3 in the original experiment
orig3Vals=1E6*combinedVals(:,3)/sum(combinedVals(:,3));
origVals=1E6*bsxfun(@rdivide,combinedVals,sum(combinedVals,1));
badSamples=false(numberOfSamples,1); badSamples(params.rna.badSamples)=true;

lowExpressionThreshold=30;
lowExpressionGenes=mean(origVals(:,~badSamples),2)<lowExpressionThreshold & sum(origVals(:,~badSamples)==0,2)>1;

spikeVals=[spike3Vals,spike3PlusMouseVals];
avgSpikeVals=mean(spikeVals,2);
deltaSpikeVals=diff(spikeVals,1,2);
lowExpressionGenes=avgSpikeVals<lowExpressionThreshold | any(spikeVals==0,2);
%%
mainExptBiolColor=[0 0 0];
mainExptTechColor=[0.8,0.8,0.8];
spikeInExptColor=[0.2,0.8,0.2];

%% Global Effect of spike In
% Compare correlation between 3 vs 3+mouse from the spike in experiment to
% distributions of correlations from
% 1) Technical variability (across chips for the same sample)
% 2) Sample Variability (across samples from the same model)


% Step I: identify which genes to use
lowExpressionThreshold=30;
lowExpressionGenes=mean(origVals(:,~badSamples),2)<lowExpressionThreshold | ...
    sum(origVals(:,~badSamples)==0,2)>1;

% Step II: Normalization function
meanGeneLevels=mean(origVals(:,~badSamples),2);
stdGeneLevels=std(origVals(:,~badSamples),[],2);
normalize =@(profile) bsxfun(@rdivide,bsxfun(@minus,profile(~lowExpressionGenes),...
    meanGeneLevels(~lowExpressionGenes)),stdGeneLevels(~lowExpressionGenes));
singleChipRPM=1E6*bsxfun(@rdivide,singleChipVals,sum(singleChipVals,1));

%Step III: Calculate technical variation
numberOfChipPairs=nchoosek(numberOfChips,2);
techVarCorrs=zeros(numberOfSamples,numberOfChipPairs);

for sampleCounter=1:numberOfSamples% Loop over samples
    pairCounter=1;
    for chip1=1:numberOfChips% Loop over all pairs of chips within sample
        for chip2=1:(chip1-1)
            techVarCorrs(sampleCounter,pairCounter)=...
                corr(normalize(singleChipRPM(:,chip1,sampleCounter)),...
                normalize(singleChipRPM(:,chip2,sampleCounter)));
            pairCounter=pairCounter+1;
        end
    end
end


%%  Make Plots
%Are a subset of genes screwed up (look at effect of tail)

%Single sample
showLabels=true;


xLim=[30,1000];
if(showLabels)
    figure;
    subplot(1,2,1);
else
    figure;
end

% Plot comparing two different chips from sample 3: i.e. technical
% variation
sampleNumber=3;
chipNumbers=[1,2];
vals=squeeze(singleChipRPM(:,chipNumbers,sampleNumber));
isLow=mean(vals,2)<30|any(vals==0,2);
vals=vals(~isLow,:);
scatter(vals(:,1),vals(:,2),50,mainExptTechColor,'filled',...
    'MarkerEdgeColor','k');
hold on;
scores=diff(vals,1,2)./mean(vals,2);
scatter(vals(scores>0.5,1),vals(scores>0.5,2),60,'w','filled',...
    'MarkerEdgeColor','r','MarkerFaceAlpha',0,'LineWidth',3,'MarkerEdgeAlpha',0.5);
xVals=linspace(xLim(1),xLim(2),100);
plot(xVals,1.25/0.75*xVals,'--r');
plot(xVals,0.75/1.25*xVals,'--k');
hold off;
set(gca,'XScale','log','YScale','log','XLim',xLim,'YLim',[30,1000]);
if(showLabels)
    set(gca,'FontSize',18,'FontName','Arial');
    xlabel('Sample 3 Chip 1 Counts');
    ylabel('Sample 3 Chip 2 Counts');
    title('Main Experiment');

else
    set(gca,'XTickLabel',[],'YTickLabel',[]);
end
axis square;

% Plot comparing sample 3 to sample 3+mouse
if(showLabels)
    subplot(1,2,2);
else
    figure;
end
vals=spikeVals;
isLow=mean(vals,2)<30|any(vals==0,2);
vals=vals(~isLow,:);
scatter(vals(:,1),vals(:,2),50,spikeInExptColor,'filled',...
    'MarkerEdgeColor','k');
hold on;
scores=diff(vals,1,2)./mean(vals,2);
scatter(vals(scores>0.5,1),vals(scores>0.5,2),60,'w','filled',...
    'MarkerEdgeColor','r','MarkerFaceAlpha',0,'LineWidth',3,'MarkerEdgeAlpha',0.5);
xVals=linspace(xLim(1),xLim(2),100);
plot(xVals,1.25/0.75*xVals,'--r');
plot(xVals,0.75/1.25*xVals,'--k');
hold off;
set(gca,'XScale','log','YScale','log','XLim',xLim,'YLim',[30,1000]);
if(showLabels)
    set(gca,'FontSize',18,'FontName','Arial');
    xlabel('Sample 3 Counts');
    ylabel('Sample (3 + 10% Mouse) Counts');
    title('Spike-In Experiment');
else
    set(gca,'XTickLabel',[],'YTickLabel',[]);
end
axis square;

