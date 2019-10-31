%% This is a script to generate Figure 2: Non-molecular
% (i.e. growth and pathology) based comparison of the PDX samples

%% Set up Paths & Data
addpath('../');
addpath('../Growth_Curves/');
addpath('../Pathology/');
addpath('../Common/Plotting/');
addpath('../Common/Plotting/ColorBrewer/');
params=GetParams({'samples','growth','pathology','figs'});
numberOfModels=length(params.samples.modelNames);
numberOfSamples=params.samples.numberOfSamples;
modelOrder=params.figs.modelOrder;% order in which models are plotted.
modelColors=params.figs.modelColors;
[~,modelOrderInverse]=sort(modelOrder);%map between initial order and that used for plotting
%% Load Growth Curve Data
% Use only the deposited tumors (i.e. 3% per model).
useOnlyDepositedTumors=true;
showText=true;

% Specify the range of tumor volumes to consider
% Some tumors never really grew, and larger ones might saturate out.
minTumorVolume=400;
maxTumorVolume=800;


% Load tumor growth data organized by model
tumorGrowth=Load_Tumor_Volumes(useOnlyDepositedTumors);
doublingTimes=cell(numberOfModels,1);
for modelCounter=1:numberOfModels
    tumorsInSizeRange=find(tumorGrowth(modelCounter).finalVol>400 &tumorGrowth(modelCounter).finalVol<800);
    doublingTimes{modelCounter}=tumorGrowth(modelCounter).doublingTime(tumorsInSizeRange);
end

%% Fig. 2b: Boxplot of Doubling Times 
f=figure;

%modelOrderInverse=1:4;
modelNumbers=arrayfun(@(x) modelOrderInverse(x)*ones(length(doublingTimes{x}),1),1:numberOfModels,'Unif',false);
h=boxplot(vertcat(doublingTimes{:}),vertcat(modelNumbers{:}),...
    'PlotStyle','traditional','Colors',modelColors);
if(useOnlyDepositedTumors)
    yLim=[5,15];
else
    yLim=[5,20];
end
set(gca,'YLim',yLim,'XTickLabel',[],'YTick',[5,10,15,20]);
set(h,{'linew'},{2})
if(~showText)
    set(gca,'YTickLabel',[]);
end
pngFile=fullfile(params.figs.saveDir,'Doubling_Times_BarPlot.png');
svgFile=regexprep(pngFile,'\.png$','\.svg');
saveas(f,pngFile,'png');
saveas(f,svgFile,'svg');
%% ANOVA to Show Different Models Have Different Doubling Times
p=anova1(vertcat(doublingTimes{:}),vertcat(modelNumbers{:}));




%% Get final tumor volumes and histopathological parameters
[pctTumor,pctStroma,pctNecrosis,heScores]=GetPathologyInfo();
finalTumorVolumes=zeros(numberOfSamples,1);
for sampleNumber=1:numberOfSamples
    info=params.samples.info(sampleNumber);
    model=info.modelNumber;
    animalID=info.mouseNum;
    mouseFlank=info.flank;
    growth=tumorGrowth(model);
    finalTumorVolumes(sampleNumber)=growth.finalVol(...
        (growth.animalIds==animalID)&(strcmp(growth.flank,mouseFlank)));
    
end

%% Fig. 2c: Make Ternary Plot of % Tumor/Necrosis/Stroma
f=figure;
TernaryPlot(pctNecrosis/100,pctTumor/100,pctStroma/100,...
    modelColors(modelOrderInverse([params.samples.info.modelNumber]),:));
if(showText)
    text(0,0,'100% Necrosis','Rotation',-60,'HorizontalAlignment','center',...
        'VerticalAlignment','Top','FontName','Arial');
    text(1000,0,'100% Tumor','Rotation',60,'HorizontalAlignment','center',...
        'VerticalAlignment','Top','FontName','Arial');
    text(500,1000*sqrt(3)/2,'100% Stroma','Rotation',0,'HorizontalAlignment','center',...
        'VerticalAlignment','Bottom','FontName','Arial');
end
pngFile=fullfile(params.figs.saveDir,'Pathology_Ternary.png');
svgFile=regexprep(pngFile,'\.png$','\.svg');
saveas(f,pngFile,'png');
saveas(f,svgFile,'svg');
%% Fig 2d: Scatter plot of %Necrosis vs Tumor Volume
f=figure;
scatter(finalTumorVolumes,pctNecrosis,100,...
    modelOrderInverse([params.samples.info.modelNumber]),'filled');
colormap(modelColors);
set(gca,'YLim',[0,100],'XLim',[400,800],'XTick',[400,600,800],...
    'YTick',[0,50,100]);
if(showText)
    xlabel('Tumor Volume (mm^3)');
    ylabel('% Necrosis');
else
    set(gca,'XTickLabel',[],'YTickLabel',[]);
end
pngFile=fullfile(params.figs.saveDir,'Necrosis_Vs_Volume.png');
svgFile=regexprep(pngFile,'\.png$','\.svg');
saveas(f,pngFile,'png');
saveas(f,svgFile,'svg');

