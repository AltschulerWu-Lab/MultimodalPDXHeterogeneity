%% Fig 3a&b: DNA/RNA Quality
addpath(genpath('../Common/Plotting/'));
addpath('../');
params=GetParams();
showLabels=true;

%% Fig.3a: DNA Quality
numberOfSamples=36;
chrData=cell(numberOfSamples,1);
ampData=cell(numberOfSamples,1);

% Loop over sampled and load DNA quality information
dnaQualDir=fullfile(params.rna.rootDir,'DNA_Processed');
for sampleCounter=1:numberOfSamples
    chrFile=fullfile(dnaQualDir,['IonXpress_' sprintf('%0.3d',sampleCounter) ...
        '_R_2017_07_06_11_23_36_user_ionproton-213-UCSF-CHPv2_Auto_user_' ...
        'ionproton-213-UCSF-CHPv2_369.chr.cov.txt']);
    ampFile=fullfile(dnaQualDir,['IonXpress_' sprintf('%0.3d',sampleCounter) ...
        '_R_2017_07_06_11_23_36_user_ionproton-213-UCSF-CHPv2_Auto_user_' ...
        'ionproton-213-UCSF-CHPv2_369.amplicon.cov.txt']);
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    chrData{sampleCounter}=readtable(chrFile);
    ampData{sampleCounter}=readtable(ampFile);
    
end

numberMappedReads=arrayfun(@(x) sum(chrData{x}.total_reads),1:numberOfSamples);
pctAssignedReads=arrayfun(@(x) 100*sum(ampData{x}.total_reads)/sum(chrData{x}.total_reads),1:numberOfSamples);


figHandle=figure;
minMappedReads=1E6;
minAssignedReads=85;

sampleQuality=numberMappedReads>minMappedReads & pctAssignedReads>minAssignedReads;
scatter(numberMappedReads,pctAssignedReads,150,...
    sampleQuality,'filled','MarkerEdgeColor','k','LineWidth',2);
colormap([0.85 0.4 0.4; 0.4 0.85 0.4 ]);
hold on;
xLim=[1,7E6];
yLim=[0,100];

xTick=[0,minMappedReads, 3E6,6E6];

plot(xLim,minAssignedReads*[1,1],'--k');
plot(minMappedReads*[1,1],yLim,'--k');
axis square;
set(gca,'FontName','Arial','FontSize',14,'YTick',0:50:125,...
   'XTick',[0,minMappedReads, 3E6,6E6],'XTickLabel',xTick*1E-6,...
   'XLim',xLim,'YLim',yLim,'XScale','linear');

set(gca,'YTick',0:50:125,'XTick',[0,minMappedReads, 3E6,6E6],...
   'XLim',xLim,'YLim',yLim,'XScale','linear');
if(showLabels)
    set(gca,'FontName','Arial','FontSize',14,'XTickLabel',xTick*1E-6);
    xlabel('Number Of Mapped Reads (\times 10^6)','FontName','Arial','FontSize',18);
    ylabel('% Reads Mapped to Amplicons','FontName','Arial','FontSize',18);
else
    set(gca,'XTickLabel',[],'YTickLabel',[]);
end

%% Fig 3b: RNA quality
rnaQuality=readtable(params.rna.qualityFile,'Delimiter','\t');
sampleQuality=ones(params.samples.numberOfSamples+1,1);
sampleQuality(params.rna.badSamples)=0;
sampleQuality=sampleQuality(1:params.samples.numberOfSamples);
figHandle=figure;
scatter(rnaQuality.Pct_Targets_Detected,rnaQuality.Mean_Read_Length,150,...
    sampleQuality,'filled','MarkerEdgeColor','k','LineWidth',2);
colormap([0.85 0.4 0.4; 0.4 0.85 0.4 ]);
hold on;
xLim=[0,60];
yLim=[50,125];
minReadLength=100;
minPctTargets=45;
plot(xLim,minReadLength*[1,1],'--k');
plot(minPctTargets*[1,1],yLim,'--k');
axis square;


set(gca,'XTick',[0,20,40,minPctTargets,60],'YTick',50:25:125,...
    'XLim',xLim,'YLim',yLim);
if(showLabels)
    set(gca,'FontName','Arial','FontSize',14);
    xlabel('Percentage of Genes Mapped','FontName','Arial','FontSize',18);
    ylabel('Average Read Length','FontName','Arial','FontSize',18);
else
    set(gca,'XTickLabel',[],'YTickLabel',[]);
end

