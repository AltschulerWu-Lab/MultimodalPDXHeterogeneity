%% Code for Fig3c: Mouse contribution to DNA

addpath('../');
addpath('../Growth_Curves/');
addpath('../Pathology/');
addpath('../DNA/');
addpath(genpath('../Common/Plotting/'));
params=GetParams({'all'});

%showText=false;

numberOfModels=length(params.samples.modelNames);
numberOfSamples=params.samples.numberOfSamples;
modelOrder=params.figs.modelOrder;% order in which models are plotted.
modelColors=params.figs.modelColors;
[~,modelOrderInverse]=sort(modelOrder);%map between initial order and that used for plotting

%% Load Data

% Load all mutation data
dnaData=load(params.dna.processedResultsFile,...
    'uniqueIDs','isMutated','eventInfo','altAlleleFraction','eventNames');
% Filtering read: in this case we don't do much
[dnaSamplesToUse,dnaMutationsToUse]=FilterDnaMutations(dnaData,...
    'keepOnlyCosmic',false,'minVAF',0,'dropNOCALL',false,...
    'mutationTypesToDrop',{'none'},'filterMouse',false);
% Get the variant allele frequency for each event
vaf=dnaData.altAlleleFraction(dnaSamplesToUse,dnaMutationsToUse);

% Identify mutations that are called in sample 37, pure mouse
isMutMouse=dnaData.isMutated(37,dnaMutationsToUse)|dnaData.altAlleleFraction(37,dnaMutationsToUse)>0;




%% Make Plot
mutColors=0.75*[0.5 0.5 1; 1 0.5 0.5];
modelNumbers=[params.samples.info.modelNumber];
figure;
axes('Position',[0.0,0.1,1,0.9]);
% Plot VAF of sites across samples (each curve is a different mutation site site)
[~,order]=sort(prctile(vaf(:,isMutMouse),95,2)); % Order samples essentially by increasing VAF

% Plot the sites mutated an un-mutated in pure mouse in different colors
plot(vaf(order,~isMutMouse),'Color',mutColors(1,:),'LineWidth',2);
hold on;
plot(vaf(order,isMutMouse),'Color',mutColors(2,:),'LineWidth',2);
hold off;

set(gca,'XLim',[0.5,length(dnaSamplesToUse)+0.5],'YLim',[0,1.01],...
    'YTick',[0,0.5,1],'TickLength',[0,0.05]);
axes('Position',[0.0,0.0,1,0.1]);

% Make colored bar at bottom denotin the model each sample came from
imagesc(modelNumbers(dnaSamplesToUse(order)));
set(gca,'XTick',[1:length(dnaSamplesToUse)]-0.5,'YTick',[],'TickLength',[0,0]);
grid on;
colormap(params.figs.modelColors(params.figs.modelOrder,:));


