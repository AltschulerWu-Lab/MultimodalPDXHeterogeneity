%% This is a script to generate Figure 4: Global variation (whole transcriptome/genome etc) in the PDX samples

figType='Main';
%% Set up Paths & Data
addpath('../');
addpath('../Growth_Curves/');
addpath('../Microscopy/');
addpath('../Pathology/');
addpath('../DNA/');
addpath(genpath('../Common/Plotting/'));
params=GetParams({'all'});

showText=false;

numberOfModels=length(params.samples.modelNames);
numberOfSamples=params.samples.numberOfSamples;
modelOrder=params.figs.modelOrder;% order in which models are plotted.
modelColors=params.figs.modelColors;
[~,modelOrderInverse]=sort(modelOrder);%map between initial order and that used for plotting

%% Load normalized data for the different assays, and calculate distances between pairs

data=cell(2,1);
data{1}=Load_Data({'RNA','DNA','RPPA','microscopy','pathway'},'tumor','normalize',true);
data{2}=Load_Data({'RNA','DNA','RPPA','microscopy','pathway'},'tumor','normalize',false);
dataTypes={'dna','rna','pathway','rppa','microscopy'};
normalizationTypes=[1,1,1,1,1];


% Note: info on different pairs e.g. same mode vs same sector is in
% GetParams
pairInfo=params.samples.pairs;
pairTypes= fieldnames(pairInfo);

corrRes=zeros(length(pairTypes),length(dataTypes));
for dataCounter=1:length(dataTypes)
    profiles=data{normalizationTypes(dataCounter)}.(dataTypes{dataCounter}).data;
    isFlagged=data{normalizationTypes(dataCounter)}.(dataTypes{dataCounter}).flaggedSamples;
    switch(dataTypes{dataCounter})
        case 'dna'
            %mutatedGenes=any(profiles(~isFlagged,:),1);
            %profiles=profiles(:,mutatedGenes);
        case 'microscopy'
            markerPrComponents=@(msNum) (1:30)+(msNum-1)*30;
            prComponentsToUse=arrayfun(markerPrComponents,params.microscopy.markerSetsToUse,'Unif',false);
            prComponentsToUse=horzcat(prComponentsToUse{:});
            profiles=profiles(:,prComponentsToUse);
            meanBlockFractions=mean(profiles,1);
            stdBlockFractions=std(profiles,[],1);
            %profiles=bsxfun(@rdivide,...
            %    bsxfun(@minus,profiles,meanBlockFractions),stdBlockFractions);
    end
  
    for p=1:length(pairTypes)
        pairs=pairInfo.(pairTypes{p});
        corrs=zeros(size(pairs,1),1);
        
        for i=1:length(corrs)
            goodPoints=~isnan(profiles(pairs(i,1),:)) & ~isnan(profiles(pairs(i,2),:)); %Only use non-NAN genes
            
            corrs(i)=corr(profiles(pairs(i,1),goodPoints)',...
                profiles(pairs(i,2),goodPoints)');
            
            if(isFlagged(pairs(i,1))||isFlagged(pairs(i,2)))
                corrs(i)=NaN;
            end
        end
        %pairCorrs.(pairTypes{p})=nanmean(corrs);
        corrRes(p,dataCounter)=nanmean(corrs);
        
    end
end

pairOrder=[1,3,4,2,6,5];
corrRes=corrRes(pairOrder,:);
pairTypes(pairOrder);
%% Fig. 4a correlation across scales and modalities 
f=figure;
imagesc(corrRes,[0,1])
colormap(gray);
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',(0:size(corrRes,2))+0.5,...
    'YTick',(0:size(corrRes,1))+0.5,...
    'GridLineStyle','-','LineWidth',3,'GridAlpha',1);
for x=1:size(corrRes,2)
    for y=1:size(corrRes,1)
        text(x,y,sprintf('%0.2f',corrRes(y,x)),...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'FontSize',20,'Color',[0.75,0.4,0.4],'FontWeight','bold');
    end
end
grid on;
% Needs to add labels for rows and columns
pngFile=fullfile(params.figs.saveDir,'Fig3_PairwiseCorr.png');
%pngFile=fullfile(params.figs.saveDir,'SuppFig_PairwiseCorr_NoNorm.png');
svgFile=regexprep(pngFile,'\.png$','\.svg');
saveas(f,pngFile,'png');
set(f,'Renderer','painters');
saveas(f,svgFile,'svg');
%% Fig 4b: Pathway Heatmap showing intra-model heterogeneity
if(strcmp(figType,'Main'))
    
    pathwayData=Load_Data('pathway','tumor');
    showText=true;
    
    cleanSamples=find(~pathwayData.pathway.flaggedSamples);'mutationTypesToDrop';
    
    pathwayRes=struct;
    pathwayRes.modelNumbers.data=params.figs.modelOrder([params.samples.info(cleanSamples).modelNumber])';
    pathwayRes.pathway.data=pathwayData.pathway.data(cleanSamples,:);
    
    if(showText)
        pathwayRes.modelNumbers.columnLabels={'Model'};
        pathwayRes.pathway.columnLabels=regexprep(pathwayData.pathway.columnLabels,'HALLMARK_','');
    else
        pathwayRes.modelNumbers.columnLabels={};
        pathwayRes.pathway.columnLabels={};
    end
    
    D=pdist(pathwayRes.pathway.data,'correlation');
    tree=linkage(D,'single');
    sampleOrder=optimalleaforder(tree,D);
    D=pdist(pathwayRes.pathway.data','correlation');
    tree=linkage(D,'average');
    colOrder=optimalleaforder(tree,D);
    
    f=figure;
    handles=DataHeatMapPlot(pathwayRes,'rowOrder',sampleOrder,'colOrders',{[1],colOrder},...
        'PlotWidths',[2,50],'colormaps',{params.figs.modelColors,gray},'dataRanges',{[],[-2,2]});
    
    if(~showText)
        handles(1).Title=title('');
        handles(2).Title=title('');
    end
    pngFile=fullfile(params.figs.saveDir,'Fig3_PathwayMat.png');
    svgFile=regexprep(pngFile,'\.png$','\.svg');
    %saveas(f,pngFile,'png');
    %set(f,'Renderer','painters');
    %saveas(f,svgFile,'svg');
    pathwayTable=table([1:length(pathwayData.pathway.columnLabels)]',...
        pathwayData.pathway.columnLabels(colOrder),'VariableNames',...
        {'Column_Number','Pathway_Name'});
    writetable(pathwayTable,'pathwayTable.xls');
    
end


%% Fig. 4c Tumorplot of IF data

sampleMapReordered=params.samples.sampleMap(:,:,params.figs.modelOrder);
sampleMapReordered=sampleMapReordered(:,params.figs.repTumorOrder,:);


if(strcmp(figType,'Main'))
    ifIntensities=load(params.microscopy.avgValueFile);
    %ifIntensities=load('../Microscopy/bgSubtractedIntensities.mat');
    
    tumorIntensities=ifIntensities.meanTumorIntensity;
    ps6If=squeeze(nanmean(tumorIntensities(:,1,3,:),4));
    bcatIf=squeeze(nanmean(tumorIntensities(:,1,2,:),4));
    %ps6If=squeeze(nanmean(tumorIntensities(:,1,2,:),4));
    ps6If=(ps6If-mean(ps6If))/std(ps6If);
    bcatIf=(bcatIf-mean(bcatIf))/std(bcatIf);
    %ps6IfTumor=ps6If(sampleMapReordered);
    f=figure;
    TumorPlot(ps6If(sampleMapReordered),-1,1,gray,[]);
    f=figure;
    TumorPlot(bcatIf(sampleMapReordered),-1,1,gray,[]);
    
    
    pngFile=fullfile(params.figs.saveDir,'Fig3_pS6_IF_TumorPlot.png');
    svgFile=regexprep(pngFile,'\.png$','\.svg');
    %saveas(f,pngFile,'png');
    set(f,'Renderer','painters');
    %saveas(f,svgFile,'svg');
   
    
end

%% Fig.4c Example pS6 Tumor Images
if(strcmp(figType,'Main'))
    %scaledImageList=load(params.microscopy.bgSubtractedImgList); 
    scaledImageList=LoadScaledImgList();
    scaledImageList=scaledImageList.scaledImageList;
    
    modelNumber=params.figs.modelOrder(2); %second row
    tumorNumber=2;%middle tumor
    
    regionNames={'Dorsal','Central','Ventral'};
    repNumber=[2,1,3];
    channelsToShow=[1,3,4];
    markerSet=1;
    markerScales=[0,30000;0,30000;0,40000;0,20000];
    %markerScales=params.microscopy.prBgSubMarkerScales;
    for region=1:3
        sampleNum=params.samples.sampleMap(region,tumorNumber,modelNumber);
        %imgList=GetImages('IF','sampleNumber',sampleNum,'markerSet',markerSet);
        %tImg=Afi2TImg(imgList{repNumber});
        %img=tImg.LoadImage('MagLevel',3);
        
        tImg=scaledImageList{markerSet,sampleNum}{repNumber(region)};
        info=GetImagingInfo(tImg.imgFileMat{1});
        mpp=info.MPP*mean(tImg.dimensions(1,:)./tImg.dimensions(3,:));
        [~,pixelRegion]=Identify_Tumor_Mask(tImg);
        img=tImg.LoadImage('MagLevel',3,'PixelRegion',pixelRegion);
        
        
        f=figure;
        
        Display_Image_With_ScaleBar(img(:,:,channelsToShow),markerScales(channelsToShow,:),mpp,...
            'scaleBarInMicrons',600,'scaleBarYPos',0.95,'scaleBarXPos',0.05,'showScaleBarText',false,...
            'Colors',{'B','R','K'});
        axis off equal;
        
        pngFile=fullfile(params.figs.saveDir,['Fig3d_IF_Image_' regionNames{region} '.png']);
        svgFile=regexprep(pngFile,'\.png$','\.svg');
        %saveas(f,pngFile,'png');
        %set(f,'Renderer','painters');
        %saveas(f,svgFile,'svg');
        
    end
end
