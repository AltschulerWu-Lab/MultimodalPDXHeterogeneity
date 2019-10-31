% Code for Figs 5a/b: Decomposition of variation across spatial scales
addpath('../');
addpath('../Microscopy/');
addpath('../Microscopy/Spatial/');
addpath(genpath('../Common/'));
params=GetParams({'microscopy','samples','figs'});
%% Load Downsampling Data: Note this assumes the calculations have already been performed
% This block will give us the response and the components of it assigned 
% to different length scale for every sample x replicate for each marker
% (we only use 2 per marker set, since we exclude the nuclear and stromal
% markers) both for the nuclear and cytoplasmic regions. For each of these
% we store the values at 10,000 randomly sampled pixels.

resultsDir= params.microscopy.downsamplingResultsDir;
numberOfReplicates=3;
numberOfSamples=params.samples.numberOfSamples;
numberOfMarkerSets=length(params.microscopy.markerSets);
markersToProfile=[2,3];% Each marker set contains 4 markers, but the 1st and 4th are DAPI/vimentin which we do not profile
markersPerSet=length(markersToProfile); 
numberOfMarkers=numberOfMarkerSets*markersPerSet;
numberOfRegions=2;% 1) nuclear + 2) cytoplasmic
numberOfPoints=10000; % This is the number of randomly sampled pixels per section

%temp=load(params.microscopy.bgSubtractedImgList);
temp=LoadScaledImgList();
imgList=temp.scaledImageList; clear('temp');

%trueData contains the actual intensity distribution of the image
trueData=cell(numberOfSamples,numberOfMarkers,numberOfReplicates,numberOfRegions);
% downsampleData contains the components of the intensity attributed to the
% different length scales
downsampleData=cell(numberOfSamples,numberOfMarkers,numberOfReplicates,numberOfRegions);
isDataPresent=false(numberOfSamples,numberOfMarkers,numberOfReplicates,numberOfRegions);
for sampleNumber=1:numberOfSamples
    for replicateNumber=1:numberOfReplicates
        resultFilename=fullfile(resultsDir,...
            ['downSampler-Sample' num2str(sampleNumber) ...
            '-Replicate-' num2str(replicateNumber) '.mat']);
        if(exist(resultFilename,'file'))
            data=load(resultFilename);
            for markerSet=1:numberOfMarkerSets
                for markerCounter=1:markersPerSet
                    markerNumber=(markerSet-1)*markersPerSet+markerCounter;
                    for regionNumber=1:numberOfRegions
                        downsampleData{sampleNumber,markerNumber,replicateNumber,regionNumber}=...
                            data.downsampleData{markerSet,markerCounter,regionNumber};
                        trueData{sampleNumber,markerNumber,replicateNumber,regionNumber}=...
                            data.trueData{markerSet,markerCounter,regionNumber};
                        if(~isempty(trueData{sampleNumber,markerNumber,replicateNumber,regionNumber}))
                            isDataPresent(sampleNumber,markerNumber,replicateNumber,regionNumber) =true;
                        end
                        
                    end
                    
                end
            end
            
            
        end
    end
end

%% Reorganize data into a single matrix per marker/region
% We reorganize the data, so as to combine all pixels for a single marker 
% (but keeping nuclear and cytoplasmic regions separate), while also
% keeping track of which model/tumor/sector etc each of these pixels came
% from. This makes it easier to break down the variation below.

downsampleMats=cell(numberOfMarkers,numberOfRegions);
trueMats=cell(numberOfMarkers,numberOfRegions);
annoMats=cell(numberOfMarkers,numberOfRegions); % This contains information about model/tumor/sector/replicate_slide that each pixel came from
modelNumbers=[params.samples.info.modelNumber];
tumorNumbers=[params.samples.info.repTumorNum];
sectorNumbers=[params.samples.info.regionNumber];

for regionNumber=1:numberOfRegions
    for markerNumber=1:numberOfMarkers
        numberProfiled=nnz(squeeze(isDataPresent(:,markerNumber,:,regionNumber)))*numberOfPoints;
        annoMats{markerNumber,regionNumber}=zeros(numberProfiled,4);
        downsampleMats{markerNumber,regionNumber}=zeros(numberProfiled,5);
        trueMats{markerNumber,regionNumber}=zeros(numberProfiled,1);
        
        dataCounter=0;
        for sampleNumber=1:numberOfSamples
            for replicateNumber=1:numberOfReplicates
                if(isDataPresent(sampleNumber,markerNumber,replicateNumber,regionNumber))
                    anno=zeros(numberOfPoints,4);
                    anno(:,1)=modelNumbers(sampleNumber);
                    anno(:,2)=tumorNumbers(sampleNumber);
                    anno(:,3)=sectorNumbers(sampleNumber);
                    anno(:,4)=replicateNumber;
                    
                    annoMats{markerNumber,regionNumber}(dataCounter+(1:numberOfPoints),:)=anno;
                    trueMats{markerNumber,regionNumber}(dataCounter+(1:numberOfPoints),:)=...
                        trueData{sampleNumber,markerNumber,replicateNumber,regionNumber};
                    downsampleMats{markerNumber,regionNumber}(dataCounter+(1:numberOfPoints),:)=...
                        downsampleData{sampleNumber,markerNumber,replicateNumber,regionNumber};
                    
                    dataCounter=dataCounter+numberOfPoints;
                end
                
            end
        end
        
        
        
    end
end


%% Perform Hierarchical Decomposition for Non-Spatial Scales
% The Spatial Hierarchical decomposition within a slide as above 
% was achieved using the spatial relationships between pixels.
% Using gaussain convolution we averaged at longer length scales and subtracted out
% average responses with residuals being passed on to lower length scales.
% We now repeat this trick (without gaussian weigting) 
% for model/tumor/section/slide: scales which show a hierarchy but lack spatial relationships. 

numberOfModels=length(params.samples.modelNames);
meanGroupMats=cell(numberOfMarkers,numberOfRegions);
deltaGroupMats=cell(numberOfMarkers,numberOfRegions);
for regionNumber=1:numberOfRegions
    parfor markerNumber=1:numberOfMarkers
        
        meanGroupMats{markerNumber,regionNumber}=zeros(size(annoMats{markerNumber,regionNumber}));
        %deltaGroupMats{markerNumber,regionNumber}=zeros(size(annoMats{markerNumber,regionNumber}));
        for modelNumber=1:numberOfModels
            isInModel=annoMats{markerNumber,regionNumber}(:,1)==modelNumber;% find pixels in model
            meanInModel=mean(trueMats{markerNumber,regionNumber}(isInModel)); % calculate mean intensity of pixels in model
            %represent first column for each pixel by the mean intensity across all pixels in the model
            meanGroupMats{markerNumber,regionNumber}(isInModel,1)=meanInModel; 
            for tumorNumber=1:3
                isInTumor=isInModel&annoMats{markerNumber,regionNumber}(:,2)==tumorNumber; % find pixels in tumor (&model)
                meanInTumor=mean(trueMats{markerNumber,regionNumber}(isInTumor));% calculate mean intensity of pixels in tumor
                %represent second column for each pixel by the mean intensity across all pixels in the model
                meanGroupMats{markerNumber,regionNumber}(isInTumor,2)=meanInTumor; % represent 
                for sectorNumber=1:3
                    isInSector=isInTumor&annoMats{markerNumber,regionNumber}(:,3)==sectorNumber;% find pixels in sector
                    meanInSector=mean(trueMats{markerNumber,regionNumber}(isInSector)); %calculate mean intensity of pixels in sector
                    %represent third column for each pixel by the mean intensity across all pixels in the sector
                    meanGroupMats{markerNumber,regionNumber}(isInSector,3)=meanInSector;
                    for repSlideNumber=1:3
                        isInSlide=isInSector&annoMats{markerNumber,regionNumber}(:,4)==repSlideNumber;
                        %represent fourth column for each pixel by the mean
                        %intensity across all pixels in the slide
                        meanInSlide=mean(trueMats{markerNumber,regionNumber}(isInSlide));
                        meanGroupMats{markerNumber,regionNumber}(isInSlide,4)=meanInSlide;
                    end
                end
                
            end
            
        end
        % We now sequentially start from the highest length scale, and
        % subtract out its average contribution and pass the residual onto
        % the lower length scales.
        deltaTemp=zeros(size(annoMats{markerNumber,regionNumber}));
        deltaTemp(:,1)=meanGroupMats{markerNumber,regionNumber}(:,1)-mean(meanGroupMats{markerNumber,regionNumber}(:,1));
        for scale=2:size(annoMats{markerNumber,regionNumber},2)
            deltaTemp(:,scale)=...
                meanGroupMats{markerNumber,regionNumber}(:,scale)-...
                meanGroupMats{markerNumber,regionNumber}(:,scale-1);
        end
         % Combine spatial and non-spatial scale downsampling
        deltaGroupMats{markerNumber,regionNumber}=horzcat(deltaTemp,...
            downsampleMats{markerNumber,regionNumber}(:,2:end));
    end
end
%% Break down variances 
% Combine both spatial and non-spatial scales in break down of variance.
% Overall variance over scales (Across all models/samples)
varExp=zeros(numberOfMarkers,8,numberOfRegions);
for regionNumber=1:numberOfRegions
    for markerNumber=1:numberOfMarkers
        varExp(markerNumber,:,regionNumber)=...
            100*var(deltaGroupMats{markerNumber,regionNumber})/var(trueMats{markerNumber,regionNumber});
    end
end

%% Fig 5b table: Plot %var across scales
showLabels=true;
markersToShow=1:6;

markerLabels=vertcat(cellfun(@(x) x(markersToProfile), params.microscopy.markerSets,'Unif',false));
markerLabels=horzcat(markerLabels{:});
scaleLabels={'Model','Tumor','Sector','Image Region','Microenvironmental','Cellular','Subcellular'};
scalesToShow=find(~strcmp(scaleLabels,'Section'));
regionNames={'Nuclear','Cytoplasmic'};

regionNumber=1; % Only nuclear
figure;
%subplot(1,2,regionNumber);
imagesc(varExp(markersToShow,scalesToShow,regionNumber),[0,50]);
colormap('gray');
if(showLabels)
    cbar=colorbar;
    cbar.Label.String='% Variance Explained';

    set(gca,'YTickLabel',markerLabels,'XTickLabel',scaleLabels,'XTickLabelRotation',30);
    title(['%' regionNames{regionNumber} ' Intensity Variance From Different Length Scales']);
else
    set(gca,'XTickLabel',[],'YTickLabel',[]);
end


%% Fig 5b: Tumor Plots to show around table
sampleMapReordered=params.samples.sampleMap(:,:,params.figs.modelOrder);
sampleMapReordered=sampleMapReordered(:,params.figs.repTumorOrder,:);

ifIntensities=load(params.microscopy.avgValueFile);
tumorIntensities=ifIntensities.meanTumorIntensity;
bcatIf=squeeze(nanmean(tumorIntensities(:,1,2,:),4));
perkIf=squeeze(nanmean(tumorIntensities(:,2,3,:),4));
bcatIf=(bcatIf-mean(bcatIf))/std(bcatIf);
perkIf=(perkIf-mean(perkIf))/std(perkIf);



f=figure;
TumorPlot(bcatIf(sampleMapReordered),-1,1,gray,[]);
f=figure;
TumorPlot(perkIf(sampleMapReordered),-1,1,gray,[]);

%% Fig 5b: Sample Whole Tumor Section Images to show around table

%scaledImageList=load(params.microscopy.bgSubtractedImgList); 
scaledImageList=LoadScaledImgList();
scaledImageList=scaledImageList.scaledImageList;

modelNumber=params.figs.modelOrder(2); %second row
tumorNumber=2;%middle tumor

regionNames={'Dorsal','Central','Ventral'};
markerNameList={'ki67','pS6','bCat'};
sampleNum=33;
for markerCounter=1:length(markerNameList)
    markerName=markerNameList{markerCounter};
    switch(markerName)

        case 'ki67'
            markerSet=2;
            markerScales=[0,30000;10000,25000;0,40000;0,20000];
            channelsToShow=[1,2,4];
            repNumber=2;
            
        case 'pS6'
            markerSet=1;
            markerScales=[0,30000;0,30000;1000,40000;0,20000];
            channelsToShow=[1,3,4];
            repNumber=2;
            
        case 'bCat'
            markerSet=1;
            markerScales=[0,30000;1000,40000;0,40000;0,20000];
            channelsToShow=[1,2,4];
            repNumber=2;
            
    end
    
    tImg=scaledImageList{markerSet,sampleNum}{repNumber};
    info=GetImagingInfo(tImg.imgFileMat{1});
    mpp=info.MPP*mean(tImg.dimensions(1,:)./tImg.dimensions(3,:));
    [~,pixelRegion]=Identify_Tumor_Mask(tImg);
    img=tImg.LoadImage('MagLevel',3,'PixelRegion',pixelRegion);
    f=figure;
    Display_Image_With_ScaleBar(img(:,:,channelsToShow),markerScales(channelsToShow,:),mpp,...
        'scaleBarInMicrons',600,'scaleBarYPos',0.95,'scaleBarXPos',0.05,'showScaleBarText',false,...
        'Colors',{'B','R','K'});
    axis off equal;
    
end
%% Fig 5b: Ki67 Zoom image
markerSet=2;

channelsToShow=[1,2,4];
repNumber=2;
tImg=scaledImageList{markerSet,sampleNum}{repNumber};
info=GetImagingInfo(tImg.imgFileMat{1});
mpp=info.MPP*mean(tImg.dimensions(1,:)./tImg.dimensions(3,:));
[~,pixelRegion]=Identify_Tumor_Mask(tImg);
img=tImg.LoadImage('MagLevel',1,'PixelRegion',{16*pixelRegion{1},16*pixelRegion{2}});


f=figure;
markerScales=[0,30000;1000,25000;0,40000;0,20000];
rowRange=[10000:10500];
colRange=[11000:11500];
Display_Image_With_ScaleBar(img(rowRange,colRange,channelsToShow),markerScales(channelsToShow,:),mpp/16,...
    'scaleBarInMicrons',20,'scaleBarYPos',0.95,'scaleBarXPos',0.05,'showScaleBarText',false,...
    'Colors',{'B','R','K'});
axis off equal;

%% Fig. 5a: Downsampling proces on a single section
markerName='pS6';
sampleNumber=33;
repNumber=2;
maskType='nuclearMask';
downsampleSigmas=[Inf,1000,100,10];

markerNum=find(strcmpi(markerLabels,markerName),1);
markerSet=ceil(markerNum/markersPerSet);
markerInSet=markersToProfile(rem(markerNum-1,markersPerSet)+1);

tImg=imgList{markerSet,sampleNumber}{repNumber};
[img,masks]=Get_Image_And_Masks(tImg); 
maskToUse=masks.(maskType);


[downsamples,~]=Hierarchical_DownSampling(img(:,:,markerInSet),maskToUse,downsampleSigmas);

areaToShow='full';
switch(areaToShow)
    case 'full'
        rowRange=1:size(img,1);
        colRange=1:size(img,2);
    case 'region'
        rowRange=5500:6500;
        colRange=6500:7500;
end
responseScalingRange=[-5000,5000];
figure;
ax=gca;
imagesc(img(:,:,markerInSet),responseScalingRange);
colormap(ax,redbluecmap);
%Display_Image(img(:,:,[1,markerInSet,4]),ax,repmat([0,35000],3,1),{'B','R','K'},[]);
%title('Raw Image');
ax.XTick=[];
ax.YTick=[];
axis square
for i=1:length(downsamples)
    figure;
    ax=gca;
    %ax=subplot(2,3,i+1);
    imagesc(downsamples{i}(rowRange,colRange),responseScalingRange);
    colormap(ax,redbluecmap);
    %colorbar;
    %title(scaleLabels(i+3));
    ax.XTick=[];
    ax.YTick=[];
end
