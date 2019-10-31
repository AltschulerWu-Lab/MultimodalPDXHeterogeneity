%% Script to calculate the average intensity of different markers
% This script loops over all 36 samples and calculates the average
% intensities for each of the 4 markers in the 4 marker sets in each of the
% 3 replicate sections. In each section mean intensities are calculated in
% the overall tissue region, and stromal, tumor(=tissue-stromal) compartments and nuclear.
% The results are stored in the file pointed to by params.microscopy.avgValueFile
%%
addpath('../');
addpath('../Growth_Curves/');
addpath('../Pathology/');
addpath(genpath('../Common/Plotting/'));
params=GetParams();
sampleMap=params.samples.sampleMap;
[numberOfRegions,numberOfTumorsPerModel,numberOfModels]=size(sampleMap);
numberOfMarkerSets=length(params.microscopy.markerSets);
numberOfMarkersPerSet=length(params.microscopy.markerSets{1});
numberOfReplicateSections=3;
numberOfSamples=params.samples.numberOfSamples;
%% Set up full list of images
useBgSubtraction=true;
tImgList=cell(numberOfSamples,numberOfMarkerSets,numberOfReplicateSections);
if(useBgSubtraction)
    %temp=load(params.microscopy.bgSubtractedImgList);
    temp=LoadScaledImgList();
    for markerSet=1:numberOfMarkerSets
        
        for sampleCounter=1:numberOfSamples
           
            for replicateCounter=1:numberOfReplicateSections
                tImgList{sampleCounter,markerSet,replicateCounter}=...
                    temp.scaledImageList{markerSet,sampleCounter}{replicateCounter};
            end
            
        end
       
    end
    outFile=params.microscopy.avgValueFile;
else
    for markerSet=1:numberOfMarkerSets
        tic;
        for sampleCounter=1:numberOfSamples
            afiFile=GetImages('IF','sampleNumber',sampleCounter,...
                'markerSet',markerSet);
            for replicateCounter=1:min(numberOfReplicateSections,length(afiFile))
                tImgList{sampleCounter,markerSet,replicateCounter}=...
                    Afi2TImg(afiFile{replicateCounter});
            end
            
        end
        toc;
    end
    outFile=params.microscopy.avgValueRawFile;
end


%% Image Intensity Calculations

meanSampleIntensity=NaN*ones(numberOfSamples,numberOfMarkerSets,...
    numberOfMarkersPerSet, numberOfReplicateSections);
meanTumorIntensity=NaN*ones(numberOfSamples,numberOfMarkerSets,...
    numberOfMarkersPerSet, numberOfReplicateSections);
meanStromalIntensity=NaN*ones(numberOfSamples,numberOfMarkerSets,...
    numberOfMarkersPerSet, numberOfReplicateSections);
meanMaskIntensity=NaN*ones(numberOfSamples,numberOfMarkerSets,...
    numberOfMarkersPerSet, numberOfReplicateSections);
pctNucPositive=NaN*ones(numberOfSamples,numberOfMarkerSets,...
    numberOfMarkersPerSet, numberOfReplicateSections);
pctPositive=NaN*ones(numberOfSamples,numberOfMarkerSets,...
    numberOfMarkersPerSet, numberOfReplicateSections);
for markerSet=1:numberOfMarkerSets
    for sampleCounter=1:numberOfSamples
        tic;
        for replicateCounter=1:numberOfReplicateSections
            tImg=tImgList{sampleCounter,markerSet,replicateCounter};
            if(~isempty(tImg))
                [tumorMask,pixelRegion]=Identify_Tumor_Mask(tImg);
                tumorMaskCropped=tumorMask(pixelRegion{1}(1):pixelRegion{1}(2),...
                    pixelRegion{2}(1): pixelRegion{2}(2));
                img=tImg.LoadImage('MagLevel',3,'pixelRegion',pixelRegion);
                %figure;
                %imagesc(img(:,:,4))
                dapiImg=img(:,:,1);dapiImg=dapiImg/max(dapiImg(:));
                dapiMask=imbinarize(dapiImg,adaptthresh(dapiImg,0.5))&tumorMaskCropped;
                vimImg=img(:,:,4);vimImg=vimImg/max(vimImg(:));
                stromalMask=imbinarize(vimImg,adaptthresh(vimImg,0.5))&tumorMaskCropped;
                closeRadius=2;
                se = strel('disk',closeRadius);
                cellMask=imclose(dapiMask&~stromalMask,se);
                tumorMask=cellMask&~stromalMask;
                
                for markerCounter=1:numberOfMarkersPerSet
                    markerIntensity=img(:,:,markerCounter);
                    mImg=markerIntensity/max(markerIntensity(:));
                    markerHighMask=imbinarize(mImg,adaptthresh(mImg,0.5))&tumorMaskCropped;
                    
                    
                    
                    meanSampleIntensity(sampleCounter,markerSet,markerCounter,replicateCounter)=...
                        mean(markerIntensity(cellMask));
                    meanTumorIntensity(sampleCounter,markerSet,markerCounter,replicateCounter)=...
                        mean(markerIntensity(tumorMask));
                    meanStromalIntensity(sampleCounter,markerSet,markerCounter,replicateCounter)=...
                        mean(markerIntensity(stromalMask));
                    meanMaskIntensity(sampleCounter,markerSet,markerCounter,replicateCounter)=...
                        mean(markerIntensity(tumorMaskCropped));
                    pctPositive(sampleCounter,markerSet,markerCounter,replicateCounter)=...
                        nnz(markerHighMask&cellMask)/nnz(cellMask);
                    pctNucPositive(sampleCounter,markerSet,markerCounter,replicateCounter)=...
                        nnz(markerHighMask&dapiMask)/nnz(dapiMask);
                end
                %         subplot(1,2,1);
                %         imagesc(ps6Img);
                %         axis off equal;
                %         subplot(1,2,2);
                %         imagesc(ps6Mask&~stromalMask);
                %         axis off equal;
                %         subplot(3,3,(r-1)*3+c);
                %         imagesc(ps6Mask&~stromalMask);
                %         axis off equal;
            end
        end
        toc;
    end
    
end


%%


save(outFile,'meanSampleIntensity',...
   'meanTumorIntensity','meanStromalIntensity','meanMaskIntensity',...
    'pctNucPositive','pctPositive');

%%
markerSet=1;
markerNumber=3;
avgVal=nanmean(meanTumorIntensity,4);
markerVals=squeeze(avgVal(:,markerSet,markerNumber));
TumorPlot(markerVals(sampleMap),min(markerVals),max(markerVals),redbluecmap,[]);
suptitle_mod(params.microscopy.markerSets{markerSet}{markerNumber},18);