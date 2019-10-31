function Downsampler(sampleNumber)
% Cluster friendly function to generate downsampling of a single sample
% (out of 36) imaged across 3 replicate sections with 4 marker sets. 
% INPUT
% sampleNum   - An integer between 1 and 36 (potentially supplied by the job-submission script), 
% 		which identifies the sample to be analyzed.
%              
% OUTPUT
% Each run generates downsampling analysis for a single sample across all 4 marker sets
% and 3 replicate sections. The results are stored as .mat files in the folder
% specified by params.microscopy.downsamplingResultsDir in GetParams.m.

cleanImage=true; %USe BG subtraction?

markersToProfile=[2,3]; % 1st and 4th markers i.e. DAPI/vimentin are less interesting
downsampleSigmas=[Inf,1000,100,10];
numberPoints=10000; % Number of randomly sampled points used to calculate variances

numberOfMarkerSets=4;
numberOfReplicateSections=3;
numberOfRegions=2; % nuclear and cytoplasmic


addpath('../../');
addpath('../');
addpath(genpath('../../Common/'));
params=GetParams('microscopy');

if(cleanImage)
    saveDir=params.microscopy.downsamplingResultsDir;
    %temp=load(params.microscopy.bgSubtractedImgList);
    temp=LoadScaledImgList(); %Uncomment the previous line and comment this one if you don't want paths changed related to rootDir
    imgList=temp.scaledImageList;
    
else
    saveDir='/work/bioinformatics/srajar/Leidos/Data/Spatial_Analysis/Raw_Downsampler/';
    temp=load(params.microscopy.rawImgList);
    imgList=temp.rawImageList;
end


for repNumber=1:numberOfReplicateSections
    saveFile=fullfile(saveDir,['downSampler-Sample' num2str(sampleNumber) '-Replicate-' num2str(repNumber) '.mat']);
    if(exist(saveFile,'file'))
        data=load(saveFile);
        downsampleData=data.downsampleData;
        trueData=data.trueData;
        isProfiled=~cellfun(@isempty,data.downsampleData);
        disp(['Sample ' num2str(sampleNumber) ' Rep Number' num2str(repNumber) ' Exists']);
    else
        downsampleData=cell(numberOfMarkerSets,length(markersToProfile),numberOfRegions);
        trueData=cell(numberOfMarkerSets,length(markersToProfile),numberOfRegions);
        isProfiled=false(size(downsampleData));
        disp(['Sample ' num2str(sampleNumber) ' Rep Number' num2str(repNumber) ' Does not Exists']);
    end
    for markerSet=1:numberOfMarkerSets
        tic;
        
        
        try
            tImg=imgList{markerSet,sampleNumber}{repNumber};
            
            if(any(any(~isProfiled(markerSet,:,:))))
                
                [img,masks]=Get_Image_And_Masks(tImg); %Local function defined below
                
                
                %Perform downsampling
                for markerCounter=1:length(markersToProfile)
                    
                    % Load marker response, and scale it to [0,1]
                    response=img(:,:,markersToProfile(markerCounter));
                    %response=response/max(max(response));
                    
                    % Perform downsampling
                    for regionCounter=1:numberOfRegions % nuclear and cytoplasmic
                        
                        if(isempty(downsampleData{markerSet,markerCounter,regionCounter}))
                            if(regionCounter==1)
                                maskToUse=masks.nuclearMask;
                            else
                                maskToUse=masks.cytoplasmicMask;
                            end
                            
                            % Perform downsampling calculations
                            [downsamples,~]=Hierarchical_DownSampling(response,maskToUse,downsampleSigmas);
                            
                            
                            % Calculate variances from randomly samples points
                            randomIdx=randsample(find(maskToUse),numberPoints);
                            reconSample=zeros(numberPoints,length(downsamples));
                            trueSample=response(randomIdx);
                            for i=1:length(downsamples)
                                reconSample(:,i)=downsamples{i}(randomIdx);
                            end
                            %variances{regionCounter}((markerSet-1)*length(markersToProfile)+markerCounter,:)=var(reconSample)./var(trueSample);
                            downsampleData{markerSet,markerCounter,regionCounter}=reconSample;
                            trueData{markerSet,markerCounter,regionCounter}=trueSample;
                            
                            save(saveFile,'trueData','downsampleData');
                            disp(['Sample ' num2str(sampleNumber) ' Marker ' num2str(markersToProfile(markerCounter)) ' x ' num2str(regionCounter) ' done!']);
                        else
                            disp(['Sample ' num2str(sampleNumber) ' Marker ' num2str(markersToProfile(markerCounter)) ' x ' num2str(regionCounter) ' already done. Skipping!']);
                        end
                    end
                    
                end
                            
                disp(['Sample ' num2str(sampleNumber) ' Marker Set ' num2str(markerSet) ' done!']);
            
            else
                disp(['Sample ' num2str(sampleNumber) ' Marker ' num2str(markerSet) ' already done. Skipping!']);
            end

        catch
            disp(['Sample ' num2str(sampleNumber) ' Marker Set ' num2str(markerSet) ' Error!']);
        end
        toc;
    end
end


end

