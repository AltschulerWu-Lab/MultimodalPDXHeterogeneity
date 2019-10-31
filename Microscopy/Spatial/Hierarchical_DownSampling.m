function [downsamples,responses] = Hierarchical_DownSampling(fullResponse,samplingMask,downsampleFactors)

% Core function to break down the variation in intensity in an image into contributions arising 
% from different spatial scales. This is done by convolving the image with a gaussian, and assigning
% the variation in this image as  that corresponding to the spatial scale defined by sigma. The residual
% image (input-convolved) is then convolved with a gaussian of smaller sigma and the process is repeated.
% INPUT
% fullResponse   - A grayscale image, specified in terms of an 2D image array, whose variation we want to analyze.
% samplingMask   - A binary mask of the same dimensions as fullResponse, with true values marking 
% 		   pixels whose intensities we want to use for the variance calculation. For example, 
% 		   a matrix with true values in nuclear pixels would be used to characterize nuclear variation.
% downSampleFactors - The list of decreasing gaussian sigmas (in units of pixels) to be used.
%              
% OUTPUT
% downsamples    - A cell array containing 2D arrays of the same size as fullResponse, containing the variation
% 		   captured at the corresponding length scale
% responses      - A cell array containing 2D arrays of the same size as fullResponse, containing the response
% 		   input to the corresponding length scale. The difference between these is the variation captured at
% 		   the length scale.
    numberOfLevels=length(downsampleFactors);
    downsamples=cell(numberOfLevels+1,1);
    responses=cell(numberOfLevels,1);
    inResponse=fullResponse;
    inResponse(~samplingMask)=0;
    inDensity=samplingMask;
    for level=1:numberOfLevels
        responses{level}=inResponse;
        [levelResponse,inDensity]=Downsample(inResponse,samplingMask,samplingMask,downsampleFactors(level));
        inResponse=inResponse-levelResponse;
        downsamples{level}=levelResponse;
        disp(['Level ' num2str(level) ' done!'])
    end
    downsamples{numberOfLevels+1}=inResponse;
end

function [outResponse,outDensity]=Downsample(inResponse,inDensity,samplingMask,downsampleSigma)
    maxDownsampleSigma=250; %becomes super slow to convolve with large kernel, so we will use a shortcut in this case
    if(islogical(inDensity))
       inDensity=double(inDensity); 
    end
    if(isinf(downsampleSigma))
        outDensity=inDensity;
        outResponse=ones(size(inResponse))*mean(inResponse(samplingMask)); % Re[place with weighted mean based on inDensity?
        outResponse(~samplingMask)=0;
    else
        
        if(downsampleSigma>maxDownsampleSigma)
           
            [outResponse,outDensity]=Downsample(inResponse,inDensity,samplingMask,10);
            
            [outResponse,outDensity]=Downsample(imresize(outResponse,0.1),imresize(outDensity,0.1),imresize(samplingMask,0.1),downsampleSigma/10);
            outResponse=imresize(outResponse,size(inResponse));
            outDensity=imresize(outDensity,size(inResponse));
            
        else
           inResponse(~samplingMask)=0;
           outDensity=imgaussfilt(inDensity,downsampleSigma);
           outResponse=imgaussfilt(inResponse,downsampleSigma)./outDensity;
           outDensity(~samplingMask)=0;
           outResponse(~samplingMask)=0;
        end
        
    end


end
