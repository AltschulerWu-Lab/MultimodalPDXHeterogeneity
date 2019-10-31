classdef TissueImageCorrected

% Convenience class that provides a standard interface to open background-corrected tissue images stored as svs files. Parameters specifying the shading correction model as determined by Build_BG_Model must be passed at initiation time. 
    properties
        bgInterp;
        shadeInterp;
	intensityMultiplier;
	intensityOffset;
        numberOfChannels;
        numberOfMagLevels;
        dimensions;
        imgFileMat;
        imgLayerMat;
        micronsPerPixel;
        channelNames;
        bitDepth;
        
    end
    
    methods
        function obj=TissueImageCorrected(tImg,bgInterp,shadeInterp)
            obj.imgFileMat=tImg.imgFileMat;
            obj.imgLayerMat=tImg.imgLayerMat;
            obj.micronsPerPixel=tImg.micronsPerPixel;
            obj.channelNames=tImg.channelNames;
            obj.bitDepth=tImg.bitDepth;
            obj.dimensions=tImg.dimensions;
            [obj.numberOfMagLevels,obj.numberOfChannels]=size(tImg.imgFileMat);
            obj.intensityMultiplier=ones(obj.numberOfChannels,1);
	    obj.intensityOffset=zeros(obj.numberOfChannels,1);			
            
            if(~iscell(bgInterp)||~isequal(size(bgInterp),size(tImg.imgFileMat))||...
                    ~all(all(cellfun(@(x) isa(x,'griddedInterpolant'),bgInterp))))
                error('Invalid bgInterp');
            else
                obj.bgInterp=bgInterp;
            end
            
            if(~iscell(shadeInterp)||~isequal(size(shadeInterp),size(tImg.imgFileMat))||...
                    ~all(all(cellfun(@(x) isa(x,'griddedInterpolant'),shadeInterp))))
                error('Invalid shadeInterp');
            else
                obj.shadeInterp=shadeInterp;
            end
            
            
        end
        
        
        function img=LoadImage(obj,varargin)
            p=inputParser;
            
            
            addParameter(p,'MagLevel',1,@(x) validateattributes(x,{'numeric'},...
                {'scalar','integer','positive','<=',obj.numberOfMagLevels}));
            addParameter(p,'Channels',1:obj.numberOfChannels,@(x) ...
                validateattributes(x,{'numeric'},{'vector','integer',...
                'positive','<=',obj.numberOfChannels}));
            addParameter(p,'PixelRegion',{}, @(x) validateattributes(x,...
                {'cell'},{'numel',2}));
            
            
            parse(p,varargin{:});
            
            %add checker function for PixelRegion
            if(isempty(p.Results.PixelRegion))
                img=zeros([obj.dimensions(p.Results.MagLevel,:),length(p.Results.Channels)]);
                pixelRegion={[1,obj.dimensions(p.Results.MagLevel,1)],...
                    [1,obj.dimensions(p.Results.MagLevel,2)]};
            else
                img=zeros([cellfun(@diff,p.Results.PixelRegion)+1,length(p.Results.Channels)]);
                pixelRegion=p.Results.PixelRegion;
            end
            for channelCounter=1:length(p.Results.Channels)
                
                channelNumber=p.Results.Channels(channelCounter);
                img(:,:,channelCounter)=imread(obj.imgFileMat{p.Results.MagLevel,channelNumber},...
                    obj.imgLayerMat(p.Results.MagLevel,channelNumber),'PixelRegion',p.Results.PixelRegion);
                
            end
            
            cleanImg=img;
            [X,Y]=meshgrid(pixelRegion{1}(1):pixelRegion{1}(2),...
               pixelRegion{2}(1):pixelRegion{2}(2));
%             [Y1,X1]=ndgrid(pixelRegion{1}(2):-1:pixelRegion{1}(1),...
%                 pixelRegion{2}(2):-1:pixelRegion{2}(1));
             %[Y,X]=ndgrid(pixelRegion{2}(1):pixelRegion{2}(2),...
             %    pixelRegion{1}(1):pixelRegion{1}(2));
            magLevel=p.Results.MagLevel;
            for channelCounter=1:length(p.Results.Channels)
                channelNumber=p.Results.Channels(channelCounter);
     
                
                temp=bsxfun(@rdivide,cleanImg(:,:,channelCounter),...
                    obj.shadeInterp{magLevel,channelNumber}(pixelRegion{2}(1):pixelRegion{2}(2)))-...
                    obj.bgInterp{magLevel,channelNumber}(X',Y');
                temp=temp*obj.intensityMultiplier(channelCounter)+obj.intensityOffset(channelCounter);
		temp(temp<0)=0;
                cleanImg(:,:,channelCounter)=temp;
            end
            img=cleanImg;
        end
        
        
    end
end
