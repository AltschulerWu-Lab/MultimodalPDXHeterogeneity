function TumorPlot(dataVals,minVal,maxVal,cmap,varargin)


p=inputParser;

addOptional(p,'modelNames',[],@(x) isempty(x) | (iscell(x) & numel(x)==4));
addOptional(p,'lineColor','k');
addOptional(p,'nanColor',[0.5,0.5,0.5]);

parse(p,varargin{:});
modelNames=p.Results.modelNames;
lineColor=p.Results.lineColor;
nanColor=p.Results.nanColor;


if(~iscell(dataVals))
    cmapIdx=ValToCmapIdx(dataVals,minVal,maxVal,cmap);
else
    cmapIdx=dataVals;
end


numberOfModels=4;
%numberOfRegions=3;
numberOfReplicateTumors=3;




if(~isempty(modelNames))
    nR=numberOfModels+1;
    nC=numberOfReplicateTumors+1;
    counter=2;
    for replicateCounter=1:numberOfReplicateTumors
        
        subaxis(nR,nC,counter,...
            'SpacingVert',0.01,...
            'SpacingHoriz',0.01,'MR',0,...
            'ML',0,'MT',0,'MB',0);
        %subplot(numberOfModels+1,numberOfReplicateTumors+1,counter);
        
        text(0,0.5,['Replicate Tumor ' num2str(replicateCounter)]);
        axis off equal;
        counter=counter+1;
    end
else
    nR=numberOfModels;
    nC=numberOfReplicateTumors;
    counter=1;
end
for modelCounter=1:numberOfModels
    if(~isempty(modelNames))
        subaxis(nR,nC,counter,...
            'SpacingVert',0.01,...
            'SpacingHoriz',0.01,'MR',0,...
            'ML',0,'MT',0,'MB',0);
        
        %subplot(numberOfModels+1,numberOfReplicateTumors+1,counter);
        text(0,0.5,modelNames{modelCounter});
        axis off equal;
        counter=counter+1;
    end
    for replicateCounter=1:numberOfReplicateTumors
        subaxis(nR,nC,counter,...
            'SpacingVert',0.01,...
            'SpacingHoriz',0.01,'MR',0,...
            'ML',0,'MT',0,'MB',0);
        %subplot(numberOfModels+1,numberOfReplicateTumors+1,counter);
        PlotTumor(squeeze(cmapIdx(:,replicateCounter,modelCounter)),cmap,...
            lineColor,nanColor);
        counter=counter+1;
    end
end




end

function cmapIdx=ValToCmapIdx(valMat,minVal,maxVal,cmap)

cmapLength=length(cmap);
rescaledVals=(valMat-minVal)./(maxVal-minVal);
rescaledVals(rescaledVals<0)=0;
rescaledVals(rescaledVals>1)=1;
cmapIdx=round(rescaledVals*(cmapLength-1))+1;
end

function PlotTumor(colorIdxs,cmap,lineColor,nanColor)
lineWidth=2;
if(~iscell(colorIdxs))
    
    topColorIdx=colorIdxs(1);
    middleColorIdx=colorIdxs(2);
    bottomColorIdx=colorIdxs(3);
    if(~isnan(topColorIdx))
        topColor=cmap(topColorIdx,:);
    else
        topColor=nanColor;
    end
    if(~isnan(middleColorIdx))
        middleColor=cmap(middleColorIdx,:);
    else
        middleColor=nanColor;
    end
    if(~isnan(bottomColorIdx))
        bottomColor=cmap(bottomColorIdx,:);
    else
        bottomColor=nanColor;
    end
    
    
else
    topColor=colorIdxs{1};
    middleColor=colorIdxs{2};
    bottomColor=colorIdxs{3};
    
end

a=2;
b=1;
t = linspace(0,pi,128);
r=a*b./sqrt((b*cos(t)).^2+(a*sin(t)).^2);
x = r.*cos(t);
y = r.*sin(t)+b/2;
patch(x,y,topColor,'LineWidth',lineWidth,'EdgeColor',lineColor);
x=[1 -1 -1 1]*a;
y=[1 1 -1 -1]*b/2;
patch(x,y,middleColor,'LineWidth',lineWidth,'EdgeColor',lineColor);
t = linspace(pi,2*pi,128);
r=a*b./sqrt((b*cos(t)).^2+(a*sin(t)).^2);
x = r.*cos(t);
y = r.*sin(t)-b/2;
patch(x,y,bottomColor,'LineWidth',lineWidth,'EdgeColor',lineColor);
axis off equal;

end
