function axisHandles = AlignedHeatmaps(dataMats,varargin)

    p=inputParser;
    
    % Data Matrices to be plotted organized as a cell array
    validateData=@(x) iscell(x)& (range(cellfun(@(y) size(y,1),x))==0)&...
       all(cellfun(@(y) isnumeric(y)||islogical(y),x)); 
    % Confirm that dataMats is a cell array containing numeric matrices with 
    % the same number of rows
    addRequired(p,'dataMats',@(x) validateData(x));
    numberOfDataMats=numel(dataMats);
 
    defaultPlotWidths=cellfun(@(x) size(x,2),dataMats);
    addParameter(p,'plotWidths',defaultPlotWidths,...
        @(x) validateattributes(x,{'numeric'},{'vector','positive'}));
    
    defaultTitles=cell(numberOfDataMats,1);
    addParameter(p,'titles',defaultTitles,...
        @(x) iscell(x) && numel(x) == numberOfDataMats);
    
    defaultColLabels=cell(numberOfDataMats,1);
    addParameter(p,'colLabels',defaultColLabels,...
        @iscell);
    
    defaultRanges=cell(numberOfDataMats,1);
    addParameter(p,'dataRanges',defaultRanges,...
        @(x) iscell(x) && numel(x) == numberOfDataMats); %should add check 
    % for 2 monotonically increasing elements (when not empty)
    
    defaultCmaps=cell(numberOfDataMats,1);
    defaultCmaps(:)={bone};
    validateColormaps=@(c) iscell(c)&& numel(c)==numberOfDataMats && ...
        all(cellfun(@(x) isValidColormap(x),c));
    addParameter(p,'colormaps',defaultCmaps,@(x) validateColormaps(x));

    defaultColLabelAngle=-90*ones(numberOfDataMats,1);
    addParameter(p,'colLabelAngle',defaultColLabelAngle,...
        @(x) validateattributes(x,{'numeric'},{'vector','numel',numberOfDataMats}));
    
    defaultColLabelPos=cell(numberOfDataMats,1);
    defaultColLabelPos(:)={'bottom'};
    addParameter(p,'colLabelPos',defaultColLabelPos,...
        @(x) iscell(x) && numel(x)==numberOfDataMats&& ...
        all(ismember(x,{'top','bottom'})));
    
    defaultGridColor='r';
    addParameter(p,'gridColor',defaultGridColor); 

    
    defaultFontSize=12;
    addParameter(p,'fontSize',defaultFontSize,...
        @(x) isscalar(x) && isnumeric(x) && x>0);
    
    
    defaultFontStyle='Arial';
    addParameter(p,'fontStyle',defaultFontStyle);
    
   
  
    parse(p,dataMats,varargin{:});
    
    cumSubplotWidth=cumsum([1;p.Results.plotWidths(:)]);
    totalSubplotWidth=sum(p.Results.plotWidths);
    %axisHandles=struct;%zeros(numberOfDataMats,1);
    for plotCounter=1:numberOfDataMats
        
        subplotRange=cumSubplotWidth(plotCounter);
        if(p.Results.plotWidths(plotCounter)>1)          
             subplotRange= subplotRange+[0,p.Results.plotWidths(plotCounter)-1];
        end
        nanColor=[0.2,0,0];
        axisHandles(plotCounter)=subplot(1,totalSubplotWidth,subplotRange);
        if(isempty(p.Results.dataRanges{plotCounter}))
            imagescNaN(dataMats{plotCounter},[],p.Results.colormaps{plotCounter},...
                nanColor);
            %imagesc(dataMats{plotCounter});
            %pcolor(double(dataMats{plotCounter}));
        else
            %imagesc(dataMats{plotCounter},p.Results.dataRanges{plotCounter});
            imagescNaN(dataMats{plotCounter},p.Results.dataRanges{plotCounter},...
                p.Results.colormaps{plotCounter},nanColor);
            %set(h,'alphadata',~isnan(dataMats{plotCounter}));
            %pcolor(double(dataMats{plotCounter}));
            %set(gca,'CLim',p.Results.dataRanges{plotCounter});
        end
        set(gca,'YTickLabel',[],'TickLength',[0,0]);
        %colormap( axisHandles(plotCounter),p.Results.colormaps{plotCounter});
        
        xTick= (0:size(dataMats{plotCounter},2))+0.5;
        yTick=(0:size(dataMats{plotCounter},1))+0.5;
        if(length(xTick)>100)
          xTick=[];               
        end
        if(length(yTick)>100)
          yTick=[];               
        end
        
        set(gca,'XTick',xTick,'YTick',yTick,...
            'XTickLabel',[],'YTickLabel',[],'GridColor','r');
        grid on;
        
        if(~isempty(p.Results.colLabels{plotCounter}))
            for j=1:numel(p.Results.colLabels{plotCounter})
                if(strcmp(p.Results.colLabelPos{plotCounter},'top'))
                    ypos= 0.5;
                    valign='middle';
                else
                    ypos= size(dataMats{plotCounter},1)-0.5*...
                        sind(p.Results.colLabelAngle(plotCounter));
                    valign='middle';
                end
                text(j,...
                    ypos,p.Results.colLabels{plotCounter}{j},...
                    'VerticalAlignment',valign,'FontSize',p.Results.fontSize,...
                    'FontName',p.Results.fontStyle,'HorizontalAlignment','left',...
                    'Rotation',p.Results.colLabelAngle(plotCounter),'Interpreter','none');
            end
        end
        if(~isempty(p.Results.titles{plotCounter}))
            title(p.Results.titles{plotCounter})
        end
    
        
    end
end

function flag=isValidColormap(cmap)
   flag=size(cmap,2)==3&all(cmap(:)>=0)&all(cmap(:)<=1);
    
end

