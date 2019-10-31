function axisHandles=DataHeatMapPlot(dataStruct,varargin)

fields=fieldnames(dataStruct);
numberOfFields=length(fields);
numberOfRows=size(dataStruct.(fields{1}).data,1);

p=inputParser;

validateStruct=@(y) all(structfun(@(x) all(ismember(...
    {'data','columnLabels'},fieldnames(x))),y));
areNRowsEqual=@(y) all(structfun(@(x) size(x.data,1),y)==numberOfRows);
addRequired(p,'dataStruct',@(x) isstruct(x) & validateStruct(x) & areNRowsEqual(x));

isOrderMat=@(x,numEl) validateattributes(x,...
    {'numeric'},{'vector','positive','integer','<=',numEl});
addParameter(p,'rowOrder',1:numberOfRows, @(x) isOrderMat(x,numberOfRows));

defaultCmaps=cell(numberOfFields,1);
defaultCmaps(:)={bone};
validateColormaps=@(c) iscell(c)&& numel(c)==numberOfFields && ...
    all(cellfun(@(x) isValidColormap(x),c));
addParameter(p,'colormaps',defaultCmaps,@(x) validateColormaps(x));

defaultFontSize=8;
addParameter(p,'fontSize',defaultFontSize,@(x) validateattributes(x,...
     {'numeric'},{'scalar','positive'}));

defaultColOrders=struct2cell(structfun(@(x) 1:size(x.data,2),dataStruct,'Unif',false));
%maxEls=cellfun(@max,defaultColOrders);
%testOrders=@(x) all(arrayfun(@(y) isOrderMat(x{y},maxEls(y)),1:length(x)));
addParameter(p,'colOrders',defaultColOrders, @(x) iscell(x)); %need to add test that it is a good order matrix

defaultPlotWidths=structfun(@(x) size(x.data,2),dataStruct);
addParameter(p,'plotWidths',defaultPlotWidths,...
    @(x) validateattributes(x,{'numeric'},{'vector','positive'}));

defaultRanges=cell(numberOfFields,1);
addParameter(p,'dataRanges',defaultRanges,...
    @(x) iscell(x) && numel(x) == numberOfFields);

parse(p,dataStruct,varargin{:});
rowOrder=p.Results.rowOrder;
colOrders=p.Results.colOrders;
plotWidths=p.Results.plotWidths;

dataMats=cell(numberOfFields,1);
columnLabels=cell(numberOfFields,1);
titles=fields;
%colormaps=cell(numberOfFields,1);
%colormaps(:)={bone};
colormaps=p.Results.colormaps;

dataRanges=cell(numberOfFields,1);
for fieldCounter=1:numberOfFields
    fieldName=fields{fieldCounter};
    temp=dataStruct.(fieldName).data(rowOrder,:);
    dataMats{fieldCounter}=temp(:,colOrders{fieldCounter});
    if(~isempty(dataStruct.(fieldName).columnLabels))
        columnLabels{fieldCounter}=dataStruct.(fieldName).columnLabels(colOrders{fieldCounter});
    end
    if(isempty(p.Results.dataRanges{fieldCounter}))
        dataRanges{fieldCounter}= [nanmin(dataMats{fieldCounter}(:)),nanmax(dataMats{fieldCounter}(:))];
    else
        dataRanges{fieldCounter}=p.Results.dataRanges{fieldCounter};
    end
end


%dataRanges=cellfun(@(x) [nanmin(x(:)),nanmax(x(:))],dataMats,'Unif',false);
colLabelAngle=-90*ones(numberOfFields,1);
%colLabelAngle(end)=90;
axisHandles=AlignedHeatmaps(dataMats,'colormaps',colormaps,...
    'colLabels',columnLabels,'titles',titles,'plotWidths',plotWidths,...
    'dataRanges',dataRanges,'colLabelAngle',colLabelAngle,'FontSize',...
    p.Results.fontSize);

end

function flag=isValidColormap(cmap)
flag=size(cmap,2)==3&all(cmap(:)>=0)&all(cmap(:)<=1);

end
