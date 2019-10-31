function [fileList,fileInfo]=GetImages(imgType,varargin)

% Convenience function to pull filenames for image files given sample requirements
% INPUT
% imgType   -   must be either 'H&E' or 'IF'
%               
% varargin  - There are various optional parameters, allowing specification of sample number, marker set and section number.             
% OUTPUT
% fileList   - a cell array containing the filenames of the afi(IF) or svs (H&E) files matching specified conditions
% fileInfo   - a table with rows corresponding to the files pulled, and columns providing inforamtion such as sampleNumber, sectionNumber, markerSet and so on.
p=inputParser;

validImageTypes={'H&E','IF'};
addRequired(p,'imgType',@(x) ismember(x,validImageTypes));

addParameter(p,'markerSet',1:4,@(x) ...
    validateattributes(x,{'numeric'},{'vector','integer',...
    'positive','<=',4}));

addParameter(p,'sampleNumber',1:36,@(x) ...
    validateattributes(x,{'numeric'},{'vector','integer',...
    'positive','<=',36}));

addParameter(p,'sectionNumber',1:100,@(x) ...
    validateattributes(x,{'numeric'},{'vector','integer',...
    'positive','<=',100}));

% addParameter(p,'repNumber',1:10,@(x) ...
%     validateattributes(x,{'numeric'},{'vector','integer',...
%     'positive','<=',10}));

validDataLocaltions={'Aperio','ForLeidos'};
addParameter(p,'dataLocation','ForLeidos',@(x) ismember(x,validDataLocaltions));


parse(p,imgType,varargin{:});

params=GetParams('microscopy');

switch(imgType)
    case 'H&E'
        heLookup=readtable(params.microscopy.heLookupFile);
        sampleNumber=cellfun(@(x) str2double(x{1}),regexp(heLookup.Slide,'_','split'));
        sectionNumber=cellfun(@(x) str2double(x{2}),regexp(heLookup.Slide,'_','split'));
        switch(p.Results.dataLocation)
            case 'Aperio'
                imgFolder=params.microscopy.aperioImgDir;
                fileName=fullfile(imgFolder,heLookup.Folder,heLookup.Filename);
            case 'ForLeidos'
                imgFolder=fullfile(params.microscopy.ToDepositDir,'Microscopy_HnE');
                fileName=fullfile(imgFolder,heLookup.Filename);
        end
     
        
        info=table(sampleNumber,sectionNumber,fileName);
        isValid=ismember(sampleNumber,p.Results.sampleNumber(:)) & ...
            ismember(sectionNumber,p.Results.sectionNumber(:));
        
        isValid(isValid)=cellfun(@(x) exist(x,'file'),fileName(isValid));
        
        fileList=fileName(isValid);
        fileInfo=info(isValid,:);
        
    case 'IF'
        
        ifLookup=readtable(params.microscopy.ifLookupFile);
        sampleNumber=cellfun(@(x) str2double(x{1}),regexp(ifLookup.Slide,'_','split'));
        sectionNumber=cellfun(@(x) str2double(x{2}),regexp(ifLookup.Slide,'_','split'));
        markerSet=cellfun(@(x) str2double(x{1}),regexp(ifLookup.MarkerSet,'\d$','match'));
        switch(p.Results.dataLocation)
            case 'Aperio'
                imgFolder=fullfile(params.microscopy.aperioImgDir,ifLookup.folder);
            case 'ForLeidos'
                 imgFolder=fullfile(params.microscopy.ToDepositDir,...
                    strcat('Microscopy_IF_Sample', ...
                    arrayfun(@num2str,sampleNumber,'Unif',false), '_MS', ...
                    arrayfun(@num2str,markerSet,'Unif',false)));
               
        end
        fileName=fullfile(imgFolder,ifLookup.filename);
        
        info=table(sampleNumber,sectionNumber,markerSet,ifLookup.comment, ...
            fileName,'VariableNames',{'sampleNumber','sectionNumber',...
            'markerSet','comments','fileName'});
        isValid=ismember(sampleNumber,p.Results.sampleNumber(:)) & ...
            ismember(sectionNumber,p.Results.sectionNumber(:))& ...
            ismember(markerSet,p.Results.markerSet(:)) &  ...
            ~cellfun(@isempty,ifLookup.folder) & ...
            ~cellfun(@isempty,ifLookup.filename);
         
        isValid(isValid)=cellfun(@(x) exist(x,'file'),fileName(isValid));

        
        fileList=fileName(isValid);
        fileInfo=info(isValid,:);
        
        
end



end
