function final_image=Display_Image_With_ScaleBar(img,color_scaling,micronsPerPixel,varargin)
% Convenience function to display an image with scale bar, while incorporating
% appropriate color handling of different channels.
p=inputParser;

addOptional(p,'axisHandle',gca,@(x) ishandle(x));

defaultColors={'B','G','R','K','Y','M'};
validColors= {'b','B','Blue','g','G','Green','r','R','Red','k',...
    'K','Gray','o','O','Orange','y','Y','Yellow','c','C','Cyan',...
    'm','M','Magenta','','-' ,'None'};
addOptional(p,'colors',defaultColors,@(x) iscell(x) & all(ismember(x,validColors)));

addOptional(p,'mask',[],@(x) isempty(x) || ...
    (islogical(x) && size(x,1)==size(img,1) && size(x,2)==size(img,2)));

addOptional(p,'scaleBarInMicrons',[],@(x) isnumeric(x)&isscalar(x));

defaultScaleBarYPos=0.1;
addOptional(p,'scaleBarYPos',defaultScaleBarYPos,@(x) validateattributes(x,{'numeric'},...
    {'scalar','<=',1,'>=',0}));
defaultScaleBarXPos=0.1;
addOptional(p,'scaleBarXPos',defaultScaleBarXPos,@(x) validateattributes(x,{'numeric'},...
    {'scalar','<=',1,'>=',0}));

defaultShowScaleBarText=true;
addOptional(p,'showScaleBarText',defaultShowScaleBarText,@(x) islogical(x));

parse(p,varargin{:});
axisHandle=p.Results.axisHandle;
colors=p.Results.colors;
mask=p.Results.mask;
scaleBarInMicrons=p.Results.scaleBarInMicrons;
showScaleBarText=p.Results.showScaleBarText;
scaleBarYPos=p.Results.scaleBarYPos;
scaleBarXPos=p.Results.scaleBarXPos;

img=double(img);
[xres,yres,nch]=size(img);
final_image=zeros(xres,yres,3);

for channel=1:nch
    img(:,:,channel)=max(min((img(:,:,channel)-color_scaling(channel,1))...
        /(color_scaling(channel,2)-color_scaling(channel,1)),1),0);
    switch colors{channel}
        case {'b','B','Blue'},    ch = 3;     wt = 1;
        case {'g','G','Green'},   ch = 2;     wt = 1;
        case {'r','R','Red'},     ch = 1;     wt = 1;
        case {'k','K','Gray'},    ch = 1:3;   wt = [1 1 1];
        case {'o','O','Orange'},  ch = 1:2;   wt = [1 .75];
        case {'y','Y','Yellow'},  ch = 1:2;   wt = [.5 .5];
        case {'c','C','Cyan'},    ch = 2:3;   wt = [1 1];
        case {'m','M','Magenta'}, ch = [1 3]; wt = [1 1];
        case {'','-' ,'None'},    ch=[];      wt=[];
    end
    for h=1:length(ch)
        if(isempty(mask))
            final_image(:,:,ch(h))=final_image(:,:,ch(h))+wt(h)*img(:,:,channel);
        else
            final_image(:,:,ch(h))=final_image(:,:,ch(h))+wt(h)*img(:,:,channel)+0.25*mask;
        end
    end
    
end

final_image=min(max(final_image,0),1);
image(final_image,'Parent',axisHandle);

if(isempty(scaleBarInMicrons))
    yBy10InMicrons=(yres*micronsPerPixel/10);
    micronsExp=floor(log10(yBy10InMicrons));
    yBy10InMicrons=round(yBy10InMicrons*10^(-micronsExp))*10^(micronsExp);
    scaleBarLength=round(yBy10InMicrons/micronsPerPixel);
    scaleBarInMicrons=yBy10InMicrons;
else
    scaleBarLength=round(scaleBarInMicrons/micronsPerPixel);
end

hold on;

plot(round(yres*scaleBarXPos)+ [1,scaleBarLength],round(xres*scaleBarYPos)*[1 1],'-w',...
    'LineWidth',5);
if(showScaleBarText)
    text(round(yres*scaleBarXPos)+1,round(xres*scaleBarYPos)-1, ...
        [num2str(scaleBarInMicrons)  '\mu'],'Color','w',...
        'HorizontalAlignment','Left','VerticalAlignment','Bottom' );
end
hold off;
%axis off;

end
