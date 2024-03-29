function image = imread2(filename,is_multi_channel,varargin)
%IMREAD2 Read image from graphics file.
%   image = IMREAD(FILENAME) reads a grayscale or color image from the file
%   specified by the string FILENAME. If the file is not in the current
%   directory, or in a directory on the MATLAB path, specify the full
%   pathname.
%
%   This function read image using the Matlab imread function except for
%   the grayscale TIFF files where tiffread2 from Francois Nedelec is used.
%
%   The return value image is an array containing the image data. If the file
%   contains a grayscale image, image is an M-by-N array. If the file contains
%   a truecolor image, image is an M-by-N-by-3 array.
%   For TIFF files containing color images that use the CMYK color space,
%   image is an M-by-N-by-4 array.
%   For TIFF files containing stack (like stk), image is an M-by-N-by-P array
%   where P is the number of Stack (hoppefully Channels)
%
%   See TIFF in the Format-Specific Information section for more
%   information.
%
%
% ------------------------------------------------------------------------------
% Copyright ©2012, The University of Texas Southwestern Medical Center
% Authors:
% Satwik Rajaram and Benjamin Pavie for the Altschuler and Wu Lab
% For latest updates, check: < http://www.PhenoRipper.org >.
%
% All rights reserved.
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details:
% < http://www.gnu.org/licenses/ >.
%
% ------------------------------------------------------------------------------
%%
rescale_param=[];
if(nargin==3)
    rescale_param=varargin{1};
else
    rescale_param=[];
end
if(~is_multi_channel)
    image=imread(filename);
    if(~isempty(rescale_param))
        A=rescale_param{1};
        B=rescale_param{2};
        image=A*image+B;
        image(image<0)=0;
    end
    return;
end
imageInfo=imfinfo(filename);
if(~strcmpi(imageInfo.Format,'tif')||strcmpi(imageInfo.ColorType,'truecolor'))
    image=imread(filename);
else
    try
        tmp=tiffread2(filename);
        if (~isImagesSameSize(tmp))
            %errordlg('This format is NOT supported');
            error('Format Not Supported! All Images in the TIFF stack MUST have the same size');
        end
        image=zeros(tmp(1).height,tmp(1).width,length(tmp));
        for i=1:length(tmp)
            %image(:,:,i)=tmp(i).data;
            image(:,:,i)=tmp(i).data;
        end
    catch %Tiffread2 do not open compressed Tiff files
        % so try a last time to use matlab library that work only for
        % single tiff channel.
        image=imread(filename);
    end
end
if(~isempty(rescale_param))
    for i=1:size(image,3)
        %If multi channel, use the same parameters RIGHT NOW
        if(is_multi_channel)
            try
                A=rescale_param{1}{1};
                B=rescale_param{1}{2};
            catch
               A=1;B=0; 
            end
                else
                    A=rescale_param{i}{1};
                    B=rescale_param{i}{2};
            end
            tmp=image(:,:,i);
            tmp=A*tmp+B;
            tmp(tmp<0)=0;
            image(:,:,i)=tmp;
        end
    end
    
    
    function isOK=isImagesSameSize(tmp)
    isOK=true;
    if(length(tmp)>1)
        height=tmp(1).height;
        width=tmp(1).width;
        for j=2:length(tmp)
            if(height==tmp(j).height&&width==tmp(j).width)
                continue;
            else
                isOK=false;
                return
            end
        end
    end
    
