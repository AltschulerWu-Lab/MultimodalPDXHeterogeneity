function data=rip_image(filenames,global_data,marker_scales,include_bg,...
    foreground_channels, analyze_channels, rescale_param)
% RIP_IMAGE PhenoRip individual images
%   DATA=RIP_IMAGE(FILENAMES,GLOBAL_DATA,MARKER_SCALES,INCLUDE_BG,...
%        FOREGROUND_CHANNELS, ANALYZE_CHANNELS, RESCALE_PARAM)  produces 
%   a phenotypic profile of a specified image based on the fractions of 
%   (previously identified) superblock types it contains.
%   
%
% rip_image arguments: 
%   FILENAMES - cell array of filenames with rows corresponding to different images; the
%   columns correspond to different channels.
%   GLOBAL_DATA - Output of identify_superblock_types that contains information
%   of colors and block types to be used in identification of superblock
%   types
%   NUMBER_OF_SUPERBLOCKS - number of superblock types
%   MARKER_SCALES - An array of size [number of channels x 2]. The first column is the min
%   value of each channel, and the second the max value. Should be the same
%   one used in identify_block_types and identify_superblock_types
%   INCLUDE_BG - a bool which determines if background pixels are used
%   FOREGROUND_CHANNELS - 
%   ANALYZE_CHANNELS - 
%   RESCALE_PARAM - 
%
%  rip_image output: DATA is a structure with fields
%    SUPERBLOCK_PROFILE - a vector containing fractions of different superblock
%    types in the image
%    BLOCK_PROFILE - avector containing the fraction of different block types
%    in the image
%    NUMBER_OF_FOREGROUND_BLOCKS - number of blocks that are considered
%    foreground in the image
%    IMAGE_SUPERBLOCK_STATES - a 2D image showing the superblock type of every
%    superblock in the image. Note: Since superblocks overlap, the central block in 
%    the superblock is used for this annotation.
%    IMAGE_IN_DISCRETE_BLOCK_STATES - a 2D image showing the block type of every
%    block in the image
%    DISTANCE_TO_SUPERBLOCK_CENTROID - for each foreground superblock, the
%    distance to the centroid of the superblock it was assigned to. This is a
%    measure of quality of the superblock assignment.

% ------------------------------------------------------------------------------
% Copyright ??2012, The University of Texas Southwestern Medical Center 
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

%% Preprocessing
if(iscell(filenames))
    isTissue=false;
elseif(isstruct(filenames))
    isTissue=true;
else
    error('Invalid Filename Variable');
end

%Load parameters from global_data (see identify_block_types.m &
%identify_superblock_types.m for description)
block_size=global_data.block_size;
cutoff_intensity=global_data.cutoff_intensity;
number_of_RGB_clusters=global_data.number_of_RGB_clusters;
number_of_block_clusters=global_data.number_of_block_clusters;
number_of_superblocks=size(global_data.superblock_centroids,1);
xres_full=global_data.xres_full;
yres_full=global_data.yres_full;
channels_per_file=global_data.channels_per_file;
A1=global_data.A2;
B1=global_data.B2;
RGB_centroids=global_data.RGB_centroids;

%Load image parameters


if(isTissue)
    number_of_channels=filenames(1).tImg.numberOfChannels;
    number_of_repeats=length(filenames);
  
else
    [number_of_repeats,number_of_channels]=size(filenames);
    number_of_channels=max(number_of_channels,channels_per_file);
end

blocks_nx=floor(xres_full/block_size);
x_offset=floor(rem(xres_full,block_size)/2)+1;
blocks_ny=floor(yres_full/block_size);
y_offset=floor(rem(yres_full,block_size)/2)+1;


%% Reading Images

% Store the block types of all the foreground blocks. Since we don't know
% how many there will be, we-preallocate to the maximum possible number
block_ids_temp=zeros(blocks_nx*blocks_ny*number_of_repeats,1);
%Count the number of foreground blocks in each training image
foreground_blocks_per_image=zeros(number_of_repeats,1);
% foreground_blocks_temp=zeros(blocks_nx*blocks_ny*number_of_repeats,1);

% For each foreground block the neighbor profile is the fractions of the different
% block types in its neighborhood (superblock). This is pre-allocated to
% maximum size
neighbor_profiles_temp=zeros(blocks_nx*blocks_ny*number_of_repeats,number_of_block_clusters+1);

% If rip_image is being called from inside PhenoRipper, it updates the
% timing bar
myhandles=getappdata(0,'myhandles');
if(~isempty(myhandles))
    %progress= get(myhandles.statusbarHandles.ProgressBar, 'Value');
    progress= myhandles.statusbarHandles.ProgressBar.getValue();
end
tStart1=tic;


block_counter=0;
superblock_counter=0;
% Loop over images corresponding to given condition (i.e., images specified
% in filenames). However, PhenoRipper passes one image at a time. so
% number_of_repeats =1 if this function is called from PhenoRipper
for image_counter=1:number_of_repeats 
    %Read and Scale Images 
    
    if(isTissue)
        
        img=read_and_scale_timg(filenames(image_counter),xres_full,yres_full,marker_scales);
        
        
    else
        
        if(channels_per_file>1)
            img=read_and_scale_image(filenames(image_counter),marker_scales,xres_full,yres_full,channels_per_file,xres_full,yres_full,rescale_param(image_counter));
        else
            img=read_and_scale_image(filenames(image_counter,:),marker_scales,xres_full,yres_full,channels_per_file,xres_full,yres_full,rescale_param(image_counter,:));
        end
    end
    number_of_channels=max(number_of_channels,channels_per_file);
    % Crop image to have integer number of blocks in each direction
    cropped_image=img(x_offset:(x_offset+blocks_nx*block_size-1),...
        y_offset:(y_offset+blocks_ny*block_size-1),1:number_of_channels);

    % Calculate intensity of cropped image and identify foreground points
    %intensity=sqrt(sum(double(cropped_image).^2,3)/number_of_channels);
    
    
    %Identify foreground pixels and store their RBG (i.e., multi-channel intensity) values
    %Define the foreground intensity mask based on the channel used to
    %define foreground
    j=0;
    for i=1:length(foreground_channels)
      if(foreground_channels(i))
        j=j+1;
        img2(:,:,j)=cropped_image(:,:,i);
      end
    end
    number_of_channels=j;
    intensity=sqrt(sum(img2.^2,3)/number_of_channels);    
    img2=[];
    %Do the analysis only on the selected channels (reset img to the selected channels)
    j=0;
    for i=1:length(analyze_channels)
      if(analyze_channels(i))
        j=j+1;
        img2(:,:,j)=cropped_image(:,:,i);
      end
    end
    cropped_image=img2;    
    number_of_channels=j;
    
    
    
    bool_foreground_image=(intensity>cutoff_intensity);
    number_of_foreground_points=sum(sum(bool_foreground_image));
    
    %Load the color profiles for the foreground points. A color profile is
    %the vector consisting of levels of pixels in the different channels (each of which lies in
    %[0,100])
    foreground_points=zeros(number_of_foreground_points,number_of_channels);% Color profiles of FG points
    for channel_counter=1:number_of_channels
        temp=squeeze(cropped_image(:,:,channel_counter)); %squeeze extracts a single channel image as a 2D array (instead of 3)
        foreground_points(:,channel_counter)=temp(bool_foreground_image);
    end
    
    
    % Identify the reduced colorstates of all the foreground pixels in
    % image. We already calculated the centroids of the reduced color
    % states, so for each pixel we need to determine which centroid its
    % color profile is closest to
    RGB_distmat=zeros(number_of_foreground_points,number_of_RGB_clusters);
    for reduced_color_state_index=1:number_of_RGB_clusters
        % Instead of looping over every foreground point, we create a
        % matrix with all rows identical, and equal to a centroid vector (in channel space coords).
        % For example if we have 3 channels, and 1000 FG points, this
        % matrix will have dimensions 1000x3, with all rows identical.
        % Subtracting this matrix from the foreground color matrix allows
        % us to calculate distance of all FG points to the centroid.
                rgb_mat=ones(size(foreground_points,1),number_of_channels);
                for i=1:number_of_channels
                    rgb_mat(:,i)=rgb_mat(:,i).*RGB_centroids(reduced_color_state_index,i);
                end
        
        %rgb_mat=repmat(RGB_centroids(reduced_color_state_index,:),number_of_foreground_points,1) ;
        RGB_distmat(:,reduced_color_state_index)=sum((foreground_points-rgb_mat).^2,2);%Calculate distances
       % RGB_distmat(:,reduced_color_state_index)=sum((foreground_points- RGB_centroids(repmat(reduced_color_state_index,number_of_foreground_points,1),:)).^2,2);
    end
    [~,foreground_pixel_states]=min(RGB_distmat,[],2); % Assign each pixel to the closest centroid
    
    % This is an image with the same x,y size as the cropped image, but
    % where each pixel has the value of the reduced color state
    % id. Background points are set to zero.
    image_in_discrete_pixel_states=zeros(block_size*blocks_nx,block_size*blocks_ny);
    image_in_discrete_pixel_states(bool_foreground_image)=foreground_pixel_states;
   
    % Find fraction of foreground points per block to identify foreground
    % blocks
    %avg_block_intensities=A1*intensity*B1/(block_size^2);
    fraction_fg_pixels_in_block=A1*double(bool_foreground_image)*B1/(block_size^2);
    %foreground_blocks=find(avg_block_intensities>cutoff_intensity);
    % Blocks that are more than 50% foreground are considered foreground blocks
    foreground_blocks=find(fraction_fg_pixels_in_block>0.5);
    
    % Describe blocks in terms of the fractions of the reduced color states
    % of constituent pixels
    reduced_color_profiles_of_blocks=zeros(length(foreground_blocks),number_of_RGB_clusters+1);
    % If backround is not used, then the weight of the background (first)
    % component is always zero
    if(include_bg)
        start_cluster=0;
    else
        start_cluster=1;
    end
    
    %  Find fractions of pixels (per block) in the different reduced color
    %  states
    for reduced_color_state_index=start_cluster:number_of_RGB_clusters
        temp=A1*(double(image_in_discrete_pixel_states==reduced_color_state_index))*B1/(block_size^2);
        reduced_color_profiles_of_blocks(:,reduced_color_state_index+1)=temp(foreground_blocks);
    end
    % Normalize so that fractions add up to one. This is needed because if
    % we are not using background points, we do not know the total number
    % of FG points in the last
    rescale_factor=sum(reduced_color_profiles_of_blocks,2);
    for i=1:length(foreground_blocks)
        reduced_color_profiles_of_blocks(i,:)=reduced_color_profiles_of_blocks(i,:)/rescale_factor(i);
    end
    
    % Calculate block types of foreground blocks:
    % Each block is characterized by a vector of
    % length (number_of_reduced_colors+1). The block type centroid vectors also
    % reside in this space. For each block we find the closest block type
    % centroid and assign it to this block type.
    block_distmat=zeros(length(foreground_blocks),number_of_block_clusters);
    for block_type_index=1:number_of_block_clusters
        temp=repmat(global_data.block_centroids(block_type_index,:),length(foreground_blocks),1);
        block_distmat(:,block_type_index) =sum((reduced_color_profiles_of_blocks-temp).^2,2);
    end
    [~,block_ids_in_image]=min(block_distmat,[],2);
    % Store the block_type ids of FG blocks in the image
    block_ids_temp(block_counter+1:block_counter+length(foreground_blocks))=block_ids_in_image;
    
    % This is an image of dim nx x ny, where nx , and ny are the number of
    % blocks in x and y directions respectively 
    % Each point has the value of the block type id
    % Background blocks are set to zero.
    image_in_discrete_block_states=zeros(blocks_nx,blocks_ny);
    image_in_discrete_block_states(foreground_blocks)=block_ids_in_image;
    data.image_in_discrete_block_states=image_in_discrete_block_states;
    
    % Calculate the superblock (neighbor) profiles:
    % For a block, this is the fractions of blocks of different types
    % that are its neighbors
    % In practice we perform this calculation by convoluting the block
    % state image with a filter
    superblock_size=3;
    h=ones(superblock_size);
    neighbor_profiles_in_image=zeros(length(foreground_blocks),number_of_block_clusters+1);
    for block_type_id=start_cluster:number_of_block_clusters
        % is_of_block_type_id is matrix of size nx x ny which is 1
        % only for blocks of type block_type_id
        is_of_block_type_id=double(image_in_discrete_block_states==block_type_id);
        % filtering with h gives the
        % number of block neighbors of chosen block type for each block
        temp=imfilter(is_of_block_type_id,h);
        % Assign fractions for neighbor profiles just for foreground blocks
        neighbor_profiles_in_image(:,block_type_id+1)=temp(foreground_blocks)/(superblock_size^2);
    end
    
    % Normalize the neighbor profiles so they add up to one
    % This is necessary when we are not counting the background blocks
    neighbor_norm=sum(neighbor_profiles_in_image,2);
    for i=1:size(neighbor_profiles_in_image,1)
        neighbor_profiles_in_image(i,:)=neighbor_profiles_in_image(i,:)/neighbor_norm(i);
    end
    
    % Do not use superblocks at the edge, since the number of neighbors is
    % different at the edges
    is_sb_used_image=false(blocks_nx,blocks_ny);%This is a bool array (nx x ny) marking the superblocks to be used
    is_sb_used_image(foreground_blocks)=true;
    is_sb_used_image(1,:)=false;is_sb_used_image(:,1)=false;
    is_sb_used_image(end,:)=false;is_sb_used_image(:,end)=false;
    is_foreground_sb=is_sb_used_image(foreground_blocks);%Reformat is_sb_used_image into a linear vector
    % if(~include_bg)
    %     is_foreground_sb=(is_foreground_sb&neighbor_profiles_in_image(:,1)==0);
    %end
    foreground_superblocks=find(is_foreground_sb);
   
    % Count the number of superblocks in the image
    foreground_blocks_per_image(image_counter)=length(foreground_superblocks);
    % Store the neighbor profiles of just the selected (i.e., not at edge)
    % foreground superblocks
    neighbor_profiles_temp(superblock_counter+1:superblock_counter+length(foreground_superblocks),:)=neighbor_profiles_in_image(foreground_superblocks,:);
 

    %count number of blocks and superblocks
    block_counter=block_counter+length(foreground_blocks);
    superblock_counter=superblock_counter+length(foreground_superblocks);
  
    %Update progress bar if function is called from inside PhenoRipper  
    if(~isempty(myhandles))
        progress=progress+1;
        %tElapsed1=toc(tStart1);
        %tElapsed=myhandles.tElapsed+tElapsed1;
        %files_analyzed=myhandles.files_analyzed+image_counter;
        %time_left=(myhandles.number_of_files-files_analyzed)*tElapsed/(files_analyzed);
        %set(myhandles.statusbarHandles.ProgressBar,...
        %    'Value',progress,'StringPainted','on', 'string',format_time(time_left));
	    %myhandles.statusbarHandles.ProgressBar.setStringPainted(true);
	    %myhandles.statusbarHandles.ProgressBar.setString(format_time(time_left));
  		%myhandles.statusbarHandles.ProgressBar.setValue(progress);
    end
end

% Use only the filled part of block ids and neighbor profiles 
% (since they were pre-allocated to max possible size)
block_ids=block_ids_temp(1:block_counter);
neighbor_profiles=neighbor_profiles_temp(1:superblock_counter,:);

% Calculate superblock types of foreground superblocks:
% Each superblock is characterized by a vector of
% length (number_of_block_types+1). The superblock type centroid vectors also
% reside in this space. For each superblock we find the closest superblock type
% centroid and assign it to this superblock type.
superblock_distmat=zeros(superblock_counter,number_of_superblocks);
for superblock_type_index=1:number_of_superblocks
    temp=repmat(global_data.superblock_centroids(superblock_type_index,:),superblock_counter,1);
    superblock_distmat(:,superblock_type_index)=sum((neighbor_profiles-temp).^2,2);
end
% Each superblock is assigned the superblock_id of the closest superblock type centroid. 
% And its distance from this centroid is stored as a measure of the quality
% of this assignment. This number is used later by PhenoRipper to colorcode
% superblocks by the confidence we have in their superblock type
[distance_to_superblock_centroid,superblock_ids]=min(superblock_distmat,[],2); 
data.distance_to_superblock_centroid=zeros(blocks_nx,blocks_ny);
data.distance_to_superblock_centroid(foreground_superblocks)=distance_to_superblock_centroid;

% This is an image of dim nx x ny, where nx , and ny are the number of
% blocks in x and y directions respectively
% Each block has the value of the superblock type id of the superblock of
% which it the central block.
% Background blocks are set to zero.
image_superblock_states=zeros(blocks_nx,blocks_ny);
image_superblock_states(foreground_blocks(is_foreground_sb))=superblock_ids;
data.image_superblock_states=image_superblock_states;

data.number_of_foreground_blocks=round(superblock_counter/number_of_repeats);

% Calculate the fractions of the different block types for the
% condition/image
data.block_profile=zeros(number_of_block_clusters,1);
for block_type_index=1:number_of_block_clusters
    data.block_profile(block_type_index)=sum(block_ids==block_type_index)/length(block_ids);
end

% Calculate the fractions of the different superblock types for the
% condition/image
data.superblock_profile=zeros(number_of_superblocks,1);
for superblock_type_index=1:number_of_superblocks
    data.superblock_profile(superblock_type_index)=sum(superblock_ids==superblock_type_index)/length(superblock_ids);
end



end

