function data=identify_superblock_types(filenames,global_data,number_of_superblocks,marker_scales,include_bg,...
    foreground_channels, analyze_channels, rescale_param)
% IDENTIFY_SUPERBLOCK_TYPES Identify PhenoRipper superblock types 
%   DATA=IDENTIFY_SUPERBLOCK_TYPES(FILENAMES,GLOBAL_DATA,...
%        NUMBER_OF_SUPERBLOCKS,MARKER_SCALES,INCLUDE_BG,...
%        FOREGROUND_CHANNELS, ANALYZE_CHANNELS, RESCALE_PARAM) 
%   Identify PhenoRipper superblock types on a sample set of images using
%   previously identified block/color types
%
% identify_superblock_type arguments: 
%   FILENAMES - cell array of filenames with rows corresponding to different images; the
%   columns correspond to different channels.
%   GLOBAL_DATA - Output of identify_block_types that contains information
%   of colors and block types to be used in identification of superblock
%   types
%   NUMBER_OF_SUPERBLOCKS - number of superblock types
%   MARKER_SCALES - An array of size [number of channels x 2]. The first column is the min
%   value of each channel, and the second the max value. Should be the same
%   one used in identify_block_types
%   INCLUDE_BG - a bool which determines if background pixels are used
%   FOREGROUND_CHANNELS - 
%   ANALYZE_CHANNELS - 
%   RESCALE_PARAM - 
%
%   identify_superblock_types output: DATA is a structure with 
%   fields output by identify_block_types and the additional fields
%   SUPERBLOCK_CENTROIDS - array of size [number_of_SUPERBLOCK_clusters x (number_of_BLOCK_clusters+1)]. 
%   Each row describes the fractions of the the different block types for a
%   distinct superblock type. The 1st column contains the fraction of background blocks
%   SUPERBLOCK_REPRESENTATIVES - a cell array containing representative images
%   of the different superblock types
%   IMAGE_IN_DISCRETE_BLOCK_STATES - a 2D image showing the block type of every
%   block in the last image analyzed
%

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
%%




%% Preprocessing
if(iscell(filenames))
    isTissue=false;
elseif(isstruct(filenames))
    isTissue=true;
else
    error('Invalid Filename Variable');
end

%Load parameters from global_data
block_size=global_data.block_size;
cutoff_intensity=global_data.cutoff_intensity;
number_of_RGB_clusters=global_data.number_of_RGB_clusters;
number_of_block_clusters=global_data.number_of_block_clusters;
xres_full=global_data.xres_full;xres_crop=global_data.xres_crop;
yres_full=global_data.yres_full;yres_crop=global_data.yres_crop;
A1=global_data.A1;
B1=global_data.B1;

number_of_superblock_representatives=3; % Number of sample images of superblock types stored for display later 

%Load image parameters
if(isTissue)
    number_of_channels=filenames(1).tImg.numberOfChannels;
    number_of_training_images=length(filenames);
  
else
    [number_of_training_images,number_of_channels]=size(filenames);
       
end


channels_per_file=global_data.channels_per_file;

blocks_nx=floor(xres_crop/block_size);
x_offset=floor(rem(xres_crop,block_size)/2)+1;%because the image may not contain an integer number of blocks
blocks_ny=floor(yres_crop/block_size);
y_offset=floor(rem(yres_crop,block_size)/2)+1;%because the image may not contain an integer number of blocks


max_training_superblocks_per_image=1000;
number_of_channels=max(number_of_channels,channels_per_file);
number_of_analysis_channels=nnz(analyze_channels);
number_of_fg_channels=nnz(foreground_channels);

%% Reading Images

% Store the block types of all the foreground blocks. Since we don't know
% how many there will be, we-preallocate to the maximum possible number
%block_ids_temp=zeros(blocks_nx*blocks_ny*number_of_training_images,1);
%Count the number of foreground blocks in each training image
foreground_blocks_per_image=zeros(number_of_training_images,1);
%foreground_blocks_temp=zeros(blocks_nx*blocks_ny*number_of_training_images,1);

% Variables used to extract superblock representatives. Instead of storing the 
% superblock images for every foreground superblock (which will be memory
% itensive), we just store image index and its (x,y) position in the image. 
% Then if a block is chosen as a superblock representative, we can use these variables to extract its image 
% As before we pre-allocate these to maximum number
image_number_of_block_temp=zeros(max_training_superblocks_per_image*number_of_training_images,1);
position_of_block_temp=zeros(max_training_superblocks_per_image*number_of_training_images,2);

% For each foreground block the neighbor profile is the fractions of the different
% block types in its neighborhood (superblock). This is pre-allocated to
% maximum size
neighbor_profiles_temp=zeros(max_training_superblocks_per_image*number_of_training_images,number_of_block_clusters+1);

block_counter=0;
superblock_counter=0;

for image_counter=1:number_of_training_images
    %Read and Scale Images    
    if(isTissue)
        
        img=read_and_scale_timg(filenames(image_counter),xres_crop,yres_crop,marker_scales);
        
        
    else
        if(channels_per_file>1)
            img=read_and_scale_image(filenames(image_counter),marker_scales,xres_full,yres_full,channels_per_file,xres_crop,yres_crop,rescale_param(image_counter));
        else
            img=read_and_scale_image(filenames(image_counter,:),marker_scales,xres_full,yres_full,channels_per_file,xres_crop,yres_crop,rescale_param(image_counter,:));
        end
    end
    
    % Crop image to have integer number of blocks in each direction
    cropped_image=img(x_offset:(x_offset+blocks_nx*block_size-1),...
        y_offset:(y_offset+blocks_ny*block_size-1),1:number_of_channels);
    
    % Calculate intensity of cropped image and identify foreground points
    %intensity=sqrt(sum(double(cropped_image).^2,3)/number_of_analysis_channels);  
    
    %Identify foreground pixels and store their RBG (i.e., multi-channel intensity) values
    %Define the foreground intensity mask based on the channel used to
    %define foreground
    intensity=sqrt(sum(cropped_image(:,:,foreground_channels~=0).^2,3)/number_of_fg_channels);
     
    %Do the analysis only on the selected channels (reset img to the selected channels)
    cropped_image=cropped_image(:,:,analyze_channels~=0);
    
    is.foreground=(intensity>cutoff_intensity);
    number_of_foreground_points=sum(sum(is.foreground));
    
    %Load the color profiles for the foreground points. A color profile is
    %the vector consisting of levels of pixels in the different channels (each of which lies in
    %[0,100])
    foreground_points=zeros(number_of_foreground_points,number_of_analysis_channels); % Color profiles of FG points
    for channel_counter=1:number_of_analysis_channels
        temp=squeeze(cropped_image(:,:,channel_counter)); %squeeze extracts a single channel image as a 2D array (instead of 3)
        foreground_points(:,channel_counter)=temp(intensity>cutoff_intensity);
    end
    
    % Identify the reduced colorstates of all the foreground pixels in
    % image. We already calculated the centroids of the reduced color
    % states, so for each pixel we need to determine which centroid its
    % color profile is closest to
    RGB_distmat=zeros(number_of_foreground_points,number_of_RGB_clusters);% Distance of FG pixels to color centroids
    for reduced_color_state_index=1:number_of_RGB_clusters
%         rgb_mat=ones(size(foreground_points,1),number_of_analysis_channels);
%         for i=1:number_of_analysis_channels
%             rgb_mat(:,i)=rgb_mat(:,i).*global_data.RGB_centroids(reduced_color_state_index,i);
%         end
        % Instead of looping over every foreground point, we create a
        % matrix with all rows identical, and equal to a centroid vector (in channel space coords).
        % For example if we have 3 channels, and 1000 FG points, this
        % matrix will have dimensions 1000x3, with all rows identical.
        % Subtracting this matrix from the foreground color matrix allows
        % us to calculate distance of all FG points to the centroid.
        rgb_mat=repmat(global_data.RGB_centroids(reduced_color_state_index,:),number_of_foreground_points,1) ;
        RGB_distmat(:,reduced_color_state_index)=sum((foreground_points-rgb_mat).^2,2);%Calculate distances
    end
    [~,foreground_pixel_states]=min(RGB_distmat,[],2);% Assign each pixel to the closest centroid
        
    % This is an image with the same x,y size as the cropped image, but
    % where each pixel has the value of the reduced color state
    % id. Background points are set to zero.
    image_in_discrete_pixel_states=zeros(block_size*blocks_nx,block_size*blocks_ny);
    image_in_discrete_pixel_states(is.foreground)=foreground_pixel_states;
    
    % Find fraction of foreground points per block to identify foreground
    % blocks
    %avg_block_intensities=A1*intensity*B1/(block_size^2);
    %foreground_blocks=find(avg_block_intensities>cutoff_intensity);
    fraction_fg_pixels_in_block=A1*double(intensity>cutoff_intensity)*B1/(block_size^2);
    % Blocks that are more than 50% foreground are considered foreground blocks
    foreground_blocks=find(fraction_fg_pixels_in_block>0.5);
    
    % Describe blocks in terms of the fractions of the reduced color states
    % of constituent pixels
    reduced_color_profiles_of_blocks=zeros(length(foreground_blocks),number_of_RGB_clusters+1);%+1 because of background state
    % If backround is not used, then the weight of the background (first) component is
    % always zero
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
    for block_cluster=1:number_of_block_clusters
        temp=repmat(global_data.block_centroids(block_cluster,:),length(foreground_blocks),1);
        block_distmat(:,block_cluster) =sum((reduced_color_profiles_of_blocks-temp).^2,2);
    end
    [~,block_ids_in_image]=min(block_distmat,[],2);
    % Store the block_type ids of FG blocks in the image
    %block_ids_temp(block_counter+1:block_counter+length(foreground_blocks))=block_ids_in_image;
    
    
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
    for block_type_id=start_cluster:number_of_block_clusters %Note: if we discard background, fraction of BG blocks is zero
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
%     if(~include_bg)
%         is_foreground_sb=(is_foreground_sb&neighbor_profiles_in_image(:,1)==0);
%     end


    
    foreground_superblocks=find(is_foreground_sb);
    selected_fg_sb=randsample(length(foreground_superblocks),min(max_training_superblocks_per_image,length(foreground_superblocks)));
    foreground_superblocks=foreground_superblocks(selected_fg_sb);
    
    
    is_selected_fg_sb=false(length(foreground_blocks),1);
    is_selected_fg_sb(foreground_superblocks)=true;
    is_foreground_sb=is_foreground_sb& is_selected_fg_sb;
    
    % Count the number of superblocks in the image
    foreground_blocks_per_image(image_counter)=length(foreground_superblocks);
    % Store the neighbor profiles of just the selected (i.e., not at edge)
    % foreground superblocks
    neighbor_profiles_temp(superblock_counter+1:superblock_counter+length(foreground_superblocks),:)=neighbor_profiles_in_image(foreground_superblocks,:);
    
    
    % Variables used to extract superblock representatives. Instead of storing the 
    % superblock images for every foreground superblock (which will be memory
    % itensive), we just store image index and its (x,y) position in the image. 
    % Then if a block is chosen as a superblock representative, we can use these variables to extract its image 
    % As before we pre-allocate these to maximum number
    image_number_of_block_temp(superblock_counter+1:superblock_counter+length(foreground_superblocks))=...
        image_counter;
    [x_pos,y_pos]=ind2sub([blocks_nx,blocks_ny],foreground_blocks);%x and y positions of block in block coordinates
    %find(fraction_fg_pixels_in_block>0.5);%repeated operation,can be speeded up
    %[x_pos,y_pos]=find(avg_block_intensities>cutoff_intensity);%repeated operation,can be speeded up
    x_pos=x_pos(is_foreground_sb); %drop edge superblocks
    y_pos=y_pos(is_foreground_sb); %drop edge superblocks
    % Save the block positions in pixel coords (conversion from block
    % coords)
    try
    position_of_block_temp(superblock_counter+1:superblock_counter+length(foreground_superblocks),:)=...
        [(x_pos'-1)*block_size+x_offset;(y_pos'-1)*block_size+y_offset]'; 
    catch err
      disp(err.message);
      %error('problem');
    end
    
    %count number of blocks and superblocks
    block_counter=block_counter+length(foreground_blocks);
    superblock_counter=superblock_counter+length(foreground_superblocks);
    
end

% Use only the filled part of block ids, image number of block, position 
% of block in image, and neighbor profiles (since they were pre-allocated 
% to max possible size)
%block_ids=block_ids_temp(1:block_counter);
image_number_of_block=image_number_of_block_temp(1:superblock_counter);
position_of_block=position_of_block_temp(1:superblock_counter,:);
neighbor_profiles=neighbor_profiles_temp(1:superblock_counter,:);

%mean_superblock_profile=mean(neighbor_profiles);
%data.mean_superblock_profile=mean_superblock_profile;

% Use k-means to cluster the neighbor profiles and identify the superblock
% types
[superblock_ids,data.superblock_centroids,~,superblock_distances]=kmeans(neighbor_profiles,number_of_superblocks...
    ,'emptyaction','singleton','start','cluster');

%Added hack start%%%%%%%%%%%
% if(size(neighbor_profiles,1)<=10000)
%   [superblock_ids,data.superblock_centroids,~,superblock_distances]=kmeans(neighbor_profiles,number_of_superblocks...
%      ,'emptyaction','singleton','start','cluster');
% else
%     chosen_superblocks=randsample(size(neighbor_profiles,1),10000);
%     neighbor_profiles_temp=neighbor_profiles(chosen_superblocks,:);
%     [~,data.superblock_centroids,~,~]=kmeans(neighbor_profiles_temp,number_of_superblocks...
%     ,'emptyaction','singleton','start','cluster');
%     superblock_distmat=zeros(size(neighbor_profiles,1),number_of_superblocks);
%     for superblock_type_index=1:number_of_superblocks
%         temp=repmat(data.superblock_centroids(superblock_type_index,:),size(neighbor_profiles,1),1);
%         superblock_distmat(:,block_type_index) =sum((neighbor_profiles-temp).^2,2);
%     end
%     [~,superblock_ids]=min(block_distmat,[],2);
% end
%Added hack END %%%%%%%%%%%%
  
% Return the fractions of block types in training images
% data.block_profile=zeros(number_of_block_clusters,1);
% for block_cluster=1:number_of_block_clusters
%    data.block_profile(block_cluster)=sum(block_ids==block_cluster)/length(block_ids);
% end

% Find the fractions of superblock types in each image, this information is
% used in selecting superblock type representatives
superblock_profiles=zeros(number_of_training_images,number_of_superblocks);
for image_number=1:number_of_training_images
    ids_in_image=superblock_ids(image_number_of_block==image_number);
    for sb_num=1:number_of_superblocks
        superblock_profiles(image_number,sb_num)=nnz(ids_in_image==sb_num);
    end
    superblock_profiles(image_number,:)=superblock_profiles(image_number,:)/...
        sum(superblock_profiles(image_number,:));
end
superblock_profiles(isnan(superblock_profiles))=eps;
%Picking representatives
% 1) For each block store filenumber and (x,y) location (remember cropping)
% 

% Optimal File Opening Steps:
% 1) Identify blocks within acceptable range of top and determine image and location
% 2) Create a FileVsBlockType matrix containing number of acceptable representatives for each superblock per file
% 3) If i is the indicator function to include a file, then minimize sum(i) such that i.*f(:,j)>number_of_reps for each j
% 4) Since all files will need to be opened, just go in order and load up all possible reps per file
%FvB=zeros(number_of_training_images,number_of_superblocks);
distance_cutoffs=zeros(number_of_superblocks,1);
locations=cell(number_of_training_images,number_of_superblocks);

% for i=1:number_of_superblocks
%   distance_cutoffs(i)=max(mink(superblock_distances(:,i),number_of_superblock_representatives));
%    
%    acceptable_blocks=find(superblock_distances(:,i)<=distance_cutoffs(i));
%    for j=1:length(acceptable_blocks)
%        block_num=acceptable_blocks(j);
%        FvB(image_number_of_block(block_num),i)=FvB(image_number_of_block(block_num),i)+1;
%        locations{image_number_of_block(block_num),i}=[locations{image_number_of_block(block_num),i};...
%            position_of_block(block_num,:)];
%    end
% end

%change units from pixels to blocks
position_of_block_in_block_coords=zeros(size(position_of_block)); 
position_of_block_in_block_coords(:,1)=(position_of_block(:,1)-x_offset)/block_size+1;
position_of_block_in_block_coords(:,2)=(position_of_block(:,2)-y_offset)/block_size+1;
%Identify superblocks not at the edge assuming we include and additional border
%of 2*block_size around the superblock
is_not_edge=(position_of_block_in_block_coords(:,1)>3)&(position_of_block_in_block_coords(:,2)>3)...
        &(position_of_block_in_block_coords(:,1)<blocks_nx-3)&(position_of_block_in_block_coords(:,2)<blocks_ny-3);
for i=1:number_of_superblocks
    % blocks_of_type=find(superblock_ids==i);
    ten_percentile_dist=prctile(superblock_distances((superblock_ids==i)&is_not_edge,i),1);
    nth_closest_distance=max(mink(superblock_distances(is_not_edge,i),number_of_superblock_representatives));
    distance_cutoffs(i)=max(ten_percentile_dist,nth_closest_distance);
    
    acceptable_blocks=find((superblock_distances(:,i)<=distance_cutoffs(i))&(superblock_ids==i)&is_not_edge);
    
    [~,image_order]=sort(superblock_profiles(:,i),'descend');
    rep_count=0;
    image_count=1;
    while((rep_count<number_of_superblock_representatives)&&(image_count<=number_of_training_images))
        image_number=image_order(image_count);
        selected_blocks=acceptable_blocks(image_number_of_block(acceptable_blocks)==image_number);
        selected_block_distances=superblock_distances(selected_blocks,i);
        [~,block_order]=sort(selected_block_distances,'ascend');
        needed_length=min(number_of_superblock_representatives-rep_count,length(selected_blocks));
        %selected_blocks=selected_blocks(1:needed_length);
        for j=1:needed_length
            block_num=selected_blocks(block_order(j));
            locations{image_number,i}=[locations{image_number,i};...
                position_of_block(block_num,:)];
        end
        rep_count=rep_count+needed_length;
        image_count=image_count+1;
    end
end


% [number_of_training_images,number_of_analysis_channels]=size(filenames);
% number_of_analysis_channels=max(number_of_analysis_channels,channels_per_file);


%included_files=find(bintprog(ones(number_of_training_images,1),-FvB',-number_of_superblock_representatives*ones(number_of_superblocks,1))>0.5);
included_files=1:number_of_training_images;
superblock_representatives=cell(number_of_superblocks,number_of_superblock_representatives);

supr_counter=zeros(number_of_superblocks,1);
for file_num=1:length(included_files)
    if(any(~cellfun('isempty',locations(file_num,:))))
        %img=zeros(xres_crop,yres_crop,number_of_analysis_channels);
        
        %img_temp=zeros(xres_crop,yres_crop,number_of_analysis_channels);
        %if((xres_full~=xres_crop)||(yres_full~=yres_crop))
        %    x1=ceil((xres_full-xres_crop)/2);
        %    y1=ceil((yres_full-yres_crop)/2);
        %    x2=x1+xres_crop-1;
        %    y2=y1+yres_crop-1;
        %else
         %   x1=1;x2=xres_full;
         %   y1=1;y2=yres_full;
        %end
        
        if(isTissue)
            
            img=read_and_scale_timg(filenames(included_files(file_num)),xres_crop,yres_crop,marker_scales);
            
            
        else
            
            
            if(channels_per_file>1)
                %img=double(imread2(cell2mat(filenames(image_counter))));
                %temp=imread2(cell2mat(filenames(included_files(file_num),1)));
                
                img=read_and_scale_image(filenames(included_files(file_num)),marker_scales,xres_full,yres_full,channels_per_file,xres_crop,yres_crop,rescale_param(included_files(file_num)));
                
                %img=temp(x1:x2,y1:y2,:);
            else
                %            for channel=1:number_of_analysis_channels
                %img(:,:,channel_counter)=imread2(cell2mat(filenames(image_counter,channel_counter)));
                %USE IMREAD FOR SINGLE CHANNEL ALWAYS
                
                img=read_and_scale_image(filenames(included_files(file_num),:),marker_scales,xres_full,yres_full,channels_per_file,xres_crop,yres_crop,rescale_param(included_files(file_num),:));
                %temp=imread(cell2mat(filenames(included_files(file_num),channel)));
                %img(:,:,channel)=temp(x1:x2,y1:y2);
                %           end
            end
        end
        %     for channel=1:number_of_analysis_channels
        %        img(:,:,channel)=imread2(cell2mat(filenames(included_files(file_num),channel)));
        %
        %     end
        for i=1:number_of_superblocks
            number_of_matches=size(locations{included_files(file_num),i},1);
            if(number_of_matches>0)
                location_info=locations{included_files(file_num),i};
                for j=1:min(number_of_matches,number_of_superblock_representatives-supr_counter(i))
                    %rep_img=zeros(3*block_size,3*block_size,number_of_analysis_channels);
                    x1=max(location_info(j,1)-3*block_size,1);
                    x2=min(location_info(j,1)+4*block_size-1,xres_crop);
                    y1=max(location_info(j,2)-3*block_size,1);
                    y2=min(location_info(j,2)+4*block_size-1,yres_crop);
                    supr_counter(i)=supr_counter(i)+1;
                    superblock_representatives{i,supr_counter(i)}=img(x1:x2,y1:y2,:);
                    
                end
            end
           
        end
    end
end



data.superblock_representatives=superblock_representatives;

end

