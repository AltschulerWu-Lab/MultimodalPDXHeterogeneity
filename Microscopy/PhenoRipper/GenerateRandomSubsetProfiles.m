function subsetProfiles=GenerateRandomSubsetProfiles(profilesSuperSet,nFGVals,...
    subsetSize,numberOfRandomizations,profileGroups)
% Function to generate random sub-samplings (based on a grouping variable) and calculate their profiles

% INPUT
% profilesSuperSet - a 2D array consisting of the full collection of profiles from 
% 		     which we are subsampling. Rows correspond to profiles from different 
% 		     areas while columns represent different components
% nFGVals          - A vector with number of elements equal to rows of profilesSuperSet, 
% 		     containing the number of foreground blocks in the corresponding area
% subsetSize       - the number of areas to be sampled to generate a profile
% numberOfRandomizations - the number of times the sampling should be performed
% profileGroups    - A vector of same size as nFGVals indicating the grouping variable for each area.
% 		     Each subsampling will select areas which share the same grouping variable.	
%
% OUTPUT
% subsetProfiles   - a 2D matrix of size numberOfRandomizationsxprofileLength, with each row correponding to the
% 			profile of a random sampling.	
[numberOfProfiles,profileLength]=size(profilesSuperSet);

if(nargin>4)
    [gNum,gName]=grp2idx(profileGroups);
    numberOfGroups=length(gName);
    gIdx=cell(numberOfGroups,1);
    for g=1:numberOfGroups
        gIdx{g}=find(gNum==g);
    end
else
    numberOfGroups=1;
end

profilesScaled=bsxfun(@times,profilesSuperSet,nFGVals);

subsetProfiles=zeros(numberOfRandomizations,profileLength);

for randCounter=1:numberOfRandomizations
    if(numberOfGroups>1)
        g=randi(numberOfGroups);
        if(subsetSize<numel(gIdx{g}))
            %profilesSelected=randsample(gIdx{g},subsetSize);
            profilesSelected=datasample(gIdx{g},subsetSize,'Replace',false);
        else
            profilesSelected=[];
        end
    else
        profilesSelected=randperm(numberOfProfiles,subsetSize);
    end
    if(~isempty(profilesSelected))
        %profileWeights=nFGVals(profilesSelected);
        %profileWeights=profileWeights(:)/sum(profileWeights);
        
        %temp=bsxfun(@times,profilesSuperSet(profilesSelected,:),profileWeights);
        %subsetProfiles(randCounter,:)=sum(temp,1);
        
        subsetProfiles(randCounter,:)=sum(profilesScaled(profilesSelected,:)/...
            sum(nFGVals(profilesSelected)),1);
        %subsetProfiles(randCounter,:)=sum(bsxfun(@times,profilesSuperSet(profilesSelected,:),profileWeights),1);
    else
         subsetProfiles(randCounter,:)=NaN;
    end
end




end
