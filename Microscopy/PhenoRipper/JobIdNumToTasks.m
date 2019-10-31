function taskNums=JobIdNumToTasks(idNum,taskSizes)
%Convenience function to map from a linear index for jobs into indices for multiple parameters being scanned in a grid
numberOfLevels=length(taskSizes);
if(numberOfLevels>1)
    taskNums=zeros(1,numberOfLevels);
    topLevelSize=prod(taskSizes(1:end-1));
    topLevelIdx=ceil(idNum/topLevelSize);
    if(topLevelIdx~=1)
        remIdx=rem(idNum-1,(topLevelIdx-1)*topLevelSize)+1;
    else
        remIdx=idNum;
    end
    taskNums(end)=topLevelIdx;
    taskNums(1:(numberOfLevels-1))=JobIdNumToTasks(remIdx,taskSizes(1:(numberOfLevels-1)));
else
    taskNums=idNum;
end

end
