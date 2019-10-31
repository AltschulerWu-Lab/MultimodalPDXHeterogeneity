function diffVals=Profile_Difference(refProb,sampleProbs,method)

%Convenience function to allow calculating the difference in PhenoRipper profiles
% between the true refProb and estimated profiles stored in sampleProbs based on
% the difference measure specified in method
	
allowedMethods={'KL','MaxFrac','AvgFrac','L1','AvgDiff','MaxDiff'};
validatestring(method,allowedMethods);

if(size(refProb,2)~=size(sampleProbs,2))
    error('Size Mismatch');
end



switch(method)
    case 'KL'
        diffVals=-sum(bsxfun(@times,bsxfun(@minus, log(sampleProbs),log(refProb)),refProb),2);
    case 'MaxFrac'
        diffVals=max(abs(bsxfun(@rdivide, sampleProbs,refProb)-1),[],2);
    case 'AvgFrac'
        diffVals=mean(abs(bsxfun(@rdivide, sampleProbs,refProb)-1),2);
    case 'L1'
        diffVals=sum(abs(bsxfun(@minus, sampleProbs,refProb)),2);
    case 'AvgDiff'
        diffVals=mean(abs(bsxfun(@minus, sampleProbs,refProb)),2);
    case 'MaxDiff'
        diffVals=max(abs(bsxfun(@minus, sampleProbs,refProb)),[],2);
    otherwise
        error('Invalid Method');
end
end
