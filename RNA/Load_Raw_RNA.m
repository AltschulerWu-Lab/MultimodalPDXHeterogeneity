function [combinedVals,singleChipVals,geneInfo,combinedData]=Load_Raw_RNA()

% Convenience function to combine RNA read count data from individual chips 
% OUTPUT
%    combinedVals   - Numeric array with rows corresponding to genes and columns to samples containing the combined number of reads from all (4) replicate chips
%    singleChipVals - numeric aray of size nuber_of_genes_x_number_of_chipsxnumber_of_samples with reads on each chip
%   geneInfo   - A table providing information on the different genes (rows correspond to genes in combinedVals).
%   combinedData      - a cell array with elements containing tables filled with the raw data for each sample.
params=GetParams({'RNA','samples'});
dataDir=params.rna.rootDir;
numberOfChips=params.rna.numberOfChips;
numberOfGenes=params.rna.numberOfGenes;
numberOfSamples=params.samples.numberOfSamples;
combinedData=cell(numberOfSamples,1);

singleChipVals=zeros(numberOfGenes, numberOfChips,numberOfSamples);
combinedVals=zeros(numberOfGenes,numberOfSamples);


%%
for sampleNumber=1:numberOfSamples
    %tic;
    sampleRnaDir=fullfile(dataDir,['RNA_Main_Sample' num2str(sampleNumber)]);
    combinedChipFile=fullfile(sampleRnaDir,['IonXpress_' sprintf('%0.3d',sampleNumber) ...
        '_CombineAlignments_CA_4_chips_combined_06-01-17_UCSF-LCM_001.amplicon.cov.txt']);
    
    combinedData{sampleNumber}=sortrows(readtable(combinedChipFile),'attributes');
    combinedVals(:,sampleNumber)=combinedData{sampleNumber}.overlaps;
    
    chipNumbers=arrayfun(@num2str,361:364,'Unif',false);
    singleChipData=cell(numberOfChips,1);
    
    for chipCounter=1:numberOfChips
        singleChipFile=dir(fullfile(sampleRnaDir,['*' chipNumbers{chipCounter} '.amplicon.cov.txt']));
        if(numel(singleChipFile)~=1)
            error('File name error');
        else
            singleChipFile=fullfile(sampleRnaDir,singleChipFile(1).name);
        end
        singleChipData{chipCounter}=readtable(singleChipFile);
        [areCommon,posInCombined]=ismember( combinedData{sampleNumber}.contig_id,...
            singleChipData{chipCounter}.contig_id);
        if(all(areCommon))
            singleChipVals(:,chipCounter,sampleNumber)=...
                singleChipData{chipCounter}.overlaps(posInCombined);
        else
            error('Missing Contigs');
        end
    end
    
    
    if(~isequal(combinedData{sampleNumber}.contig_id,combinedData{1}.contig_id))
        error('Mismatch in Gene List across samples');
    end
    
    %toc;
    %disp(['Sample ' num2str(sampleNumber) ' done!']);
end

geneName=cell2table(cellfun(@(x) x{2}, ...
    regexp(combinedData{1}.attributes,'=','split'),'Unif',false),...
    'VariableNames',{'geneName'});
geneInfo=[geneName,combinedData{1}(:,[1,2,3,4,6])];

end
