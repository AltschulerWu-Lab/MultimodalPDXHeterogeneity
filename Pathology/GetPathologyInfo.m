function [pctTumor,pctStroma,pctNecrosis,heScores]=GetPathologyInfo()

%A convenience function to return the percentage of tumor/stroma/necrosis for the 36 samples based on pathologist characterization.	
params=GetParams('pathology');

heScores=readtable(params.pathology.scottHEScoringFile);
%heScores=heScores(1:36,:);

pctNecrosis=heScores.Pct_Necrosis;
pctTumor=heScores.Pct_Tumor;
pctStroma=heScores.Pct_Stroma;

end
