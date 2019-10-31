# This is an R-script to perform GSVA analysis on the RNAseq data using the hallmark gene set downloaded  from msigdb
# Please change the filenames to point to the locations of the relevant files in your data set
# Load libraries ----
#library(CMSclassifier)
library(XLConnect)
#library("org.Hs.eg.db")
library(plyr)
library(GSVA)
library(GSVAdata)

# Load RNA Data -----
rnaDataDir='/work/bioinformatics/srajar/Leidos/Data/Info_Final/';#CHANGE THIS
rnaDataFile=file.path(rnaDataDir,'normalized_reads.xls');#CHANGE THIS
rnaData<-readWorksheetFromFile(rnaDataFile,sheet=1);
rnaExpression<-as.matrix(rnaData[,3:39]);
geneSymbols<-as.matrix(rnaData[,1]);
rownames(rnaExpression)<-geneSymbols;

# Load GSEA gene sets ----
geneSetDir=rnaDataDir;
hallmarkGeneSetFile=file.path(geneSetDir,'h.all.v6.0.symbols.gmt');
hallmarkGmt<-getGmt(hallmarkGeneSetFile,geneIdType = SymbolIdentifier());
# Run GSVA ----
gsvaHallmarkRes<-gsva(rnaExpression,hallmarkGmt,min.sz=1,max.sz=500,mx.diff=FALSE,kcdf='Poisson');
write.table(gsvaHallmarkRes,'gsvaHallmark_rnaseq.txt',quote = FALSE,sep='\t',col.names = NA,row.names=TRUE)
