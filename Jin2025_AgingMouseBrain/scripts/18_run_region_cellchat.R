#!/usr/bin/env Rscript
.libPaths(c("C:/Users/57265/AppData/Local/R/win-library/4.4", .libPaths()))
suppressPackageStartupMessages({library(CellChat);library(Matrix);library(data.table)})
args<-commandArgs(trailingOnly=TRUE); if(length(args)<2) stop("region age")
region<-args[1]; age<-args[2]; root<-"D:/vscode/Jin2025_AgingMouseBrain"
inp<-file.path(root,"data/cellchat_region",region,age)
out<-file.path(root,"results/cellchat_region",region,age);dir.create(out,recursive=TRUE,showWarnings=FALSE)
sink(file.path(out,"run.log"),split=TRUE);on.exit(sink(),add=TRUE);set.seed(20260705)
counts<-as(readMM(gzfile(file.path(inp,"counts.mtx.gz"))),"dgCMatrix")
genes<-fread(file.path(inp,"genes.csv"))[[1]];meta<-fread(file.path(inp,"metadata.csv"),data.table=FALSE)
rownames(counts)<-genes;colnames(counts)<-meta$cell_label;rownames(meta)<-meta$cell_label
lib<-Matrix::colSums(counts);norm<-counts%*%Diagonal(x=1e4/lib);norm@x<-log1p(norm@x)
cc<-createCellChat(norm,meta=meta,group.by="cell_type_major");cc@DB<-CellChatDB.mouse
cc<-subsetData(cc);cc<-identifyOverExpressedGenes(cc);cc<-identifyOverExpressedInteractions(cc)
cc<-computeCommunProb(cc,type="triMean",raw.use=TRUE,population.size=FALSE,nboot=100,seed.use=20260705)
cc<-filterCommunication(cc,min.cells=10);cc<-computeCommunProbPathway(cc);cc<-aggregateNet(cc)
cc<-tryCatch(netAnalysis_computeCentrality(cc,slot.name="netP"),error=function(e)cc)
saveRDS(cc,file.path(out,paste0("cellchat_",region,"_",age,".rds")),compress=FALSE)
lr<-subsetCommunication(cc,slot.name="net");pw<-subsetCommunication(cc,slot.name="netP")
fwrite(lr,file.path(out,"lr_interactions.csv"));fwrite(pw,file.path(out,"pathway_interactions.csv"))
fwrite(as.data.frame(as.table(cc@net$count)),file.path(out,"network_count.csv"))
fwrite(as.data.frame(as.table(cc@net$weight)),file.path(out,"network_weight.csv"))
focus<-lr[grepl("TGFB|CX3|CSF|IL34|COMPLEMENT|MIF|SPP1|GRN|PSAP|GAS6|PROS1|AXL|MERTK|CCL|CXCL|APOE|TREM",paste(lr$pathway_name,lr$interaction_name_2,lr$ligand,lr$receptor),ignore.case=TRUE),]
fwrite(focus,file.path(out,"microglia_focus_interactions.csv"))
fwrite(data.frame(region=region,age=age,cells=ncol(norm),lr=nrow(lr),pathways=length(cc@netP$pathways),total_weight=sum(cc@net$weight)),file.path(out,"summary.csv"))
cat("DONE",region,age,ncol(norm),nrow(lr),length(cc@netP$pathways),"
")