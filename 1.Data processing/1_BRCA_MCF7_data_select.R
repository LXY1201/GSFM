#compare expressions between TCGA tumors  and cell lines; detect inconsistent samples; visualize the samples

library("RColorBrewer")
library("gplots")
library("ggplot2")
library("DESeq2")
library("RMySQL")
library("ROCR")
library(sparcl)

cancer="BRCA"
cancer_name = 'Breast cancer'
comparison_gene_set = "varying5k"
cutoff = 0.05
num_varying_genes = 5000
cell="MCF7"
proj = "TCGA-BRCA"
########
#functions
########
is_outlier <- function(cancer, cell_line_one_tumor_anno){
  #compute correlation between tumors and dz related cell line as well as other cell lines
  ccle_tcga = read.csv("ccle/ccle_cellline_tcga_mapping_updated.csv", stringsAsFactors=F)
  ccle_tcga_target_celllines = ccle_tcga$CCLE.name[ccle_tcga$tcga.tumor == cancer]
  cor_same = cell_line_one_tumor_anno$cor[cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines]
  cor_others = cell_line_one_tumor_anno$cor[!(cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines)]
  p = wilcox.test(cor_same, cor_others, alternative = 'greater')  
  return(p$p.value)
}

is_outlier3 <- function(cancer, cell_line_one_tumor_anno){
  ccle_tcga = read.csv("ccle/ccle_cellline_tcga_mapping_updated.csv", stringsAsFactors=F)
  ccle_tcga_target_celllines = ccle_tcga$CCLE.name [ccle_tcga$tcga.tumor == cancer]
  cor_same = cell_line_one_tumor_anno$cor[cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines]
  cor_others = cell_line_one_tumor_anno$cor[!(cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines)]
  pred <- prediction(c(cor_same, cor_others), c(rep(T, length(cor_same)), rep(F, length(cor_others))))
  perf <- performance(pred, measure = "auc") 
  
  return(perf@y.values[[1]])
}

bean.plot <- function(tumor_cell_cor, cancer, type){
  group = sapply(1:nrow(tumor_cell_cor), function(id){
    if (tumor_cell_cor[id, "histology_class"] == "others"){
      paste(tumor_cell_cor[id, "sample_id"], 1)
    }else{
      paste(tumor_cell_cor[id, "sample_id"], 2)
    }
  })
  
  tumor_cell_cor$group = group
  
  pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", cancer, "_", type, "_beanplot.pdf",sep=""))
  par(mar = c(10, 4, 4, 2) + 0.1, las=2)
  beanplot(cor ~ group, data = tumor_cell_cor, ll = 0.15,
           main = "", ylab = "correlation", side = "both",
           border = NA, col = list(c("grey","green","green"), c("grey", "red","red")), 
           beanlines = "median", overalllin="median", method="jitter")
  legend("bottomleft", fill = c("green", "red"),
         legend = c("others", paste(cancer_name,"related")))
  dev.off()
}

box.plot <- function(tumor_cell_cor, cancer, type){
  pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", cancer, "_", type, "boxplot.pdf",sep=""))
  #random select 5 samples
  p <- ggplot(tumor_cell_cor, aes(factor(sample_id), cor, color= histology_class))
  print(p +   geom_boxplot(outlier.colour = "grey", notch=F, outlier.shape = NA) + geom_jitter() +   
          ylab("correlation") +
          xlab("sample") +  theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

plot.sample <- function(tumor_cell_all_subset, cancer, sample1, type){
  dir.create(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", "samples_", cancer,  "_", type, sep=""), showWarnings = FALSE)
  tumor_cell_all_outliers_sample = subset(tumor_cell_all_subset, as.character(sample_id) == sample1)
  top_tumor_types = as.character(get.top.cor.tumor.type(tumor_cell_all_outliers_sample, sample1))
  
  tumor_cell_cor= subset(tumor_cell_all_outliers, tumor_type %in% unique(c(cancer, top_tumor_types)))
  
  #order based on median cor
  tumor_cell_cor_merged = aggregate(cor ~ tumor_type_name, tumor_cell_cor, median)
  tumor_cell_cor_merged = tumor_cell_cor_merged[order(tumor_cell_cor_merged$cor, decreasing=T), ]
  tumor_cell_cor$tumor_type_name = factor(tumor_cell_cor$tumor_type_name, levels = tumor_cell_cor_merged$tumor_type_name)
  
  pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", "samples_", cancer, "_", type, "/", sample1, ".pdf",sep=""))
  
  p <- ggplot(tumor_cell_cor, aes(tumor_type_name, cor))
  print(p +   geom_boxplot(outlier.colour = "grey", notch=F, outlier.shape = NA) + geom_jitter() +   
          ylab("correlation") +
          xlab("") +  theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

get.top.cor.tumor.type <- function(tumor_cell_all_outliers_sample, sample1){
  tumor_types = unique(tumor_cell_all_outliers_sample$tumor_type)
  
  p_values = NULL
  for (tumor_type in tumor_types){
    
    ccle_tcga_target_celllines = tumor_cell_all_outliers_sample$CCLE.name[tumor_cell_all_outliers_sample$tumor_type == tumor_type]
    cor_same = tumor_cell_all_outliers_sample$cor[tumor_cell_all_outliers_sample$CCLE.name %in% ccle_tcga_target_celllines]
    cor_others = tumor_cell_all_outliers_sample$cor[!(tumor_cell_all_outliers_sample$CCLE.name %in% ccle_tcga_target_celllines)]
    p = wilcox.test(cor_same, cor_others, alternative = 'greater')  
    p_values = c(p_values, p$p.value)
  }
  
  tumor_type_p = data.frame(tumor_type = tumor_types, p = p_values)
  tumor_type_p = tumor_type_p[order(tumor_type_p$p),]
  return (tumor_type_p$tumor_type[1:5])
}


############################
#MAIN
###########################
###########################

###########
##process RNASEQ from TCGA
#quality control
#cross validation; validate with independent sets
#1.data_from_gdac
raw.data.header = read.csv(paste('raw/gdac/',cancer,'/rnaseq/gdac.broadinstitute.org_', cancer, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2014051800.0.0/', cancer, '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep=""), sep='\t', nrow=1, check.names=F)

raw.data = read.csv(paste('raw/gdac/',cancer,'/rnaseq/gdac.broadinstitute.org_', cancer, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2014051800.0.0/', cancer, '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep=""), sep='\t', skip=1, check.names=F)
genes = as.character(raw.data[,1])
GeneID = sapply(genes, function(id){
  as.character(unlist(strsplit(id, "\\|"))[2])
})


countTable = raw.data[,colnames(raw.data) == 'scaled_estimate']
raw.data.header = colnames(raw.data.header)[colnames(raw.data) == 'scaled_estimate']
colnames(countTable) = raw.data.header 
rownames(countTable) = GeneID


if (!file.exists(cancer)){
  dir.create(cancer)
}
if (!file.exists(paste(cancer, "/tumor_cell_line", sep=""))){
  dir.create(paste(cancer, "/tumor_cell_line", sep=""))
}

##compare with expression in ccle
load("ccle/ccle_meta_updated.RData")
load("ccle/ccle_gene_updated.RData")

iqr_gene = apply(ccle_gene[,-c(1)], 1, IQR)
gene_ids = ccle_gene$GeneID


varying_genes = gene_ids[(order(iqr_gene, decreasing=T))][1:num_varying_genes] 
selected_genes = varying_genes

ccle_gene = subset(ccle_gene, GeneID %in% selected_genes)    

tumor_sig = data.frame(GeneID, countTable)
cell_line_tumor = merge( ccle_gene, tumor_sig, by="GeneID")
cell_line_tumor_cor = cor(cell_line_tumor[,-c(1)], method="spearman")

tumors = colnames(tumor_sig)[-1]
celllines = colnames(ccle_gene)[-1]

tumor_outliers = NULL
tumor_cell_all = data.frame()
for (tumor in tumors){
  cell_line_one_tumor_cor = cell_line_tumor_cor[tumor, celllines ] 
  cell_line_one_tumor_cor = data.frame(sample = names(cell_line_one_tumor_cor), cor=cell_line_one_tumor_cor)
  cell_line_one_tumor_anno = merge(ccle_meta, cell_line_one_tumor_cor, by.y="sample", by.x="CCLE.name")
  cell_line_one_tumor_anno = subset(cell_line_one_tumor_anno, select=c("Cell.line.primary.name", "CCLE.name", "cor", "Histology","Hist.Subtype1" ,"Site.Primary"))  
  cell_line_one_tumor_anno$barcode = paste(unlist(strsplit(tumor, '\\.')), collapse="-")
  tags = unlist(strsplit(as.character(tumor), "\\."))
  cell_line_one_tumor_anno$patient_id =   paste(tags[1:3], collapse="-")
  cell_line_one_tumor_anno$sample_id =   paste(tags[1:4], collapse="-")
  #tumor_outliers = c(tumor_outliers, is_outlier(cancer, cell_line_one_tumor_anno))
  cell_line_one_tumor_anno$outlier = is_outlier(cancer, cell_line_one_tumor_anno)
  tumor_cell_all = rbind(tumor_cell_all, cell_line_one_tumor_anno)
}

tumor_cell_all = tumor_cell_all[order(tumor_cell_all$outlier,tumor_cell_all$barcode, tumor_cell_all$cor, decreasing=T),]

if (!file.exists(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", sep=""))){
  dir.create(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", sep=""))
}
write.csv(tumor_cell_all, paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", "tumor_cell_all", cancer, ".csv", sep=""))

###
#detect tumor samples that are not correlated to cancer cell lines 
tumor_cell_all_p = aggregate(outlier ~ patient_id + sample_id, tumor_cell_all, min)
tumor_cell_all_p_adj = p.adjust(tumor_cell_all_p$outlier, "fdr")
#cutoff: by default 0.05; could use the cutoff by randomly sampling; if so, I think we don't need to adjust p value
patient_outliers =  unique(tumor_cell_all_p$patient_id[tumor_cell_all_p$outlier > cutoff])
sample_outliers =  unique(tumor_cell_all_p$sample_id[tumor_cell_all_p$outlier > cutoff])

write(patient_outliers, paste(cancer,"/tumor_cell_line/", comparison_gene_set, "/", "outlier_", cancer, "_outlier.txt", sep=""), ncolumn=1)

#visualize outliers
ccle_tcga = read.csv("ccle/ccle_cellline_tcga_mapping_updated.csv")
ccle_tcga_target_celllines = ccle_tcga$CCLE.name[ccle_tcga$tcga.tumor == cancer]

histology_class <- sapply(tumor_cell_all$CCLE.name, function(cell_line){
  if (cell_line %in% ccle_tcga_target_celllines){
    paste(cancer, "related")
  }else{
    "others"
  }
})

tumor_cell_all$histology_class <-  histology_class

tumor_cell_all_outliers = tumor_cell_all[ as.character(tumor_cell_all$sample_id) %in% sample(as.character(sample_outliers), min(10, length(sample_outliers))),]
bean.plot(tumor_cell_all_outliers, cancer, "bad")
#box.plot(tumor_cell_all_outliers, cancer, "bad")

good_samples = sample(unique(tumor_cell_all$sample_id[tumor_cell_all$outlier < cutoff]), min(10, length(unique(tumor_cell_all$sample_id[tumor_cell_all$outlier < cutoff]))))
tumor_cell_all_good = tumor_cell_all[ as.character(tumor_cell_all$sample_id) %in% as.character(good_samples),]
bean.plot(tumor_cell_all_good, cancer, "good")
#box.plot(tumor_cell_all_good, cancer, "good")

#find top 5 enriched cell type
tumor_cell_mapping = read.csv("ccle/ccle_cellline_tcga_mapping_updated.csv")
tumor_cell_mapping = subset(tumor_cell_mapping, select=c("CCLE.name", "tcga.tumor", "tumor_type_name"))
names(tumor_cell_mapping) = c("CCLE.name", "tumor_type", "tumor_type_name")

tumor_cell_all_outliers = merge(tumor_cell_all_outliers, tumor_cell_mapping,  by.x="CCLE.name")
samples = unique(tumor_cell_all_outliers$sample_id)
for(sample1 in samples){
  plot.sample(tumor_cell_all_outliers, cancer, sample1, "bad")
}

tumor_cell_all_good = merge(tumor_cell_all_good, tumor_cell_mapping,  by.x="CCLE.name")
samples = unique(tumor_cell_all_good$sample_id)
for(sample1 in samples){
  plot.sample(tumor_cell_all_good, cancer, sample1, "good")
}

############
#clinical features
#use raw count instead of scaled estimate, 
raw.data.header = read.csv(paste('raw/gdac/',cancer, '/rnaseq/gdac.broadinstitute.org_', cancer, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2014051800.0.0/', cancer, '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep=""), sep='\t', nrow=1, check.names=F)
raw.data = read.csv(paste('raw/gdac/',cancer, '/rnaseq/gdac.broadinstitute.org_', cancer, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2014051800.0.0/', cancer, '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep=""), sep='\t', skip=1, check.names=F)
countTable = raw.data[,colnames(raw.data) == 'raw_count']
countTable = raw.data[,colnames(raw.data) == 'raw_count']
raw.data.header = colnames(raw.data.header)[colnames(raw.data) == 'raw_count']
colnames(countTable) = raw.data.header 
rownames(countTable) = GeneID

table(tumor_cell_all$Cell.line.primary.name)
tumor_cell_all_cancer <- tumor_cell_all[tumor_cell_all$Cell.line.primary.name == cell,]
tumor_cell_all_cancer <- tumor_cell_all_cancer[tumor_cell_all_cancer$outlier< cutoff,]

library(stringr)
dat = read.table(paste0(proj,".htseq_counts.tsv.gz"),check.names = F,row.names = 1,header = T)
tumor_cell_all_BRCA <- tumor_cell_all_cancer[tumor_cell_all_cancer$sample_id%in%colnames(dat),]
sample_id <- tumor_cell_all_BRCA$sample_id
patient_id <- tumor_cell_all_BRCA$patient_id
tumor_type <- ifelse(as.numeric(str_sub(tumor_cell_all_BRCA$sample_id,14,15)) < 10,'tumor','non-tumor')
tcga_barcode <-tumor_cell_all_BRCA$barcode 

data <- dat[,colnames(dat)%in%tumor_cell_all_BRCA$sample_id]

save(dat,file=paste(cancer,"/TCGA_BRCA_dat.Rdata", sep=""))

