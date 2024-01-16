#1.LINCS_BRCA_VIPER
#1.1 imput data
exp <- rio::import("../1.data process/Drug/LINCS_MCF7_drug.csv")
rownames(exp) <- exp$V1
exp <- exp[,-1]
exp <- t(exp)
exp <- as.data.frame(exp)
library('org.Hs.eg.db')
#1.2 Map symbol to entrezid
entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(exp), columns = 'ENTREZID', keytype ='SYMBOL' )
entrezid <- entrezid[!is.na(entrezid$ENTREZID),]
idx <- intersect(rownames(exp),entrezid$SYMBOL)
exp <- exp[rownames(exp)%in%idx,]
exp$symbol <- rownames(exp)
exp <- merge(exp,entrezid,by.x="symbol",by.y="SYMBOL")
rownames(exp) <- exp$ENTREZID
exp <- exp[,2:6232]

#1.3 Viper data preparation
exp1 <- lapply(exp, as.numeric)
exp2 <- as.data.frame(exp1)
rownames(exp2) <- rownames(exp)
colnames(exp2) <- colnames(exp)
exprs <- as.matrix(exp2)

library("convert")
minimalSet <- ExpressionSet(assayData = exprs) 

library(aracne.networks)
data(package="aracne.networks")$results[, "Item"]
data(regulonbrca)

#1.4 viper
library(viper)
vpres <- viper(minimalSet, regulonbrca, verbose = FALSE)
dim(vpres)
data_drug <- vpres@assayData[["exprs"]]
data_drug <- as.data.frame(data_drug)

#1.5 Map entrezid to symbol
entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(data_drug), columns = 'SYMBOL', keytype ='ENTREZID' )
entrezid <- entrezid[!is.na(entrezid$SYMBOL),]
idx <- intersect(rownames(data_drug),entrezid$ENTREZID)
data_drug <- data_drug[rownames(data_drug)%in%idx,]
data_drug$ENTREZID <- rownames(data_drug)

data_drug <- merge(data_drug,entrezid,by="ENTREZID")
rownames(data_drug) <- data_drug$SYMBOL
data_drug <- data_drug[,2:6232]
write.csv(data_drug,"Viper/LINCS_BRCA_VIPER_Drug.csv")

#2.TCGA_BRCA_VIPER
load("../1.data process/Patient/TCGA_BRCA_exp.Rdata")
exp <- as.data.frame(scale(exp))
# 2.1 Map symbol to entrezid
entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(exp), columns = 'ENTREZID', keytype ='SYMBOL' )
entrezid <- entrezid[!is.na(entrezid$ENTREZID),]
idx <- intersect(rownames(exp),entrezid$SYMBOL)
exp <- exp[rownames(exp)%in%idx,]
exp$symbol <- rownames(exp)
exp <- merge(exp,entrezid,by.x="symbol",by.y="SYMBOL")
rownames(exp) <- exp$ENTREZID
exp <- exp[,2:1134]

#2.2 Viper data preparation
exp <- as.data.frame(scale(exp))
exp1 <- lapply(exp, as.numeric)
exp2 <- as.data.frame(exp1)
rownames(exp2) <- rownames(exp)
colnames(exp2) <- colnames(exp)

exprs <- as.matrix(exp2)
minimalSet <- ExpressionSet(assayData = exprs) 

library(aracne.networks)
data(package="aracne.networks")$results[, "Item"]
data(regulonbrca)

#2.3 Viper
vpres <- viper(minimalSet, regulonbrca, verbose = FALSE)
dim(vpres)
data_patient <- vpres@assayData[["exprs"]]
data_patient <- as.data.frame(data_patient)

#2.5 Map entrezid to symbol
entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(data_patient), columns = 'SYMBOL', keytype ='ENTREZID' )
entrezid <- entrezid[!is.na(entrezid$SYMBOL),]
idx <- intersect(rownames(data_patient),entrezid$ENTREZID)
data_patient <- data_patient[rownames(data_patient)%in%idx,]
data_patient$ENTREZID <- rownames(data_patient)
data_patient <- merge(data_patient,entrezid,by="ENTREZID")
rownames(data_patient) <- data_patient$SYMBOL
data_patient <- data_patient[,2:1134]

write.csv(data_patient,file = "Viper/TCGA_BRCA_VIPER_Patient.csv")

#3.Patient_Viper_data
exp <- rio::import("Viper/TCGA_BRCA_VIPER_Patient.csv")
rownames(exp) <- exp$V1
exp <- exp[,-1]

exp <- as.data.frame(scale(exp))
#3.1 group_list
library(tidyverse)
table(str_sub(colnames(exp),14,15))
group_list = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
group_list = factor(group_list,levels = c("normal","tumor"))
table(group_list)

#3.2 limma analysis
library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
patient_sig <- deg[,1,drop=F]

#3.3 Patient data value range (-1,1)
# Function to scale scores
.S <- function(scores) {
  p <- max(scores)
  q <- min(scores)
  ifelse(scores == 0, 0, ifelse(scores > 0, scores / p, -scores / q))
}

scaled_score <- .S(patient_sig[,1])
scaled_score <- as.data.frame(scaled_score)
rownames(scaled_score) <- rownames(patient_sig)
scaled_score_patient <- scaled_score


#4. Viper_Drug_data
drug <- rio::import("Viper/LINCS_BRCA_VIPER_Drug.csv")
rownames(drug) <- drug$V1
drug <- drug[,-1]

#4.1 Drug data value range (-1,1)
data_drug <- list()
count=0
for (i in 1:length(drug)) {
  count=count+1
  print(count)
  x <- drug[i]
  scaled_score <- .S(x[,1])
  scaled_score <- as.data.frame(scaled_score)
  names(scaled_score) <-colnames(x)
  data_drug[[i]] <- scaled_score
  
}
data_drug_score <- do.call(cbind,data_drug)
rownames(data_drug_score) <- rownames(drug)
scaled_score_drug <- data_drug_score
#5.save
save(scaled_score_patient,scaled_score_drug,file = "Viper/BRCA_drug_patient_Viper.Rdata")



load("Viper/BRCA_drug_patient_Viper.Rdata")

#6.RS
code_dir <- "RS-master/"
cancer <- "BRCA"
#core functions
source(paste(code_dir, "core_functions.R", sep=""))
library("ROCR")

lincs_signatures <- scaled_score_drug
gene.list <- rownames(lincs_signatures)
sig.ids <- colnames(scaled_score_drug)
dz_signature <- scaled_score_patient

#Up
dz_signature$Protein <- rownames(dz_signature)
dz_signature_top <- dz_signature[order(dz_signature$scaled_score,decreasing = TRUE),]
dz_signature_top$up_down <- "up"
dz_signature_top$up_down[dz_signature_top$scaled_score<0] <- "down"
dz_signature_top$Protein<- rownames(dz_signature_top)
dz_genes_up <- subset(dz_signature_top,up_down=="up",select="Protein")

#Down
dz_signature_bottom <- dz_signature[order(dz_signature$scaled_score),]
dz_signature_bottom$up_down <- "up"
dz_signature_bottom$up_down[dz_signature_bottom$scaled_score<0] <- "down"
dz_signature_bottom$Protein<- rownames(dz_signature_bottom)
dz_genes_down <- subset(dz_signature_bottom,up_down=="down",select="Protein")

#compute RS
#only choose the top 200 genes
max_gene_size <- 125
if (nrow(dz_genes_up) > max_gene_size){
  dz_genes_up <- data.frame(GeneID= dz_genes_up[1:max_gene_size,])
}
if (nrow(dz_genes_down) > max_gene_size){
  dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,])
}
dz_genes_down$up_down <- "down"
dz_genes_up$up_down <- "up"

dz_cmap_scores <- NULL
count <- 0
for (exp_id in sig.ids) {
  count <- count + 1
  print(count)
  cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
  colnames(cmap_exp_signature) <- c("ids","rank")
  dz_cmap_scores <- c(dz_cmap_scores, cmap_score_new(dz_genes_up,dz_genes_down,cmap_exp_signature))
}

#random scores
N_PERMUTATIONS <- 10000 #default 100000
random_sig_ids <- sample(1:ncol(lincs_signatures),N_PERMUTATIONS,replace=T)
count <- 0
random_cmap_scores <- NULL
for (expr_id in random_sig_ids){
  count <- count + 1
  print(count)
  cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
  colnames(cmap_exp_signature) <- c("ids","rank")
  
  random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
  rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
  rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
  random_cmap_scores <- c(random_cmap_scores, cmap_score_new(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature))
}

p <- sapply(dz_cmap_scores, function(score){
  sum(random_cmap_scores < score)/length(random_cmap_scores)
})

padj <- p.adjust(p, "fdr")
results <- data.frame(id = sig.ids, RS = dz_cmap_scores, p, padj)
results <- results[order(results$RS),]

#save
save(results,file = "Viper/BRCA_results_Viper.Rdata")
write.csv(results,file = "Viper/BRCA_score_Viper_all.csv")

#7 RS——IC50
#7.1.导入all_drug
data <- rio::import("Viper/BRCA_score_Viper_all.csv")
library(dplyr)
data <- data[,-1]
data$id <- toupper(data$id)

#7.2 Convert drug names and get intersection
drug <- rio::import("Data/MCF7_24_10_ID_NEW.csv")
drug$pert_iname <- toupper(drug$pert_iname)
drug <- drug[,c(1,2)]
idx <- intersect(data$id,drug$pert_iname) 
exp1 <- merge(data,drug,by.x = "id",by.y = "pert_iname")
exp1 <- exp1[!duplicated(exp1$id),]
exp2 <- data[!data$id%in%idx,]
exp2$pert_id <- exp2$id
exp <- rbind(exp1,exp2)

#7.3 Import IC50 from CHEMBL database
ic50 <- rio::import("Data/MCF7_IC50.csv")
ic50$id <- toupper(ic50$id)
idy <- intersect(exp$id,ic50$id) 
drug_ic50 <- merge(exp,ic50,by="id")
names(drug_ic50)[2] <- "RS"
drug_ic50_all <- subset(drug_ic50,select=c("id","RS","standard_value"))
save(drug_ic50_all,file = "Viper/BRCA_drug_ic50_all_Viper.Rdata")
write.csv(drug_ic50_all,file = "Viper/BRCA_drug_ic50_all_Viper.csv")

#7.4 Import FDA drugs
fda <- rio::import("Data/FDA.csv")
idn <- intersect(exp$id,fda$pert_iname) 
fda_all_RS <- merge(exp,fda,by.x = "id",by.y = "pert_iname")
write.csv(fda_all_RS,file = "Viper/BRCA_MCF7_Drug_FDA_RS_Viper.csv")
save(fda_all_RS,file = "Viper/BRCA_fda_all_RS_Viper.Rdata")


#8 RS-IC50
load("Viper/BRCA_drug_ic50_all_Viper.Rdata")
library(ggplot2)
library(cowplot)

cancer <- "BRCA"
drug_activity_RS <- drug_ic50_all
drug_activity_RS <- subset(drug_activity_RS,select=c("id","RS","standard_value"))
plot(drug_activity_RS$RS, log(drug_activity_RS$standard_value, 10))
drug_activity_RS <- drug_activity_RS[order(drug_activity_RS$RS),]
plot(drug_activity_RS$RS, log(drug_activity_RS$standard_value, 10))
cor_test <- cor.test(drug_activity_RS$RS, log(drug_activity_RS$standard_value, 10),method="spearman", exact=FALSE)
cor_test
lm_cmap_ic50 <- lm(RS ~ log(standard_value, 10), drug_activity_RS)
drug_activity_RS$activity <- "Effective"
drug_activity_RS$activity[drug_activity_RS$standard_value>10000] <- "Ineffective"
efficacy_test <- t.test(drug_activity_RS$RS[drug_activity_RS$activity == "Effective"], drug_activity_RS$RS[drug_activity_RS$activity == "Ineffective"])


pdf(paste( "fig/", cancer, "_TCGA_LINCS_MCF7_Viper_RS-IC50", ".pdf", sep=""))

p <- ggplot(drug_activity_RS, aes(RS, log(drug_activity_RS$standard_value, 10)  )) +  theme_bw()  +
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) +  geom_point(aes(color = activity), size=3) +  # 在这里设置颜色
  scale_color_manual(values = c("Effective" = "#DC756D", "Ineffective" = "#6BBEC3"),guide = "none") +
  annotate("text", label = paste(cancer, ",", "MCF7", ",","TCGA",sep=""),
           x = -0.1, y = 7.5, size = 6, colour = "black") +
  annotate("text", label = paste("r=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), 
           x = -0.1, y = 7, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab(expression("RS"[Viper])) + guides(shape="none", size="none") +
  ylab(expression(paste("log"[10], "(IC"[50], ") nm"))) + coord_cartesian(xlim = c(-1.3, 1), ylim=c(0, 7.5))
bar_plot =  ggplot(drug_activity_RS, aes(activity, RS,colour  = activity  )) + geom_boxplot() +  coord_flip()+
  xlab("") + ylab(expression("RS"[Viper])) +
  scale_color_manual(values = c("Effective" = "#DC756D", "Ineffective" = "#6BBEC3"),guide = "none")+
  annotate("text", label = paste("P=", format(efficacy_test$p.value, digit=3, scientific=T), sep=""), 
           x = 1.5, y = 0.8, size = 6, colour = "black") + theme(legend.position ="none")+
  theme_bw()+theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18)) 

ggdraw() +
  draw_plot(p, 0, .25, 1, .75) +
  draw_plot(bar_plot, 0, 0, 1, .25) 

dev.off()
