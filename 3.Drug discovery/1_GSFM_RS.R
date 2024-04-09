#1.patient_GSFM
exp <- rio::import("../../DOCK/project/Sample_output/BRCA/TCGA_BRCA_patient_GSFM.csv")
rownames(exp) <- exp$V1
exp <- exp[,-1]
exp <- as.data.frame(t(scale(exp)))
#1.1 group
library(tidyverse)
table(str_sub(colnames(exp),14,15))
group_list = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
group_list = factor(group_list,levels = c("normal","tumor"))
table(group_list)

#1.2  limma
library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
patient_sig <- deg[,1,drop=F]
patient_sig[is.na(patient_sig)] <- 0

#1.3 rank (-1,1)
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

#2.drug_GSFM
drug <- rio::import("../../DOCK/project/Sample_output/BRCA/LINCS_MCF7_drug_GSFM.csv")
rownames(drug) <- drug$V1
drug <- drug[,-1]
drug <- as.data.frame(t(drug))

#2.1 rank (-1,1)
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
#3.save
save(scaled_score_patient,scaled_score_drug,file = "GSFM/BRCA_drug_patient_GSFM.Rdata")

#4.RS
load("GSFM/BRCA_drug_patient_GSFM.Rdata")
code_dir <- "RGES-master/"
cancer <- "BRCA"
#core functions

source(paste(code_dir, "core_functions.R", sep=""))
library("ROCR")

lincs_signatures <- scaled_score_drug
gene.list <- rownames(lincs_signatures)
sig.ids <- colnames(scaled_score_drug)
dz_signature <- scaled_score_patient

dz_signature$GSFM <- rownames(dz_signature)
dz_signature_top <- dz_signature[order(dz_signature$scaled_score,decreasing = TRUE),]
dz_signature_top$up_down <- "up"
dz_signature_top$up_down[dz_signature_top$scaled_score<0] <- "down"
dz_signature_top$GSFM<- rownames(dz_signature_top)
dz_genes_up <- subset(dz_signature_top,up_down=="up",select="GSFM")

dz_signature_bottom <- dz_signature[order(dz_signature$scaled_score),]
dz_signature_bottom$up_down <- "up"
dz_signature_bottom$up_down[dz_signature_bottom$scaled_score<0] <- "down"
dz_signature_bottom$GSFM<- rownames(dz_signature_bottom)
dz_genes_down <- subset(dz_signature_bottom,up_down=="down",select="GSFM")

#compute RS
#only choose the top 200 genes
max_gene_size <- 125
if (nrow(dz_genes_up) > max_gene_size){
  dz_genes_up <- data.frame(GSFM= dz_genes_up[1:max_gene_size,])
}
if (nrow(dz_genes_down) > max_gene_size){
  dz_genes_down <- data.frame(GSFM=dz_genes_down[1:max_gene_size,])
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
results <- data.frame(id = sig.ids, cmap_score = dz_cmap_scores, p, padj)
results <- results[order(results$cmap_score),]

results$cmap_score <- .S(results$cmap_score)
#4. save BRCA_HALLMARK_RS
save(results,file = "GSFM/BRCA_results_GSFM.Rdata")

write.csv(results,file = "GSFM/BRCA_score_GSFM_all.csv")


#5.RS—IC50
#5.1 import all_drug
data <- rio::import("GSFM/BRCA_score_GSFM_all.csv")#6231
library(dplyr)
data <- data[,-1]
data$id <- toupper(data$id)

#5.2
drug <- rio::import("Data/MCF7_24_10_ID_NEW.csv")
drug$pert_iname <- toupper(drug$pert_iname)
drug <- drug[,c(1,2)]

#5.3
idx <- intersect(data$id,drug$pert_iname) 
exp1 <- merge(data,drug,by.x = "id",by.y = "pert_iname")
exp1 <- exp1[!duplicated(exp1$id),]
exp2 <- data[!data$id%in%idx,]
exp2$pert_id <- exp2$id
exp <- rbind(exp1,exp2)
write.csv(exp,"out/BRCA_MCF7_allnames_GSFM.csv")

#5.4 import IC50 from CHEMBL
ic50 <- rio::import("Data/MCF7_IC50.csv")
ic50$id <- toupper(ic50$id)
idy <- intersect(exp$id,ic50$id) #86
drug_ic50 <- merge(exp,ic50,by="id")
names(drug_ic50)[2] <- "RS"
drug_ic50_all <- subset(drug_ic50,select=c("id","RS","standard_value"))
save(drug_ic50_all,file = "GSFM/BRCA_drug_ic50_all_GSFM.Rdata")
write.csv(drug_ic50_all,file = "GSFM/BRCA_drug_ic50_all_GSFM.csv")


#5.5  import FDA drug
fda <- rio::import("Data/FDA.csv")
fda$pert_iname <- toupper(fda$pert_iname)
idz <- intersect(fda$pert_iname,data$id) 
drug_fda <- merge(drug_ic50,fda,by.x = "id",by.y = "pert_iname")
names(drug_fda)[2] <- "RS"
drug_fda <- drug_fda[order(drug_fda$RS),]
save(drug_fda,file = "GSFM/BRCA_drug_fda_GSFM.Rdata")
write.csv(drug_fda,file = "GSFM/BRCA_drug_fda_GSFM.csv")


#5.6fda-drug-all_RS
idn <- intersect(exp$id,fda$pert_iname) 
fda_all_RS <- merge(exp,fda,by.x = "id",by.y = "pert_iname")
names(fda_all_RS)[2] <- "RS"
fda_all_RS <- fda_all_RS[order(fda_all_RS$RS),]
write.csv(fda_all_RS,file = "GSFM/BRCA_MCF7_Drug_FDA_RS_GSFM.csv")
save(fda_all_RS,file = "GSFM/BRCA_fda_all_RS_GSFM.Rdata")


#6 corrlation
load("GSFM/BRCA_drug_ic50_all_GSFM.Rdata")
library(ggplot2)
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

library(cowplot)


pdf(paste( "fig/", cancer, "_TCGA_LINCS_MCF7_GSFM_RS-IC50",".pdf", sep=""))

p <- ggplot(drug_activity_RS, aes(RS, log(drug_activity_RS$standard_value, 10)  )) +  theme_bw()  +
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) +  
  geom_point(aes(color = activity), size=3) +  # 在这里设置颜色
  scale_color_manual(values = c("Effective" = "#DC756D", "Ineffective" = "#6BBEC3"),guide = "none")+ 
  annotate("text", label = paste(cancer, ",", "MCF7", ",","TCGA",sep=""),
           x = -0.65, y = 7.5, size = 6, colour = "black") +

  annotate("text", label = paste("r=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), 
           x = -0.65, y = 7, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab(expression("RS"[GSFM])) + guides(shape="none", size="none") +
  ylab(expression(paste("log"[10], "(IC"[50], ") nm"))) + coord_cartesian(xlim = c(-1, -0.35), ylim=c(0, 7.5))
bar_plot =  ggplot(drug_activity_RS, aes(activity, RS,colour  = activity  )) + geom_boxplot() +  coord_flip()+
  xlab("") + ylab(expression("RS"[GSFM])) +
  scale_color_manual(values = c("Effective" = "#DC756D", "Ineffective" = "#6BBEC3"),guide = "none")+
  annotate("text", label = paste("P=", format(efficacy_test$p.value, digit=3, scientific=T), sep=""), 
           x = 1.5, y = -0.25, size = 6, colour = "black") + theme(legend.position ="none")+
  theme_bw()+theme(legend.position ="bottom", 
                   axis.text=element_text(size=18), 
                   axis.title=element_text(size=18),
                   axis.text.x = element_text(color = "black")) 

ggdraw() +
  draw_plot(p, 0, .25, 1, .75) +
  draw_plot(bar_plot, 0, 0, 1, .25) 

dev.off()


