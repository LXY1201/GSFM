#1.Import TCGA_BRCA patient matrix
load("../1.data process/Patient/TCGA_BRCA_exp.Rdata")
exp <- as.data.frame(exp)

#1.2 group_list
library(tidyveRSe)
table(str_sub(colnames(exp),14,15))
group_list = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
group_list = factor(group_list,levels = c("normal","tumor"))
table(group_list)

#1.3 limma analysis
library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
patient_sig <- deg[,1,drop=F]

#1.4 Patient data value range (-1,1)
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

#2 Drug data value range (-1,1)
drug <- rio::import("../1.data process/Drug/LINCS_MCF7_drug.csv")
rownames(drug) <- drug$V1
drug <- drug[,-1]
drug <- as.data.frame(t(drug))

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
#save data
save(scaled_score_patient,scaled_score_drug,file = "Gene/BRCA_drug_patient_gene.Rdata")



load("Gene/BRCA_drug_patient_gene.Rdata")
#3.RS
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
dz_signature$Gene <- rownames(dz_signature)
dz_signature_top <- dz_signature[order(dz_signature$scaled_score,decreasing = TRUE),]
dz_signature_top$up_down <- "up"
dz_signature_top$up_down[dz_signature_top$scaled_score<0] <- "down"
dz_signature_top$Gene<- rownames(dz_signature_top)
dz_genes_up <- subset(dz_signature_top,up_down=="up",select="Gene")

#Down
dz_signature_bottom <- dz_signature[order(dz_signature$scaled_score),]
dz_signature_bottom$up_down <- "up"
dz_signature_bottom$up_down[dz_signature_bottom$scaled_score<0] <- "down"
dz_signature_bottom$Gene<- rownames(dz_signature_bottom)
dz_genes_down <- subset(dz_signature_bottom,up_down=="down",select="Gene")

#compute RS
#only choose the top 200 genes
max_gene_size <- 125
if (nrow(dz_genes_up) > max_gene_size){
  dz_genes_up <- data.frame(GeneID= dz_genes_up[1:max_gene_size,])
}
if (nrow(dz_genes_down) > max_gene_size){
  dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,])
}
dz_genes_up_m <- dz_genes_up
dz_genes_down_m <- dz_genes_down

dz_genes_down_m$up_down <- "down"
dz_genes_up_m$up_down <- "up"


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
N_PERMUTATIONS <- 10000 
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

write.csv(results,file = "Gene/BRCA_score_Gene_all.csv")


#4.imput data
data <- rio::import("Gene/BRCA_score_Gene_all.csv")
library(dplyr)
data <- data[,-1]
data$id <- toupper(data$id)

#4.2 Convert drug names and get inteRSection
drug <- rio::import("Data/MCF7_24_10_ID_NEW.csv")
drug$pert_iname <- toupper(drug$pert_iname)
drug <- drug[,c(1,2)]
idx <- inteRSect(data$id,drug$pert_iname) 
exp1 <- merge(data,drug,by.x = "id",by.y = "pert_iname")
exp1 <- exp1[!duplicated(exp1$id),]
exp2 <- data[!data$id%in%idx,]
exp2$pert_id <- exp2$id
exp <- rbind(exp1,exp2)


#4.3 Import IC50 from CHEMBL database
ic50 <- rio::import("Data/MCF7_IC50.csv")
ic50$id <- toupper(ic50$id)
idy <- inteRSect(exp$id,ic50$id) 
drug_ic50 <- merge(exp,ic50,by="id")
drug_ic50_all <- subset(drug_ic50,select=c("id","RS","standard_value"))
save(drug_ic50_all,file = "Gene/BRCA_drug_ic50_all_Gene.Rdata")
write.csv(drug_ic50_all,file = "Gene/BRCA_drug_ic50_all_Gene.csv")

#4.4 Import FDA drugs
fda <- rio::import("Data/FDA.csv")
idn <- inteRSect(exp$id,fda$pert_iname) 
fda_all_RS <- merge(exp,fda,by.x = "id",by.y = "pert_iname")
write.csv(fda_all_RS,file = "Gene/BRCA_MCF7_Drug_FDA_RS_Gene.csv")
save(fda_all_RS,file = "Gene/BRCA_fda_all_RS_Gene.Rdata")


#5 RS-IC50 corrlation
load("Gene/BRCA_drug_ic50_all_Gene.Rdata")
library(ggplot2)
library(cowplot)


cancer <- "BRCA"
drug_activity_RS <- drug_ic50_all
drug_activity_RS <- subset(drug_activity_RS,select=c("id","RS","standard_value"))
plot(drug_activity_RS$RS, log(drug_activity_RS$standard_value, 10))
drug_activity_RS <- drug_activity_RS[order(drug_activity_RS$RS),]
cor_test <- cor.test(drug_activity_RS$RS, log(drug_activity_RS$standard_value, 10),method="spearman", exact=FALSE)
cor_test
lm_cmap_ic50 <- lm(RS ~ log(standard_value, 10), drug_activity_RS)

drug_activity_RS$activity <- "Effective"
drug_activity_RS$activity[drug_activity_RS$standard_value>10000] <- "Ineffective"
efficacy_test <- t.test(drug_activity_RS$RS[drug_activity_RS$activity == "Effective"], drug_activity_RS$RS[drug_activity_RS$activity == "Ineffective"])


pdf(paste( "fig/", cancer, "_TCGA_LINCS_MCF7_Gene_RS-IC50", ".pdf", sep=""))

p <- ggplot(drug_activity_RS, aes(RS, log(drug_activity_RS$standard_value, 10)  )) +  theme_bw()  +
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) +geom_point(aes(color = activity), size=3) +  # 在这里设置颜色
  scale_color_manual(values = c("Effective" = "#DC756D", "Ineffective" = "#6BBEC3"),guide = "none") +
  annotate("text", label = paste(cancer, ",", "MCF7", ",","TCGA",sep=""),
           x = -0.15, y = 7.5, size = 6, colour = "black") +
  
  annotate("text", label = paste("r=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), 
           x = -0.15, y = 7, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab(expression("RS"[Gene])) + guides(shape="none", size="none") +
  ylab(expression(paste("log"[10], "(IC"[50], ") nm"))) + coord_cartesian(xlim = c(-0.65, 0.41), ylim=c(0, 7.5))
bar_plot =  ggplot(drug_activity_RS, aes(activity, RS,colour  = activity  )) + geom_boxplot() +  coord_flip()+
  xlab("") + ylab(expression("RS"[Gene])) +
  scale_color_manual(values = c("Effective" = "#DC756D", "Ineffective" = "#6BBEC3"),guide = "none")+
  annotate("text", label = paste("P=", format(efficacy_test$p.value, digit=3, scientific=T), sep=""), 
           x = 1.5, y = 0.35, size = 6, colour = "black") + theme(legend.position ="none")+
  theme_bw()+theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18)) 

ggdraw() +
  draw_plot(p, 0, .25, 1, .75) +
  draw_plot(bar_plot, 0, 0, 1, .25) 

dev.off()
