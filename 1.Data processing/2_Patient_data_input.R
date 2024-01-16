load("BRCA/TCGA_BRCA_dat.Rdata")
library(tidyverse)
group_list = ifelse(as.numeric(str_sub(colnames(dat),14,15)) < 10,'tumor','normal')
group_list = factor(group_list,levels = c("normal","tumor"))
table(group_list)

#1.Reverse log
dat = as.matrix(2^dat - 1)
dat[1:4,1:4]

#2.Use apply to convert to an integer matrix
exp = apply(dat, 2, as.integer)
exp[1:4,1:4] 
rownames(exp) = rownames(dat)
exp[1:4,1:4]

# 3.Expression matrix row name ID conversion
library(stringr)
head(rownames(exp))
library(AnnoProbe)
rownames(exp) = str_split(rownames(exp),"\\.",simplify = T)[,1]
head(rownames(exp))
re = annoGene(rownames(exp),ID_type = "ENSEMBL");head(re)
library(tinyarray)
exp = trans_array(exp,ids = re,from = "ENSEMBL",to = "SYMBOL")
exp[1:4,1:4]

#4.Filter criteria: pcg gene
exp <- as.data.frame(exp)
load("all_gtf.rdata")
idx <- intersect(pcg,rownames(exp))
exp <- exp[rownames(exp)%in%idx,]
exp <- log(exp+1)

save(exp,file = "Patient/TCGA_BRCA_exp.Rdata")
write.csv(exp,file = "Patient/TCGA_BRCA_patient.csv")

exp1 <- as.data.frame(t(scale(exp)))
write.csv(exp1,file = "../2.data transformation/project/Sample_input/BRCA/TCGA_BRCA_patient.csv")

