#MCF7
#1.loading data
load("lincs.Rdata") 

#2.LINCS_MCF7_drug_data
library(stringr)
datacolnames <- as.data.frame(datacolnames)
datacolnames$drug=trimws(str_split(datacolnames$V1,'__',simplify = T)[,1])
datacolnames$cell=trimws(str_split(datacolnames$V1,'__',simplify = T)[,2])
datacolnames$tct=trimws(str_split(datacolnames$V1,'__',simplify = T)[,3])
table(datacolnames$cell)

colnames(data) <- datacolnames$cell
data1 <-data[colnames(data)=="MCF7"] 
datacolnames1 <- datacolnames[datacolnames$cell=="MCF7",]
colnames(data1) <- datacolnames1$drug

#3.Expression matrix row name ID conversion
library('org.Hs.eg.db')
entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(data1), columns = 'SYMBOL', keytype ='ENTREZID' )

data1$ENTREZID <- rownames(data1)
data2 <- merge(data1,entrezid,by="ENTREZID")
data2 <- data2[!duplicated(data2$SYMBOL),]
data2 <- data2[!is.na(data2$SYMBOL),]
data3 <- data2[,2:6232]
rownames(data3) <- data2$SYMBOL

#4.Filter criteria
load("all_gtf.rdata")
idx <- intersect(pcg,rownames(data3))
data4 <- data3[rownames(data3)%in%idx,]

#Input data is z-score
#The row name is the sample name, and the column name is the gene name.
data4_t <-as.data.frame(t(data4)) 

write.csv(data4_t,file = "Drug/LINCS_MCF7_drug.csv")
write.csv(data4_t,file = "../2.data transformation/project/Sample_input/BRCA/LINCS_MCF7_drug.csv")