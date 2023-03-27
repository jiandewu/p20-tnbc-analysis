# source("http://bioconductor.org/biocLite.R")
# biocLite("rgl")
# biocLite("limma")
# biocLite("edgeR")
library(limma)
library("rgl")
library(gplots)
# library(edgeR)


setwd(Sys.getenv("BLIA_BLIS_DATA_DIR", unset = "."))
hdat <- read.csv(file = "BLIA_Adipose.txt", row.names = 1, sep = "\t")
pheno <- read.csv(file = "BLIA_Adipose_pheno.txt", sep = "\t")

#---hdat <- apply(hdat, 2, function(x) log2(x+1))

# hdat=removeBatchEffect(hdat,batch=c("0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1"
# ))
#---boxplot(hdat)
#---colm <- colMeans(hdat)
#---diff <- mean(colm[1:54])-mean(colm[55:316])
#---hdat[,c(1:54)] <- hdat[,c(1:54)]-diff
# boxplot(hdat)


rownames(pheno) <- pheno$sample
colnames(hdat) <- rownames(pheno)
rownames(pheno) == colnames(hdat)

##### PCA
sample.groups <- c(rep("red", 54), rep("blue",262))
pc <- princomp(hdat)
#---plot(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", main = "PCA")
#---plot3d(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", zlab="Component 3", type = "s", size = 1.5)

pc$loadings

samples <- factor(pheno$class)
design <- model.matrix(~0 + samples)
colnames(design) <- levels(samples)
design

fit <- lmFit(hdat, design)
names(fit)

fit1 <- contrasts.fit(fit, c(-1,1))


fit1 <- eBayes(fit1)
volcanoplot(fit1,  main = "BLIA vs Adipose")

write.table(topTable(fit1, number=dim(fit1)[1]), file = "result/BLIA_Adipose_p.txt", sep="\t", quote = F, row.names=T)
write.table(hdat, file = "result/BLIA_Adipose_exp.txt", sep="\t", quote = F, row.names=T)

nrow(topTable(fit1, coef=1, number=54675, lfc=2, p.value = 0.05))
probeset.list1 <- topTable(fit1, coef=1, number=54675, lfc=2, p.value = 0.05)
head(probeset.list1)
write.table(probeset.list1, file = "result/BLIA_Adipose_lfc2_p0.05.txt", sep="\t", quote = F, row.names=T)

genes <- rownames(probeset.list1)
head(genes)
sign_exp <- hdat[genes,]
head(sign_exp)
write.table(sign_exp, file = "result/BLIA_Adipose_lfc2_p0.05_exp.txt", sep="\t", quote = F, row.names=T)

