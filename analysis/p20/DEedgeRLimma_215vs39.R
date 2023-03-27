library(limma)
library(edgeR)


options(stringsAsFactors = FALSE)

setwd(Sys.getenv("P20_DATA_DIR", unset = "."))
rawdata <- read.table("258samples_v3.txt", header = TRUE, sep="\t")
row.names(rawdata) <- rawdata[,1]
rawdata <- rawdata[,-1]

# delete bilate samples
drops1 <- if (nzchar(Sys.getenv("P20_BILATERAL_EXCLUSIONS_FILE"))) readLines(Sys.getenv("P20_BILATERAL_EXCLUSIONS_FILE")) else character(0)
rawdata <- rawdata[, !(names(rawdata) %in% drops1)] # all samples

data1 <- DGEList(counts = data.matrix(rawdata[-c(1:11),])) # all samples

obeselevel <- rawdata["OBESE",]
obeselevel[obeselevel==1] = "Obese"
obeselevel[obeselevel==0] = "Lean"
obese <- as.factor(obeselevel)
racelevel <- rawdata["AA",]
racelevel[racelevel==1] = "AA"
racelevel[racelevel==0] = "EA"
race <- as.factor(racelevel)
stagelevel <- rawdata["LATE",]
stagelevel[stagelevel==1] = "Stage234"
stagelevel[stagelevel==0] = "Stage1"
stage <- as.factor(stagelevel)

# label concordant with the Vanderbilt algorithm
# Optional discordant-sample IDs are loaded from P20_DISCORDANT_SAMPLES_FILE.
drops2 <- if (nzchar(Sys.getenv("P20_DISCORDANT_SAMPLES_FILE"))) readLines(Sys.getenv("P20_DISCORDANT_SAMPLES_FILE")) else character(0)
classifylevel <- rawdata["OBESE",]
classifylevel[colnames(rawdata) %in% drops2] = "Discordant"
classifylevel[!colnames(rawdata) %in% drops2] = "TNBC"
classify <- as.factor(classifylevel)

data1$samples$obese <- obese
data1$samples$race <- race
data1$samples$stage <- stage
data1$samples$classify <- classify

table(data1$samples$classify,data1$samples$race)
#----------------#
# AA  EA
# Lean   27  50
# Obese 102  75
#----------------#

# Transformations from the raw-scale
seqdepth_mean <- mean(data1$samples$lib.size) * 1e-6
seqdepth_median <- median(data1$samples$lib.size) * 1e-6

summary(data1$samples$lib.size * 1e-6)
#----------------#
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9896 15.9975 19.1382 19.3924 23.1632 36.3612 
#----------------#
nsamples <- ncol(data1)

table(rowSums(data1$counts==0)==nsamples)
#----------------#
# FALSE  TRUE 
# 25511  2767 
#----------------#

# Plotting the distribution log-CPM values shows that 
# a sizeable proportion of genes within each sample are 
# either unexpressed or lowly-expressed with log-CPM values 
# that are small or negative

lcpm.cutoff <- log2(10/seqdepth_median + 2/seqdepth_mean)
library(RColorBrewer)
col <- brewer.pal(nsamples, "Paired")

lcpm <- cpm(data1, log=TRUE, prior.count=2)
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.5), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(data1), text.col=col, bty="n")

# Removing genes that are lowly expressed
keep.exprs <- filterByExpr(data1, group=classify)
data2 <- data1[keep.exprs,, keep.lib.sizes=FALSE]
dim(data2)
## [1] 16912   254

lcpm <- cpm(data2, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.5), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(data1), text.col=col, bty="n")

# Normalising gene expression distributions
par(mfrow=c(1,2))
lcpm <- cpm(data2, log=TRUE, prior.count=2)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="Log-cpm")
data2 <- calcNormFactors(data2, method = "TMM")
# data2$samples$norm.factors
#----------------#
#  [1] 1.1005944 1.0213063 0.9484682 1.0938885 0.9135321 1.1018828 1.0813619 0.9746213 0.9708049 1.0128725
# [11] 0.9602972 1.0116042 1.0162326 1.0180576 1.0575765 0.9646913 1.0454093 1.0034687 0.9942733 1.0220867
# [21] 1.0738821 1.0687922 1.0386444 1.0766185 1.0205136 0.9873789 0.9493396 0.9206837 1.0122073 1.0469196
# [31] 1.0697226 1.0397318 0.9758840 0.9220561 0.9134998 0.9835488 1.0012017 0.9562158 1.0324136 1.0505211
# [41] 0.9432645 1.0054619 0.9092097 0.8814829 0.9797460 0.9249077 0.9916127 0.9844348
#----------------#
lcpm <- cpm(data2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Normalised data",ylab="Log-cpm")

# Unsupervised clustering of samples
lcpm <- cpm(data2, log=TRUE)
par(mfrow=c(1,2))
col.stage <- stage
levels(col.stage) <-  brewer.pal(nlevels(col.stage), "Set1")
col.stage <- as.character(col.stage)
col.race <- race
levels(col.race) <-  brewer.pal(nlevels(col.race), "Set2")
col.race <- as.character(col.race)
plotMDS(lcpm, labels=stage, col=col.stage)
title(main="A. Stage Category")
plotMDS(lcpm, labels=race, col=col.race, dim=c(3,4))
title(main="B. Race")


# design <- model.matrix(~0+stage+race)
# design <- model.matrix(~0+stage)
design <- model.matrix(~0+classify)
colnames(design) <- gsub("classify", "", colnames(design))
colnames(design) <- gsub("race", "", colnames(design))
contr.matrix <- makeContrasts(
  # This will test for whether the log-expression of group Lean is equal to 
  # the average of groups ObeseI, ObeseII and ObeseIII.
  # Stage1.vs.Stage234 = Stage1-Stage234, 
  TNBC.vs.Discordant = TNBC-Discordant,
  levels = colnames(design))
contr.matrix


par(mfrow=c(1,2))
v <- voom(data2, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
dt <- decideTests(efit)
summary(dt)
#----------------#
# AA.vs.EA
# Down        16
# NotSig    19270
# Up          9
#----------------#

# tfit <- treat(vfit, lfc=1)
# summary(decideTests(tfit))


plotMD(efit, column=1, status=dt[,1], main=colnames(efit)[1], xlim=c(-3,13))

# vennDiagram(dt[,c(1,4)], circle.col=c("turquoise", "salmon"))
# 
# write.fit(efit, dt, file="results.txt")

# de.common <- which(dt[,1]!=0 & dt[,4]!=0)
de.treat <- which(dt[,1]!=0)
Ctrl.vs.Treat <- topTreat(efit, coef=1, n=Inf)
significant.gene <- row.names(Ctrl.vs.Treat[Ctrl.vs.Treat$adj.P.Val<0.05,])
write.table(Ctrl.vs.Treat[Ctrl.vs.Treat$adj.P.Val<0.05,], file = "215/215vs38.pval005.txt", sep="\t", quote = F, row.names=T)

up.gene <- row.names(Ctrl.vs.Treat[Ctrl.vs.Treat$adj.P.Val<0.05&Ctrl.vs.Treat$logFC>0,])
#----------------#
# [1] "TBC1D3C" "USP17L8" "SPINK2"  "FOLH1B"  "SPATA4"  "UGT2B17" "WIF1"    "PLXNB2"  "UNC45A"
#----------------#
down.gene <- row.names(Ctrl.vs.Treat[Ctrl.vs.Treat$adj.P.Val<0.05&Ctrl.vs.Treat$logFC<0,])
#----------------#
#  [1] "LOC284581" "USP17L4"   "MTO1"      "SMG1P5"    "TMEM184C"  "FAM107B"   "CLLU1-AS1" "FSHR"      "TMEM47"    "SPATC1L"   "NFYA"      "USP17L1"  
# [13] "ATF7IP"    "CFAP298"   "FUT8"      "GSTT2B"   
#----------------#

# Volcano draw
# volcanoplot(efit, coef=4,  main = "Lean.vs.ObeseIII", highlight=8, names = row.names(efit$p.value), hl.col="blue", xlim=c(-3,3), ylim=c(0,5))
# abline(v=c(-1,1), h=3, lty=3)
library(ggplot2)
logFC.cutoff <- 0.5
adj.P.Val.cutoff <- 0.05
Volcano.color <- as.factor(ifelse(Ctrl.vs.Treat$adj.P.Val<adj.P.Val.cutoff & 
                                    abs(Ctrl.vs.Treat$logFC)>=logFC.cutoff,
                                  ifelse(Ctrl.vs.Treat$logFC > logFC.cutoff,'Up','Down'),'Not'))

p <- ggplot(data=Ctrl.vs.Treat, aes(x=logFC, y=-log10(adj.P.Val), col=Volcano.color)) +
  geom_point(alpha=0.8, size = 1.8) +
  theme_bw(base_size = 15) +
  labs(x="log2 (fold change)",y="-log10 (adj.p-value)",title = "Molecular Confirmed vs Molecular Misclassified ") +
  geom_hline(yintercept = -log10(adj.P.Val.cutoff), lty=4,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff), lty=4,col="black",lwd=0.6) +
  scale_color_manual(name = "", values = c("red", "green", "gray"), limits = c("Up", "Down", "Not")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
plot(p)
Vcol.title <- "Molecular Confirmed vs Molecular Misclassified"

# Heatmap draw
library(gplots)
up.gene <- row.names(Ctrl.vs.Treat[Ctrl.vs.Treat$adj.P.Val<0.05&Ctrl.vs.Treat$logFC>2,])
down.gene <- row.names(Ctrl.vs.Treat[Ctrl.vs.Treat$adj.P.Val<0.05&Ctrl.vs.Treat$logFC<(-2),])
# Lean.vs.ObeseIII.topgenes <- row.names(Lean.vs.ObeseIII)[1:length(de.ObeseIII)]
Ctrl.vs.Treat.topgenes <- c(up.gene,down.gene)
i <- which(row.names(v) %in% Ctrl.vs.Treat.topgenes)
j <- c(which(v$targets$classify == "TNBC"), which(v$targets$classify == "Discordant"))
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,j], scale="row", Rowv = TRUE, Colv=FALSE, dendrogram="none",
          key = TRUE, keysize = 1.5, density.info="none",
          labRow=row.names(v)[i], labCol=classify[j], 
          col=mycol, trace="none", 
          margin=c(5,10), lhei=c(2,10), main = "Confirmed vs Misclassified ")
