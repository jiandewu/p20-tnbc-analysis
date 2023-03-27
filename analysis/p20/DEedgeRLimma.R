library(limma)
library(edgeR)


options(stringsAsFactors = FALSE)

setwd(Sys.getenv("P20_DATA_DIR", unset = "."))
rawdata <- read.table("P20_Annotated_Export_10032019.txt", header = TRUE, sep="\t")
## Individual-level data example omitted from the public repository.
setwd(Sys.getenv("P20_DATA_DIR", unset = "."))

# Find error symbols
sort(rawdata[,1])[c(1:27)]
#----------------#
# [1] "1-Dec"  "1-Mar"  "1-Mar"  "1-Sep"  "10-Mar" "10-Sep" "11-Mar" "11-Sep" "12-Sep" "14-Sep" "2-Mar"  "2-Mar" 
# [13] "2-Sep"  "3-Mar"  "3-Sep"  "4-Mar"  "4-Sep"  "5-Mar"  "5-Sep"  "6-Mar"  "6-Sep"  "7-Mar"  "7-Sep"  "8-Mar" 
# [25] "8-Sep"  "9-Mar"  "9-Sep" 
#----------------#

# Find duplicate rows
n_occur <- data.frame(table(rawdata[,1]))
sort(rawdata[rawdata[,1] %in% n_occur$Var1[n_occur$Freq > 1],1])
#----------------#
# [1] "1-Mar"        "1-Mar"        "2-Mar"        "2-Mar"        "LOC110117498" "LOC110117498" "SUGT1P4-STRA"
# [8] "SUGT1P4-STRA"
#----------------#

# Rename duplicate rows
row.names(rawdata) <- ave(as.character(rawdata[,1]), rawdata[,1], FUN=function(x) if (length(x)>1) paste0(x[1], '(', seq_along(x), ')') else x[1])

data1 <- DGEList(counts = rawdata[c(-1:-3,-nrow(rawdata)+1,-nrow(rawdata)),-1])

grouplevel <- rawdata[2,-1]
grouplevel[grouplevel==1] = "Lean"
grouplevel[grouplevel==2] = "ObeseI"
grouplevel[grouplevel==3] = "ObeseII"
grouplevel[grouplevel==4] = "ObeseIII"
group <- as.factor(grouplevel)
racelevel <- rawdata[1,-1]
racelevel[racelevel>1] = "Black"
racelevel[racelevel==1] = "White"
race <- as.factor(racelevel)
pool <- as.factor(rawdata[3,-1])
data1$samples$group <- group
data1$samples$race <- race
data1$samples$pool <- pool

table(data1$samples$group,data1$samples$race)
#----------------#
#          Black White
# Lean         4    10
# ObeseI       8     9
# ObeseII      8     2
# ObeseIII     2     5
#----------------#

# Transformations from the raw-scale
seqdepth_mean <- mean(data1$samples$lib.size) * 1e-6
seqdepth_median <- median(data1$samples$lib.size) * 1e-6

summary(data1$samples$lib.size * 1e-6)
#----------------#
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.069   3.008   6.251   6.527   9.616  14.074 
#----------------#
nsamples <- ncol(data1)

table(rowSums(data1$counts==0)==nsamples)
#----------------#
# FALSE  TRUE 
# 23382  4589 
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
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(data1), text.col=col, bty="n")

# Removing genes that are lowly expressed
keep.exprs <- filterByExpr(data1, group=group)
data2 <- data1[keep.exprs,, keep.lib.sizes=FALSE]
dim(data2)
## [1] 15937    48

lcpm <- cpm(data2, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
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
data2$samples$norm.factors
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
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.race <- race
levels(col.race) <-  brewer.pal(nlevels(col.race), "Set2")
col.race <- as.character(col.race)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Obesity Category")
plotMDS(lcpm, labels=race, col=col.race, dim=c(3,4))
title(main="B. Race")


# design <- model.matrix(~0+group+race)
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
colnames(design) <- gsub("race", "", colnames(design))
contr.matrix <- makeContrasts(
  # This will test for whether the log-expression of group Lean is equal to 
  # the average of groups ObeseI, ObeseII and ObeseIII.
  Lean.vs.Obese = Lean-(ObeseI+ObeseII+ObeseIII)/3, 
  Lean.vs.ObeseI = Lean-ObeseI, 
  Lean.vs.ObeseII = Lean-ObeseII, 
  Lean.vs.ObeseIII = Lean-ObeseIII, 
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
#        Lean.vs.Obese Lean.vs.ObeseI Lean.vs.ObeseII Lean.vs.ObeseIII
# Down               0              0               0               22
# NotSig         15937          15937           15937            15899
# Up                 0              0               0               16
#----------------#

# tfit <- treat(vfit, lfc=1)
# summary(decideTests(tfit))


plotMD(efit, column=4, status=dt[,4], main=colnames(efit)[4], xlim=c(-3,13))

# vennDiagram(dt[,c(1,4)], circle.col=c("turquoise", "salmon"))
# 
# write.fit(efit, dt, file="results.txt")

# de.common <- which(dt[,1]!=0 & dt[,4]!=0)
de.ObeseIII <- which(dt[,4]!=0)
Lean.vs.ObeseIII <- topTreat(efit, coef=4, n=Inf)
up.gene <- row.names(Lean.vs.ObeseIII[Lean.vs.ObeseIII$adj.P.Val<0.05&Lean.vs.ObeseIII$logFC>0,])
#----------------#
#  [1] "ZNF841"    "EPM2AIP1"  "KIF27"     "ZNF586"    "ZNF808"    "LINC01061" "ZMYM2"     "ZNF224"    "ZNF587"   
# [10] "GUSBP11"   "ZNF582"    "CCNL2"     "DND1"      "ZNF432"    "ZNF700"    "ZNF70"    
#----------------#
down.gene <- row.names(Lean.vs.ObeseIII[Lean.vs.ObeseIII$adj.P.Val<0.05&Lean.vs.ObeseIII$logFC<0,])
#----------------#
#  [1] "PCMT1"     "BROX"      "CDC40"     "GNAI3"     "SNX7"      "GTF2H5"    "HEBP2"     "VTA1"      "ELAVL2"   
# [10] "SNAPIN"    "DCLRE1B"   "STXBP3"    "KATNA1"    "RSU1"      "ZNHIT3"    "PRUNE1"    "NEK7"      "COG2"     
# [19] "TATDN3"    "AIDA"      "CAPZA1"    "LRRC37A6P"
#----------------#

# Volcano draw
# volcanoplot(efit, coef=4,  main = "Lean.vs.ObeseIII", highlight=8, names = row.names(efit$p.value), hl.col="blue", xlim=c(-3,3), ylim=c(0,5))
# abline(v=c(-1,1), h=3, lty=3)
library(ggplot2)
logFC.cutoff <- 1
adj.P.Val.cutoff <- 0.05
Volcano.color <- as.factor(ifelse(Lean.vs.ObeseIII$adj.P.Val<adj.P.Val.cutoff & 
                                    abs(Lean.vs.ObeseIII$logFC)>=logFC.cutoff,
                                  ifelse(Lean.vs.ObeseIII$logFC > logFC.cutoff,'Up','Down'),'Not'))

p <- ggplot(data=Lean.vs.ObeseIII, aes(x=logFC, y=-log10(adj.P.Val), col=Volcano.color)) +
  geom_point(alpha=0.8, size = 0.8) +
  theme_bw(base_size = 15) +
  labs(x="log2 (fold change)",y="-log10 (adj.p-value)",title = "Lean.vs.ObeseIII") +
  geom_hline(yintercept = -log10(adj.P.Val.cutoff), lty=4,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff), lty=4,col="black",lwd=0.6) +
  scale_color_manual(name = "", values = c("green", "red", "gray"), limits = c("Up", "Down", "Not")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
plot(p)

# Heatmap draw
library(gplots)
# Lean.vs.ObeseIII.topgenes <- row.names(Lean.vs.ObeseIII)[1:length(de.ObeseIII)]
Lean.vs.ObeseIII.topgenes <- c(up.gene,down.gene)
i <- which(row.names(v) %in% Lean.vs.ObeseIII.topgenes)
j <- c(which(v$targets$group == "Lean"), which(v$targets$group == "ObeseIII"))
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,j], scale="row", Rowv = TRUE, Colv=FALSE, dendrogram="none",
          key = TRUE, keysize = 1.5, density.info="none",
          labRow=row.names(v)[i], labCol=group[j], 
          col=mycol, trace="none", 
          margin=c(8,6), lhei=c(2,10), main = "Lean vs ObeseIII Heatmap")
