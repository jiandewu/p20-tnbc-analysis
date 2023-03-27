library(limma)
library(edgeR)


options(stringsAsFactors = FALSE)

setwd(Sys.getenv("P20_DATA_DIR", unset = "."))
rawdata <- read.table("258samples.txt", header = TRUE, sep="\t")

# delete some samples
# rawdata <- rawdata[, !(colnames(rawdata) %in% c("CONTROLLED_SAMPLE_ID"))]
# rawdata[3, (colnames(rawdata) %in% c("CONTROLLED_SAMPLE_ID"))] =30

# Find error symbols
##sort(rawdata[,1])[c(1:27)]
#----------------#
# [1] "1-Dec"  "1-Mar"  "1-Mar"  "1-Sep"  "10-Mar" "10-Sep" "11-Mar" "11-Sep" "12-Sep" "14-Sep" "2-Mar"  "2-Mar" 
# [13] "2-Sep"  "3-Mar"  "3-Sep"  "4-Mar"  "4-Sep"  "5-Mar"  "5-Sep"  "6-Mar"  "6-Sep"  "7-Mar"  "7-Sep"  "8-Mar" 
# [25] "8-Sep"  "9-Mar"  "9-Sep" 
#----------------#

# Find duplicate rows
##n_occur <- data.frame(table(rawdata[,1]))
##sort(rawdata[rawdata[,1] %in% n_occur$Var1[n_occur$Freq > 1],1])
#----------------#
# [1] "1-Mar"        "1-Mar"        "2-Mar"        "2-Mar"        "LOC110117498" "LOC110117498" "SUGT1P4-STRA"
# [8] "SUGT1P4-STRA"
#----------------#

# delete bilate samples
rawdata <- rawdata[,c(-245,-244,-149,-148)]


# Rename duplicate rows
##row.names(rawdata) <- ave(as.character(rawdata[,1]), rawdata[,1], FUN=function(x) if (length(x)>1) paste0(x[1], '(', seq_along(x), ')') else x[1])

row.names(rawdata) <- rawdata[,1]

# delete AA
##rawdata <- rawdata[,rawdata[1,]==1]
# delete EA
rawdata <- rawdata[,rawdata[1,]>1]
rawdata <- rawdata[,-1]

# data1 <- DGEList(counts = rawdata[c(-1:-3),-1])
data1 <- DGEList(counts = rawdata[c(-1:-3),])

obeselevel <- rawdata[2,]
obeselevel[obeselevel==1] = "Lean"
obeselevel[obeselevel==2] = "ObeseI"
obeselevel[obeselevel==3] = "ObeseII"
obeselevel[obeselevel==4] = "ObeseIII"
obese <- as.factor(obeselevel)
racelevel <- rawdata[1,]
racelevel[racelevel>1] = "AA"
racelevel[racelevel==1] = "EA"
race <- as.factor(racelevel)
agelevel <- rawdata[3,]
agelevel[agelevel>=65] = "Old"
agelevel[agelevel>=45 & agelevel<65] = "Middle"
agelevel[agelevel<45] = "Young"
age <- as.factor(agelevel)
data1$samples$obese <- obese
data1$samples$race <- race
data1$samples$age <- age

table(data1$samples$obese,data1$samples$race)
#----------------#
# AA EA
# Lean     27 50
# ObeseI   42 39
# ObeseII  31 20
# ObeseIII 29 16
#----------------#

# Transformations from the raw-scale
seqdepth_mean <- mean(data1$samples$lib.size) * 1e-6
seqdepth_median <- median(data1$samples$lib.size) * 1e-6

summary(data1$samples$lib.size * 1e-6)
#--All Data------#
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9896 15.9975 19.1382 19.3924 23.1632 36.3612 
#----------------#
#--EA Data-------#
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.32   15.36   18.38   19.03   22.80   35.98 
#--AA Data-------#
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9896 16.5911 19.5577 19.7418 23.4279 36.3612 

nsamples <- ncol(data1)

# Check how many rows are all 0s
table(rowSums(data1$counts==0)==nsamples)
#--All Data------#
# FALSE  TRUE 
# 25511  2767 
#----------------#
#--EA Data-------#
# FALSE  TRUE 
# 25059  3219 
#----------------#
# FALSE  TRUE 
# 25164  3114 

#########################################################################################
# for filter process.
group <- obese
#########################################################################################

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
data1.filted <- data1[keep.exprs,, keep.lib.sizes=FALSE]
dim(data1.filted)
# [1] 18089   254
# [1] 18183   125
# [1] 17904   129

lcpm <- cpm(data1.filted, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(data1), text.col=col, bty="n")
#########################################################################################

# Normalising gene expression distributions
par(mfrow=c(1,2))
lcpm <- cpm(data1.filted, log=TRUE, prior.count=2)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="Log-cpm")

data1.filted.norm <- calcNormFactors(data1.filted, method = "TMM")
# data2$samples$norm.factors
#----------------#
#  [1] 1.1005944 1.0213063 0.9484682 1.0938885 0.9135321 1.1018828 1.0813619 0.9746213 0.9708049 1.0128725
# [11] 0.9602972 1.0116042 1.0162326 1.0180576 1.0575765 0.9646913 1.0454093 1.0034687 0.9942733 1.0220867
# [21] 1.0738821 1.0687922 1.0386444 1.0766185 1.0205136 0.9873789 0.9493396 0.9206837 1.0122073 1.0469196
# [31] 1.0697226 1.0397318 0.9758840 0.9220561 0.9134998 0.9835488 1.0012017 0.9562158 1.0324136 1.0505211
# [41] 0.9432645 1.0054619 0.9092097 0.8814829 0.9797460 0.9249077 0.9916127 0.9844348
#----------------#

lcpm <- cpm(data1.filted.norm, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Normalised data",ylab="Log-cpm")

# Unsupervised clustering of samples
lcpm <- cpm(data1.filted.norm, log=TRUE)
par(mfrow=c(1,2))
col.obese <- obese
levels(col.obese) <-  brewer.pal(nlevels(col.obese), "Set1")
col.obese <- as.character(col.obese)
col.race <- race
levels(col.race) <-  brewer.pal(nlevels(col.race), "Set2")
col.race <- as.character(col.race)
plotMDS(lcpm, labels=obese, col=col.obese)
title(main="A. Obesity Category")
plotMDS(lcpm, labels=race, col=col.race, dim=c(3,4))
title(main="B. Race")


#########################################################################################
# for filter process.
group <- obese
#########################################################################################
design <- model.matrix(~0+obese+age)
# design <- model.matrix(~0+group)
colnames(design) <- gsub("obese", "", colnames(design))
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

data2 <- data1.filted.norm
par(mfrow=c(1,2))
v <- voom(data2, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
dt <- decideTests(efit)
summary(dt)
#----------------#
# Lean.vs.Obese Lean.vs.ObeseI Lean.vs.ObeseII Lean.vs.ObeseIII
# Down               0              0               0                3
# NotSig         18089          18089           18089            18083
# Up                 0              0               0                3
#----------------#
# Lean.vs.Obese Lean.vs.ObeseI Lean.vs.ObeseII Lean.vs.ObeseIII
# Down               0              0               0                2
# NotSig         18183          18183           18183            18180
# Up                 0              0               0                1
#----------------#
# Lean.vs.Obese Lean.vs.ObeseI Lean.vs.ObeseII Lean.vs.ObeseIII
# Down               1              0               0                0
# NotSig         17903          17904           17904            17904
# Up                 0              0               0                0

# tfit <- treat(vfit, lfc=1)
# summary(decideTests(tfit))


plotMD(efit, column=1, status=dt[,1], main=colnames(efit)[1], xlim=c(-3,13))

# vennDiagram(dt[,c(1,4)], circle.col=c("turquoise", "salmon"))
# 
# write.fit(efit, dt, file="results.txt")

# de.common <- which(dt[,1]!=0 & dt[,4]!=0)
de.Obese <- which(dt[,1]!=0)
Lean.vs.Obese <- topTreat(efit, coef=1, n=Inf) # coef = summary column
up.gene <- row.names(Lean.vs.Obese[Lean.vs.Obese$adj.P.Val<0.05&Lean.vs.Obese$logFC>0,])
#----------------#
# [1] "MUSK"   "XIST"   "GPCPD1"
#----------------#
# [1] "ITPR1"

down.gene <- row.names(Lean.vs.Obese[Lean.vs.Obese$adj.P.Val<0.05&Lean.vs.Obese$logFC<0,])
#----------------#
# [1] "STK32C"        "CNPY4"         "CCZ1P-OR7E38P"
#----------------#
# [1] "PXDNL"  "STK32C"
#----------------#
# [1] "INO80B"


# Volcano draw
# volcanoplot(efit, coef=4,  main = "Lean.vs.ObeseIII", highlight=8, names = row.names(efit$p.value), hl.col="blue", xlim=c(-3,3), ylim=c(0,5))
# abline(v=c(-1,1), h=3, lty=3)
library(ggplot2)
logFC.cutoff <- 0.1
adj.P.Val.cutoff <- 0.05
Volcano.color <- as.factor(ifelse(Lean.vs.Obese$adj.P.Val<adj.P.Val.cutoff & 
                                    abs(Lean.vs.Obese$logFC)>=logFC.cutoff,
                                  ifelse(Lean.vs.Obese$logFC > logFC.cutoff,'Up','Down'),'Not'))

p <- ggplot(data=Lean.vs.Obese, aes(x=logFC, y=-log10(adj.P.Val), col=Volcano.color)) +
  geom_point(alpha=0.8, size = 1.8) +
  theme_bw(base_size = 15) +
  labs(x="log2 (fold change)",y="-log10 (adj.p-value)",title = "Lean.vs.ObeseIII") +
  geom_hline(yintercept = -log10(adj.P.Val.cutoff), lty=4,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff), lty=4,col="black",lwd=0.6) +
  scale_color_manual(name = "", values = c("red", "green", "gray"), limits = c("Up", "Down", "Not")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
plot(p)

# Heatmap draw
library(gplots)
# Lean.vs.ObeseIII.topgenes <- row.names(Lean.vs.ObeseIII)[1:length(de.ObeseIII)]
Lean.vs.Obese.topgenes <- c(up.gene,down.gene)
i <- which(row.names(v) %in% Lean.vs.Obese.topgenes)
j <- c(which(v$targets$group == "Lean"), which(v$targets$group == "ObeseI"), 
       which(v$targets$group == "ObeseII"), which(v$targets$group == "ObeseIII"))
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,j], scale="row", Rowv = TRUE, Colv=FALSE, dendrogram="none",
          key = TRUE, keysize = 1.5, density.info="none",
          labRow=row.names(v)[i], labCol=group[j], 
          col=mycol, trace="none", 
          margin=c(16,12), lhei=c(2,10), main = "Lean vs Obese Heatmap")
