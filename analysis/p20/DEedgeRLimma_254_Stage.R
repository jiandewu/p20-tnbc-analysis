library(limma)
library(edgeR)


options(stringsAsFactors = FALSE)

setwd(Sys.getenv("P20_DATA_DIR", unset = "."))
rawdata <- read.table("258samples_v3.txt", header = TRUE, sep="\t")
# rawdata <- rawdata[, !(colnames(rawdata) %in% c("CONTROLLED_SAMPLE_ID"))]
# delete bilate samples
# rawdata <- rawdata[,c(-245,-244,-149,-148)]
drops <- if (nzchar(Sys.getenv("P20_BILATERAL_EXCLUSIONS_FILE"))) readLines(Sys.getenv("P20_BILATERAL_EXCLUSIONS_FILE")) else character(0)
rawdata <- rawdata[, !(names(rawdata) %in% drops)] # all samples
row.names(rawdata) <- rawdata[,1]
rawdata <- rawdata[,-1]

# rawdata <- rawdata[, rawdata["AA",]==1] # AA samples
# rawdata <- rawdata[, rawdata["AA",]==0] # EA samples
rawdata <- rawdata[, rawdata["LATE",]==1] # Late Stage samples


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
stagelevel[stagelevel==1] = "Stage2-4"
stagelevel[stagelevel==0] = "Stage1"
stage <- as.factor(stagelevel)
data1$samples$obese <- obese
data1$samples$race <- race
data1$samples$stage <- stage

table(data1$samples$stage,data1$samples$race)
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
keep.exprs <- filterByExpr(data1, group=race)
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
design <- model.matrix(~0+race)
colnames(design) <- gsub("stage", "", colnames(design))
colnames(design) <- gsub("race", "", colnames(design))
contr.matrix <- makeContrasts(
  # This will test for whether the log-expression of group Lean is equal to 
  # the average of groups ObeseI, ObeseII and ObeseIII.
  # Stage1.vs.Stage234 = Stage1-Stage234, 
  AA.vs.EA = AA-EA,
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
# Down        122
# NotSig    16675
# Up          115
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
write.table(Ctrl.vs.Treat[Ctrl.vs.Treat$adj.P.Val<0.05,], file = "Stage1.vs.Stage234/Stage234.AA_EA.pval005.txt", sep="\t", quote = F, row.names=T)

up.gene <- row.names(Ctrl.vs.Treat[Ctrl.vs.Treat$adj.P.Val<0.05&Ctrl.vs.Treat$logFC>0,])
#----------------#
# [1] "RNFT1-DT"         "TBC1D3C"          "NUTM2F"           "NPIPB15"          "SORD2P"           "C17orf97"         "NUTM2E"           "CCL3L3"
# [9] "TBC1D3P1-DHX40P1" "FOLH1B"           "SLC35E2B"         "CCL3L1"           "RASA4B"           "PLA2G4C"          "RPS26"            "L1TD1"
# [17] "ANKRD20A9P"       "CBWD3"            "SCAMP5"           "ZNF285"           "TAS2R43"          "NOTCH2NLC"        "ASAH2B"           "TSPAN10"
# [25] "DDX11L1"          "CUTALP"           "RPS26P11"         "WIF1"             "C4BPA"            "FBXO39"           "AMIGO1"           "HPR"
# [33] "GPR89B"           "SNORA38"          "RPL14"            "POLR2J3"          "JSRP1"            "ANKRD36BP2"       "HERC2P3"          "TMEM232"
# [41] "IZUMO1"           "TNFRSF17"         "FZD8"             "FOXS1"            "LRRIQ3"           "SHLD2P1"          "SCARNA1"          "TRIM61"
# [49] "SCARNA23"         "LOC441455"        "TUBA3E"           "STAP1"            "SNORA24"          "PCDHA11"          "IGHV1OR15-1"      "LOC102724760"
# [57] "BOLA2-SMG1P6"     "SNORA26"          "TUBB8"            "HERPUD1"          "SNORA63"          "GJB5"             "SLC22A16"         "XKR8"
# [65] "UPK3BL2"          "DPYSL4"           "GSTM1"            "OSER1"            "CRMP1"            "SOX30"            "SNAP91"           "LINC01530"
# [73] "CHST5"            "PRH1-PRR4"        "LOC100506388"     "KCNQ5"            "FAM155B"          "LOC102724624"     "IGLL5"            "POLR2J2"
# [81] "NCMAP"            "RPL10A"           "POP5"             "SNORA79"          "GSTM3"            "TMEM255A"         "DHRS4L2"          "GSTO2"
# [89] "ZNF561-AS1"       "PLRG1"            "GSTM2P1"          "SPDYE2"           "SPDYE2B"          "CDH4"             "UPK3BL1"          "ERO1B"
# [97] "ATP8A1"           "MTMR8"            "MFSD8"            "TMA16"            "CRISPLD1"         "CDKL4"            "FAM155A"          "SNORD116-17"
# [105] "SNORD116-19"      "TOX"              "ZNF883"           "SEC14L2"          "CDKN2B"           "TTLL2"            "GGT2"             "MAP9"
# [113] "GUSBP2"           "CCDC74B"          "WDR45"
#----------------#
down.gene <- row.names(Ctrl.vs.Treat[Ctrl.vs.Treat$adj.P.Val<0.05&Ctrl.vs.Treat$logFC<0,])
#----------------#
# [1] "LOC284581"     "SNORD3B-1"     "SNORD3B-2"     "HERC2P2"       "LRRC37A2"      "NPIPA3"        "EP400P1"       "ADAL"          "NOTCH2NLR"     "TBC1D3I"
# [11] "TBC1D3H"       "SPATC1L"       "HTATSF1P2"     "PKDREJ"        "NPIPA2"        "ARL17B"        "SPDYE17"       "SMG1P5"        "USP17L1"       "PM20D1"
# [21] "LRRC37A3"      "KANSL1"        "GOLGA6A"       "USP17L4"       "TBC1D3L"       "RGPD2"         "HERC2P9"       "LOC729737"     "ERV3-1-ZNF117" "MAPK8IP1P2"
# [31] "ARL17A"        "CCDC144A"      "USP32P1"       "HERC2"         "GTF2IP4"       "ENGASE"        "FAM157B"       "NPIPA5"        "ASAH2"         "RGPD1"
# [41] "FAM106C"       "SNORA25"       "SMG1P6"        "ZNF514"        "ERV3-1"        "NPIPB3"        "NBPF26"        "PPIL3"         "NWD2"          "SPDYE10P"
# [51] "USP6"          "WASHC1"        "LOC101926935"  "LOC100132062"  "LOC100132287"  "GUCY1B2"       "ZNF117"        "ZSCAN29"       "SMG1P2"        "JDP2"
# [61] "LINC01001"     "SPATA20"       "CDH26"         "LILRB3"        "MARCO"         "AIFM3"         "HRNR"          "ZNF600"        "ZNF587"        "CDRT1"
# [71] "RPL23AP7"      "DHX58"         "TP53BP1"       "HEATR3"        "CCDC122"       "TBC1D3D"       "PPFIA1"        "TRIM66"        "POLR1A"        "CLEC4G"
# [81] "MZT2A"         "TMEM236"       "DLG5"          "FBXL19-AS1"    "NPIPB12"       "BMS1P22"       "CEP152"        "PLIN5"         "GOLGA6L9"      "CD209"
# [91] "COG4"          "CCDC146"       "ATAD3C"        "DYNC2I1"       "GTF2F1"        "OR1L4"         "LRMDA"         "STYXL1"        "MRC1"          "MEIS3P1"
# [101] "CCL13"         "C4orf19"       "MTIF3"         "POM121"        "ATG7"          "LOC646652"     "ACCS"          "CARD9"         "APIP"          "DND1"
# [111] "RAD51L3-RFFL"  "LOC649352"     "PRG4"          "LOC101929141"  "NEURL2"        "AGAP13P"       "HK3"           "STAG3L2"       "GPX2"          "CCDC144B"
# [121] "NPEPL1"        "TAF3"
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
  labs(x="log2 (fold change)",y="-log10 (adj.p-value)",title = "AA vs EA") +
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
up.gene <- row.names(Ctrl.vs.Treat[Ctrl.vs.Treat$adj.P.Val<0.05&Ctrl.vs.Treat$logFC>1,])
down.gene <- row.names(Ctrl.vs.Treat[Ctrl.vs.Treat$adj.P.Val<0.05&Ctrl.vs.Treat$logFC<(-0.01),])
# Lean.vs.ObeseIII.topgenes <- row.names(Lean.vs.ObeseIII)[1:length(de.ObeseIII)]
Ctrl.vs.Treat.topgenes <- c(up.gene,down.gene)
i <- which(row.names(v) %in% Ctrl.vs.Treat.topgenes)
j <- c(which(v$targets$stage == "Stage1"), which(v$targets$stage == "Stage234"))
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,j], scale="row", Rowv = TRUE, Colv=FALSE, dendrogram="none",
          key = TRUE, keysize = 1.5, density.info="none",
          labRow=row.names(v)[i], labCol=stage[j], 
          col=mycol, trace="none", 
          margin=c(5,10), lhei=c(2,10), main = "Stage1 vs Stage234")
