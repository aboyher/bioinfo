library(ballgown)
library(dplyr)
library(genefilter)
library(FactoMineR)
library(ggplot2)

##change ballgowndir to whereever the ballgown folders are
ballgowndir = "/Users/aboyher/danforth/projects/cmd2/rnaseq10/tme7_v0.5j_0_chromo/ballgown"
setwd(ballgowndir)
#############################

sample_pattern = "ballgown"
data_dir = "."
pdata = read.csv("phenodata.csv",header =T, stringsAsFactors = F)
bg = ballgown(dataDir = data_dir, samplePattern = sample_pattern, pData = pdata)
#bg = ballgown(dataDir = data_dir, samplePattern = sample_pattern)
fpkm = texpr(bg, 'all')
colnames(fpkm)
fpkm = dplyr::select(fpkm,-t_id,-strand,-start,-end,-length,-num_exons)
fpkm = fpkm[, !grepl("^cov", names(fpkm))]
fpkm_table_name = "rnaseq10_tme7_v0.5j_0_chromo_fpkm.csv"
write.csv(fpkm, fpkm_table_name, row.names = F)
df = dplyr::select(fpkm, starts_with("FPKM"))
df = df[rowSums(df) > 0, ]
reps = as.data.frame(t(df))
pca <- PCA(reps, graph = F)
pca_df <- data.frame(
  "Samples" = rownames(reps),
  "PC1" = pca$ind$coord[, 1],
  "PC2" = pca$ind$coord[, 2]
)

####Code will vary by experiment
splitstring = "FPKM.ballgown_stringtie_tme7_v0.5j_0_chromo_rnaseq10_"
SPLIT = strsplit(row.names(pca_df), splitstring)
SPLIT = sapply(SPLIT, "[", 2)
SPLIT = sapply(strsplit(SPLIT, "_190625"), "[",1)
names = SPLIT
plant = gsub('.{1}_', '', names)
condition = sapply(strsplit(SPLIT, "_"), "[", 2)
pca_df$reps <- names
pca_df$samples <- condition
pca_df$plant <- plant
varexp <- signif(c(pca$eig[1, 2], pca$eig[2, 2]), 4)
plottitle = "RNAseq10 TME7_v0.5j_0_Chromo PCA"
plotsubtitle = ""
p = ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(data = aggregate(cbind(PC1, PC2) ~ samples, pca_df, mean),aes(color = samples),size = 5) +
  geom_point(aes(color = samples)) +
  # stat_ellipse(aes(fill=condition,color=condition),geom = "polygon",alpha=0.25, level = 0.8)+
  stat_ellipse(aes(fill=condition,color=condition),geom = "polygon",alpha=0.25)+
  # geom_polygon(aes(fill = pca_df$samples),show.legend = F,alpha = .3) +
  geom_vline(xintercept = 0,linetype = "dashed") +
  geom_hline(yintercept = 0,linetype = "dashed") +
  theme_light() +
  xlab(paste("PC1 (", varexp[1], "%)", sep = "")) +
  ylab(paste("PC2 (", varexp[2], "%)", sep = "")) +
  theme_light() +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=20, hjust=0)) +
  theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16)) +
  theme(legend.text= element_text(family = "Trebuchet MS", color="#666666", face="bold", size=13)) +
  theme(panel.border = element_rect(colour = "gray60",fill = NA,size = 1,linetype = 1)) +
  labs(title=plottitle, subtitle=plotsubtitle)+
  theme(aspect.ratio = 1)

p
#plot_filename = "rnaseq10_tme7_v0.5j_0_chromo_pca_cluster.png"
#ggsave(plot_filename, p)


bg_gexpr = gexpr(bg)
bg_gexpr
colnames(bg_gexpr)
rg = stattest(bg, feature="gene", meas="FPKM", covariate="condition", getFC = T)
rg$lf2c = log2(rg$fc)
rg_sig = filter(rg, pval < 0.05)
max(rg_sig$lf2c)
min(rg_sig$lf2c)
rg_sig = filter(rg_sig, abs(lf2c) > 2)
nrow(rg_sig)
write.csv(rg_sig, "rnaseq10_siggenes_pval.0.05_lf2c_2.csv", row.names = F)

rg_Manes = filter(rg, grepl("Manes", id))
max(rg_Manes$lf2c)
min(rg_Manes$lf2c)

rt = stattest(bg, feature = "transcript", meas = "FPKM", covariate = 'condition', getFC = T)
nrow(rt)
rt$lf2c = log2(rt$fc)
max(rt$lf2c)
min(rt$lf2c)
rt_sig = filter(rt, pval< 0.05)
nrow(rt_sig)
rt_sig = filter(rt_sig, abs(lf2c) > 2)
nrow(rt_sig)
bg_texpr = texpr(bg, 'all')
head(bg_texpr)
t_table = unique(dplyr::select(bg_texpr, t_id, t_name, gene_id, gene_name, chr, start, end))
rt_sig_merge = merge(rt_sig, t_table, by.x = c("id"), by.y=c("t_id"))
head(rt_sig_merge)
rt_sig_merge = select(rt_sig_merge, id, t_name, gene_id, gene_name, chr, start,end,pval,qval,fc,lf2c)
head(rt_sig_merge)
write.csv(rt_sig_merge, "rnaseq10_sigtranscripts_pval.0.1_lf2c_2.csv", row.names = F)

###Plotting

bg_filt = subset(bg, "rowVars(texpr(bg)) > 1", genomesubset = T)
bg_filt
results_transcripts = stattest(bg_filt, feature="transcript", covariate = "condition", getFC = T, meas = "FPKM")
dim(results_transcripts)
results_genes = stattest(bg_filt, feature = "gene", covariate = "condition", getFC = T, meas = "FPKM")
results_transcripts = data.frame(geneNames = geneNames(bg_filt), geneIDs = geneIDs(bg_filt), results_transcripts)
results_transcripts$mean = rowMeans(texpr(bg_filt))
p = ggplot(results_transcripts, aes(log2(mean), log2(fc), colour= qval<0.05)) + 
  scale_color_manual(values=c("#999999", "#FF0000")) + 
  geom_point() +
  geom_hline(yintercept=0)
p

nrow(rt)
rt_m = merge(rt, t_table, by.x=c("id"), by.y = c("t_id"))
rt_m = rt_m[!is.na(rt_m$pval),]
nrow(rt_m)
rt_m$Significance = "Not Significant"
rt_m[rt_m$pval < 0.05, "Significance"] = "P_Value < 0.05"
rt_m[abs(rt_m$lf2c) > 2, "Significance"] = "Absolute Log2FoldChange > 2"
rt_m[rt_m$pval < 0.05 & abs(rt_m$lf2c) >2, "Significance"] = "Absolute Log2FoldChange > 2 and Pvalue < 0.05"

p = ggplot(rt_m, aes(lf2c, -log10(pval))) + 
  geom_point(aes(colour=Significance)) + 
  scale_color_manual(values=c("#E43000","#53CB00","#999999", "#009DCB")) +
  theme_light() +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=20, hjust=0)) +
  theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16)) +
  theme(legend.text= element_text(family = "Trebuchet MS", color="#666666", face="bold", size=13)) +
  theme(panel.border = element_rect(colour = "gray60",fill = NA,size = 1,linetype = 1)) +
  labs(title="RNAseq10 Volcano Plot", subtitle="TME7_V0.5.3a") +
  labs(x = "Log2FoldChange") +
  geom_vline(xintercept = 2,linetype = "dashed") +
  geom_vline(xintercept = -2,linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed") +
  theme(aspect.ratio = 1)
p  
ggsave("rnaseq10_tme7_v0.5.3a_volcano.png",p)

rt_m = merge(rt, t_table, by.x=c("id"), by.y = c("t_id"))
rt_m = rt_m[!is.na(rt_m$qval),]
rt_m$Significance = "Not Significant"
rt_m[rt_m$qval < 0.05, "Significance"] = "Q_Value < 0.05"
rt_m[abs(rt_m$lf2c) > 2, "Significance"] = "Absolute Log2FoldChange > 2"
rt_m[rt_m$qval < 0.05 & abs(rt_m$lf2c) >2, "Significance"] = "Absolute Log2FoldChange > 2 and Qvalue < 0.05"

p = ggplot(rt_m, aes(lf2c, -log10(qval))) + 
  geom_point(aes(colour=Significance)) + 
  scale_color_manual(values=c("#E43000","#53CB00","#999999", "#009DCB")) +
  theme_light() +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=20, hjust=0)) +
  theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16)) +
  theme(legend.text= element_text(family = "Trebuchet MS", color="#666666", face="bold", size=13)) +
  theme(panel.border = element_rect(colour = "gray60",fill = NA,size = 1,linetype = 1)) +
  labs(title="RNAseq10 Volcano Plot", subtitle="TME7_v0.5.3a") +
  labs(x = "Log2FoldChange") +
  geom_vline(xintercept = 2,linetype = "dashed") +
  geom_vline(xintercept = -2,linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed") +
  theme(aspect.ratio = 1)
p  
ggsave("rnaseq10_tme7_v0.5.3a_volcano_adjusted_q.png",p)


########## Manes.12G063800 (Super-Scaffold_26|quiver|pilon:7404779-7407442) and Manes.12G063900 (Super-Scaffold_26|quiver|pilon:7408027-7411111)
rt = stattest(bg, feature = "transcript", meas = "FPKM", covariate = 'condition', getFC = T)
bg_texpr = texpr(bg, 'all')
t_table = unique(dplyr::select(bg_texpr, t_id, t_name, gene_id, gene_name, chr, start, end))
rt = merge(rt, t_table, by.x = c("id"), by.y=c("t_id"))
head(rt)
rt_filt = filter(rt, grepl("Super-Scaffold_26", chr)
                 
                 