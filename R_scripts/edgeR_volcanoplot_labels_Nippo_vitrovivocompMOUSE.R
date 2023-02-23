library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
library(data.table)
library(edgeR)



#################### Remove 3rd rep ####################
counts = read.delim('./data/working_matrices/Nippo_Allsamples_02-13-23_counts_matrix_rm_col_filabv1TPM.txt',header=TRUE,row.names="gene_id")




## Make data frame into matrix
genes = counts$gene_id
counts$gene_id = NULL
counts = as.matrix(counts)
rownames(counts) = genes
head(counts)

counts_0 = counts[rowSums(counts[,-1])>0,]
dim(counts_0) # 37890 by 24


########## Load in metadata ####################

samples_1<-read.delim('Metadata_samples_2-20-23.txt',header=TRUE)
(samples_1) #### look at samples file to make there is no NA or empty rows.
samples_1$sample_name2 = gsub("-","\\.",samples_1$sample_name)


# sanity check that all samples in metadata are in colnames of data
table(samples_1$sample_name2 %in% colnames(counts_0))


genelist = genes

#################### edgeR ####################
# edgeR uses COUNTS

# Set group
group <- samples_1$stage8_sampletype


y <- DGEList(counts=counts_0, group=group)


keep <- filterByExpr(y) 
y <- y[keep, , keep.lib.sizes=FALSE]

y <- estimateDisp(y) # 15437 genes left

# if the pair is c("A","B") then the comparison is B - A, 
# so genes with positive log-fold change are up-regulated in group B compared with group A 
# Calculate DE between B = FI and A = ctrl


### VT comparison all to day 0 Mouse - MUST CHANGE TITLES ###

##CHANGE COLORS!! Change orchid1 to plum1 royalblue to orchid1

###"plum1", "grey","orchid1"


edger_C8 <- exactTest(y, pair=c("VT_0d_mFem","VT_1d_mFem"))$table

edger_C9 <- exactTest(y, pair=c("VT_0d_mMale","VT_1d_mMale"))$table

edger_C10 <- exactTest(y, pair=c("VT_0d_mFem","VT_2d_mFem"))$table

edger_C11 <- exactTest(y, pair=c("VT_0d_mMale","VT_2d_mMale"))$table

edger_C12 <- exactTest(y, pair=c("VT_0d_mFem","VT_3d_mFem"))$table

edger_C13 <- exactTest(y, pair=c("VT_0d_mMale","VT_3d_mMale"))$table

edger_C14 <- exactTest(y, pair=c("VT_0d_mFem","VT_5d_mFem"))$table

edger_C15  <- exactTest(y, pair=c("VT_0d_mMale","VT_5d_mMale"))$table

edger_C16 <- exactTest(y, pair=c("VT_0d_mFem","VT_7d_mFem"))$table

edger_C17 <- exactTest(y, pair=c("VT_0d_mMale","VT_7d_mMale"))$table

#################### Volcano plots ####################
library(ggrepel)
library(ggpubr)

## C8
# add a column of NAs
edger_C8$diffexpressed <- "NS"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
edger_C8$diffexpressed[edger_C8$logFC > 1 & edger_C8$PValue < 0.05] <- "Up"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
edger_C8$diffexpressed[edger_C8$logFC < -1 & edger_C8$PValue < 0.05] <- "Down"

table(edger_C8$diffexpressed)

edger_C8$delabel <- NA
edger_C8$delabel[which(rownames(edger_C8) %in% genelist)] <- rownames(edger_C8)[which(rownames(edger_C8) %in% genelist)]

volcano_C8 = ggplot(data=edger_C8, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(box.padding = 6,size = 5,color="black") +
  scale_color_manual(values=c("plum1", "grey","orchid1")) +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  xlim(-15,15) +ylim(0,100) +  ggtitle("Mouse-derived day 0 Vs. day 1 Females") +
  geom_text(size=6, x=-10, y=90, label=paste(table(edger_C8$diffexpressed)[1]), color="black") +
  geom_text(size=6,x=10, y=90, label=paste(table(edger_C8$diffexpressed)[3]), color="black") +
  theme(plot.title = element_text(size=18),legend.position = "none",axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16), axis.title.x = element_text(size=16),axis.text.x = element_text(size=16))


## C9
# add a column of NAs
edger_C9$diffexpressed <- "NS"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
edger_C9$diffexpressed[edger_C9$logFC > 1 & edger_C9$PValue < 0.05] <- "Up"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
edger_C9$diffexpressed[edger_C9$logFC < -1 & edger_C9$PValue < 0.05] <- "Down"

table(edger_C9$diffexpressed)

edger_C9$delabel <- NA
edger_C9$delabel[which(rownames(edger_C9) %in% genelist)] <- rownames(edger_C9)[which(rownames(edger_C9) %in% genelist)]

volcano_C9 = ggplot(data=edger_C9, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(box.padding = 6,size = 5,color="black") +
  scale_color_manual(values=c("plum1", "grey","orchid1")) +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  xlim(-15,15) +ylim(0,100) +  ggtitle("Mouse-derived day 0 Vs. day 1 Males") +
  geom_text(size=6, x=-10, y=90, label=paste(table(edger_C9$diffexpressed)[1]), color="black") +
  geom_text(size=6,x=10, y=90, label=paste(table(edger_C9$diffexpressed)[3]), color="black") +
  theme(plot.title = element_text(size=18),legend.position = "none",axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16), axis.title.x = element_text(size=16),axis.text.x = element_text(size=16))

## C10
# add a column of NAs
edger_C10$diffexpressed <- "NS"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
edger_C10$diffexpressed[edger_C10$logFC > 1 & edger_C10$PValue < 0.05] <- "Up"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
edger_C10$diffexpressed[edger_C10$logFC < -1 & edger_C10$PValue < 0.05] <- "Down"

table(edger_C10$diffexpressed)

edger_C10$delabel <- NA
edger_C10$delabel[which(rownames(edger_C10) %in% genelist)] <- rownames(edger_C10)[which(rownames(edger_C10) %in% genelist)]

volcano_C10 = ggplot(data=edger_C10, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(box.padding = 6,size = 5,color="black") +
  scale_color_manual(values=c("plum1", "grey","orchid1")) +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  xlim(-15,15) +ylim(0,100) +  ggtitle("Mouse-derived day 0 Vs. day 2 Females") +
  geom_text(size=6, x=-10, y=90, label=paste(table(edger_C10$diffexpressed)[1]), color="black") +
  geom_text(size=6,x=10, y=90, label=paste(table(edger_C10$diffexpressed)[3]), color="black") +
  theme(plot.title = element_text(size=18),legend.position = "none",axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16), axis.title.x = element_text(size=16),axis.text.x = element_text(size=16))


## C11
# add a column of NAs
edger_C11$diffexpressed <- "NS"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
edger_C11$diffexpressed[edger_C11$logFC > 1 & edger_C11$PValue < 0.05] <- "Up"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
edger_C11$diffexpressed[edger_C11$logFC < -1 & edger_C11$PValue < 0.05] <- "Down"

table(edger_C11$diffexpressed)

edger_C11$delabel <- NA
edger_C11$delabel[which(rownames(edger_C11) %in% genelist)] <- rownames(edger_C11)[which(rownames(edger_C11) %in% genelist)]

volcano_C11 = ggplot(data=edger_C11, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(box.padding = 6,size = 5,color="black") +
  scale_color_manual(values=c("plum1", "grey","orchid1")) +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  xlim(-15,15) +ylim(0,100) +  ggtitle("Mouse-derived day 0 Vs. day 2 Males") +
  geom_text(size=6, x=-10, y=90, label=paste(table(edger_C11$diffexpressed)[1]), color="black") +
  geom_text(size=6,x=10, y=90, label=paste(table(edger_C11$diffexpressed)[3]), color="black") +
  theme(plot.title = element_text(size=18),legend.position = "none",axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16), axis.title.x = element_text(size=16),axis.text.x = element_text(size=16))

## C12
# add a column of NAs
edger_C12$diffexpressed <- "NS"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
edger_C12$diffexpressed[edger_C12$logFC > 1 & edger_C12$PValue < 0.05] <- "Up"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
edger_C12$diffexpressed[edger_C12$logFC < -1 & edger_C12$PValue < 0.05] <- "Down"

table(edger_C12$diffexpressed)

edger_C12$delabel <- NA
edger_C12$delabel[which(rownames(edger_C12) %in% genelist)] <- rownames(edger_C12)[which(rownames(edger_C12) %in% genelist)]

volcano_C12 = ggplot(data=edger_C12, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(box.padding = 6,size = 5,color="black") +
  scale_color_manual(values=c("plum1", "grey","orchid1")) +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  xlim(-15,15) +ylim(0,100) +  ggtitle("Mouse-derived day 0 Vs. day 3 Females") +
  geom_text(size=6, x=-10, y=90, label=paste(table(edger_C12$diffexpressed)[1]), color="black") +
  geom_text(size=6,x=10, y=90, label=paste(table(edger_C12$diffexpressed)[3]), color="black") +
  theme(plot.title = element_text(size=18),legend.position = "none",axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16), axis.title.x = element_text(size=16),axis.text.x = element_text(size=16))


## C13
# add a column of NAs
edger_C13$diffexpressed <- "NS"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
edger_C13$diffexpressed[edger_C13$logFC > 1 & edger_C13$PValue < 0.05] <- "Up"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
edger_C13$diffexpressed[edger_C13$logFC < -1 & edger_C13$PValue < 0.05] <- "Down"

table(edger_C13$diffexpressed)

edger_C13$delabel <- NA
edger_C13$delabel[which(rownames(edger_C13) %in% genelist)] <- rownames(edger_C13)[which(rownames(edger_C13) %in% genelist)]

volcano_C13 = ggplot(data=edger_C13, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(box.padding = 6,size = 5,color="black") +
  scale_color_manual(values=c("plum1", "grey","orchid1")) +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  xlim(-15,15) +ylim(0,100) +  ggtitle("Mouse-derived day 0 Vs. day 3 Males") +
  geom_text(size=6, x=-10, y=90, label=paste(table(edger_C13$diffexpressed)[1]), color="black") +
  geom_text(size=6,x=10, y=90, label=paste(table(edger_C13$diffexpressed)[3]), color="black") +
  theme(plot.title = element_text(size=18),legend.position = "none",axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16), axis.title.x = element_text(size=16),axis.text.x = element_text(size=16))

## C14
# add a column of NAs
edger_C14$diffexpressed <- "NS"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
edger_C14$diffexpressed[edger_C14$logFC > 1 & edger_C14$PValue < 0.05] <- "Up"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
edger_C14$diffexpressed[edger_C14$logFC < -1 & edger_C14$PValue < 0.05] <- "Down"

table(edger_C14$diffexpressed)

edger_C14$delabel <- NA
edger_C14$delabel[which(rownames(edger_C14) %in% genelist)] <- rownames(edger_C14)[which(rownames(edger_C14) %in% genelist)]

volcano_C14 = ggplot(data=edger_C14, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(box.padding = 6,size = 5,color="black") +
  scale_color_manual(values=c("plum1", "grey","orchid1")) +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  xlim(-15,15) +ylim(0,100) +  ggtitle("Mouse-derived day 0 Vs. day 5 Females") +
  geom_text(size=6, x=-10, y=90, label=paste(table(edger_C14$diffexpressed)[1]), color="black") +
  geom_text(size=6,x=10, y=90, label=paste(table(edger_C14$diffexpressed)[3]), color="black") +
  theme(plot.title = element_text(size=18),legend.position = "none",axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16), axis.title.x = element_text(size=16),axis.text.x = element_text(size=16))


## C15
# add a column of NAs
edger_C15$diffexpressed <- "NS"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
edger_C15$diffexpressed[edger_C15$logFC > 1 & edger_C15$PValue < 0.05] <- "Up"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
edger_C15$diffexpressed[edger_C15$logFC < -1 & edger_C15$PValue < 0.05] <- "Down"

table(edger_C15$diffexpressed)

edger_C15$delabel <- NA
edger_C15$delabel[which(rownames(edger_C15) %in% genelist)] <- rownames(edger_C15)[which(rownames(edger_C15) %in% genelist)]

volcano_C15 = ggplot(data=edger_C15, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(box.padding = 6,size = 5,color="black") +
  scale_color_manual(values=c("plum1", "grey","orchid1")) +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  xlim(-15,15) +ylim(0,100) +  ggtitle("Mouse-derived day 0 Vs. day 5 Males") +
  geom_text(size=6, x=-10, y=90, label=paste(table(edger_C15$diffexpressed)[1]), color="black") +
  geom_text(size=6,x=10, y=90, label=paste(table(edger_C15$diffexpressed)[3]), color="black") +
  theme(plot.title = element_text(size=18),legend.position = "none",axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16), axis.title.x = element_text(size=16),axis.text.x = element_text(size=16))

## C16
# add a column of NAs
edger_C16$diffexpressed <- "NS"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
edger_C16$diffexpressed[edger_C16$logFC > 1 & edger_C16$PValue < 0.05] <- "Up"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
edger_C16$diffexpressed[edger_C16$logFC < -1 & edger_C16$PValue < 0.05] <- "Down"

table(edger_C16$diffexpressed)

edger_C16$delabel <- NA
edger_C16$delabel[which(rownames(edger_C16) %in% genelist)] <- rownames(edger_C16)[which(rownames(edger_C16) %in% genelist)]

volcano_C16 = ggplot(data=edger_C16, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(box.padding = 6,size = 5,color="black") +
  scale_color_manual(values=c("plum1", "grey","orchid1")) +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  xlim(-15,15) +ylim(0,100) +  ggtitle("Mouse-derived day 0 Vs. day 7 Females") +
  geom_text(size=6, x=-10, y=90, label=paste(table(edger_C16$diffexpressed)[1]), color="black") +
  geom_text(size=6,x=10, y=90, label=paste(table(edger_C16$diffexpressed)[3]), color="black") +
  theme(plot.title = element_text(size=18),legend.position = "none",axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16), axis.title.x = element_text(size=16),axis.text.x = element_text(size=16))


## C17
# add a column of NAs
edger_C17$diffexpressed <- "NS"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
edger_C17$diffexpressed[edger_C17$logFC > 1 & edger_C17$PValue < 0.05] <- "Up"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
edger_C17$diffexpressed[edger_C17$logFC < -1 & edger_C17$PValue < 0.05] <- "Down"

table(edger_C17$diffexpressed)

edger_C17$delabel <- NA
edger_C17$delabel[which(rownames(edger_C17) %in% genelist)] <- rownames(edger_C17)[which(rownames(edger_C17) %in% genelist)]

volcano_C17 = ggplot(data=edger_C17, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(box.padding = 6,size = 5,color="black") +
  scale_color_manual(values=c("plum1", "grey","orchid1")) +
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  xlim(-15,15) +ylim(0,100) +  ggtitle("Mouse-derived day 0 Vs. day 7 Males") +
  geom_text(size=6, x=-10, y=90, label=paste(table(edger_C17$diffexpressed)[1]), color="black") +
  geom_text(size=6,x=10, y=90, label=paste(table(edger_C17$diffexpressed)[3]), color="black") +
  theme(plot.title = element_text(size=18),legend.position = "none",axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16), axis.title.x = element_text(size=16),axis.text.x = element_text(size=16))



pdf("./DE/edgeR_volcano_MouseVTvsday0_males_lfc1_pval0.05_labels.pdf",width = 18, height =8)
ggarrange (volcano_C9,volcano_C11, volcano_C13, volcano_C15, volcano_C17) 

#WTKO Comparisons volcano_C1,volcano_C2,volcano_C3,volcano_C4hr,volcano_C5

# VT comparisons 
#day 0 volcano_C6,volcano_C7
#day 1 volcano_C8, volcano_C9 
#day 2 volcano_C10, volcano_C11
#day 3 volcano_C12, volcano_C13 
#day 5 volcano_C14, volcano_C15
#day 7 volcano_C16, volcano_C17
#day1-7 females
#volcano_C8,volcano_C10, volcano_C12, volcano_C14, volcano_C16

#day1-7 males
#volcano_C9,volcano_C11, volcano_C13, volcano_C15, volcano_C17
dev.off()

################### write out DEGs table ################

##WTvKO

write.table (edger_C1, "DE/DEGs_Lst/Nippo_WTvKO_fem_5d_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)
write.table (edger_C2, "DE/DEGs_Lst/Nippo_WTvKO_fem_7d_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)

write.table (edger_C3, "DE/DEGs_Lst/Nippo_WTvKO_male_5d_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)
write.table (edger_C4, "DE/DEGs_Lst/Nippo_WTvKO_male_7d_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)
write.table (edger_C5, "DE/DEGs_Lst/Nippo_WTvKO_L4_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)

##VitroMvR Timcourse

write.table (edger_C6, "DE/DEGs_Lst/Nippo_VT_MvR_Fem_d0_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)
write.table (edger_C7, "DE/DEGs_Lst/Nippo_VT_MvR_Male_d0_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)

write.table (edger_C8, "DE/DEGs_Lst/Nippo_VT_MvR_Fem_d1_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)
write.table (edger_C9, "DE/DEGs_Lst/Nippo_VT_MvR_Male_d1_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)

write.table (edger_C10, "DE/DEGs_Lst/Nippo_VT_MvR_Fem_d2_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)
write.table (edger_C11, "DE/DEGs_Lst/Nippo_VT_MvR_Male_d2_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)

write.table (edger_C12, "DE/DEGs_Lst/Nippo_VT_MvR_Fem_d3_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)
write.table (edger_C13, "DE/DEGs_Lst/Nippo_VT_MvR_Male_d3_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)

write.table (edger_C14, "DE/DEGs_Lst/Nippo_VT_MvR_Fem_d5_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)
write.table (edger_C15, "DE/DEGs_Lst/Nippo_VT_MvR_Male_d4_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)

write.table (edger_C16, "DE/DEGs_Lst/Nippo_VT_MvR_Fem_d7_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)
write.table (edger_C17, "DE/DEGs_Lst/Nippo_VT_MvR_Male_d7_DEOutput_matrix_org_sex_tpmabv1_.txt", row.names = TRUE, col.names = T, sep = '\t', quote = F)

