# Heatmap 
        
library(ggplot2)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(viridis)
library(data.table)
library(magick)

#setwd("C:/Users/prdrm/Desktop/Research_Mortazavi/Adler_Nippostrongylus_Project/")

####### first upload the TPM file ##########
tpm_0 <- read.delim('data/filtered/Nippo_Allsamples_02-13-23_tpm_matrix_col_filtabv1.txt',header=TRUE,row.names="gene_id")
head(tpm_1)

meta = read.delim('Metadata_samples_2-20-23.txt',header=TRUE, stringsAsFactors = FALSE)
head(meta)


meta$sample_name = gsub(" ",".",meta$sample_name)
meta$sample_name = gsub("-",".",meta$sample_name)

# Check that metadata matches counts
table(colnames(tpm_0) == meta$sample_name)


head(tpm_0)
dim(tpm_0)

# sanity check that all samples in metadata are in colnames of data
table(meta$sample_name %in% colnames(tpm_0))

############ Vitro sets day 0 ONLY ("VT_0d.*", meta$stage8_sampletype)
########### vitro timecourse ("VT.*d", meta$stage8_sampletype)


tpm_1=tpm_0[,grep ("VT.*d", meta$stage8_sampletype)] ## subset the data matrix 

head (tpm_1) 

#subset meta

meta_1=meta[match(colnames(tpm_1), meta$sample_name),]  ## subset the metadata to match the data
head(meta_1)

# sanity check that subsetted matrix and subsetted metadata match
table(meta_1$sample_name == colnames(tpm_1)) 


#################### process###############################
# Remove rows where max TPM is less than 0 - 936 left
tpm_1 = tpm_1 [rowSums(tpm_1[,-1])>0,]
dim(tpm_1) # 21334(filt abv 1 TPM) by 213 samples



#################### Heatmap and clustering #################### 
tpm_0_log2_scaled=as.data.frame(t(apply(log2(tpm_1+1), MARGIN = 1, FUN = function(x){(x-min(x))/(max(x)-min(x))})))

tpm_0_log2_scaled=as.matrix(tpm_0_log2_scaled)


######## Heatmap of all genes #########
# Randomly sample rows, using full matrix crashes bc it's too big
# Keep samples in order
#pdf(file = "Heatmap_log2tpm_KOWT_rowScaled.pdf", width =10, height = 10)
Heatmap(name = "log2(TPM+1)", as.matrix(tpm_0_log2_scaled), # sample 5,000 genes 
        column_labels = as.character(colnames(tpm_0_log2_scaled)), col=viridis(100),
        show_row_names = FALSE, show_row_dend = FALSE, cluster_columns = FALSE, # don't cluster samples
        show_column_names=TRUE, column_title = "Sample", column_names_gp = gpar(fontsize = 12),
        row_title = "genes", use_raster=T)

#dev.off()


##### Add just genes of intrest ####


endocanna = c("NBR_0001270801","NBR_0001219601")

blood_dig = c("NBR_0000892601","NBR_0001095601")


ESP = c()

genes_1= c(endocanna,blood_dig)



tpm_0_log2_scaled_subset = tpm_0_log2_scaled[rownames(tpm_0_log2_scaled) %in% genes_1,]
tpm_0_log2_scaled_subset = tpm_0_log2_scaled_subset[match(genes, rownames(tpm_0_log2_scaled_subset)),]


#pdf(file = "Heatmap_log2tpm_rowScaled_genes_ofinterest.pdf", width =10, height = 10)
Heatmap(name = "log2(TPM+1)", as.matrix(tpm_0_log2_scaled_subset), 
        column_labels = as.character(colnames(tpm_0_log2_scaled_subset)), col=viridis(100),
        show_row_names = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, # don't cluster samples
        show_column_names=TRUE, column_title = "Sample", column_names_gp = gpar(fontsize = 12),
        row_title = "genes", use_raster=T)
#dev.off()

