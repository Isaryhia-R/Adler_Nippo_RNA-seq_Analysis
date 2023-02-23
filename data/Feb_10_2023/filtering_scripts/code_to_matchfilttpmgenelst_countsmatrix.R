filtgenes <- read.table("genes_filtabv1TPM_col.txt", header = T, stringsAsFactors = FALSE)
head(filtgenes)
counts <- read.table("Nippo_Allsamples_02-13-23_efflength_matrix_rm_nm_col.txt", header = T, stringsAsFactors = FALSE)
head(counts)
filtgene.counts  <- filtgenes[match(counts$gene_id, filtgenes$gene_id), ]
head(filtgene.counts)
combined.data <- cbind(filtgene.counts, counts)# adding the V2 here allows it to grab the corresponding gene in the row with the genes isolated
head (combined.data)
finaldata <- na.omit(combined.data)
write.table (finaldata, "Nippo_Allsamples_02-13-23_efflength_matrix_rm_col_filabv1TPM.txt", row.names = F, col.names = T, sep = '\t', quote = F)

