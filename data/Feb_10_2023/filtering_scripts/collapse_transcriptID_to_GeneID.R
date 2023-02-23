## Read the gene and transcript correspondence information
gene.info <- read.table("gene_id_trans_id_Nippo.txt", header = T, stringsAsFactors = F)
head (gene.info)

## Read your text file with trans and gene id file
counts <- read.table("Nippo_Allsamples_02-13-23_efflength_matrix_rm_nm.txt", header = T, stringsAsFactors = F)
head(counts)

## Replace the transcript IDs in TPM file with the gene IDs
counts$gene_id <- gene.info$gene_id[match(counts$gene_id, gene.info$trans_id)]
head(counts)

## Collapse the rows with the same gene ID
require(plyr)
collaped.counts <- ddply(counts, "gene_id", numcolwise(sum))
head(collaped.counts)
write.table(collaped.counts, "Nippo_Allsamples_02-13-23_efflength_matrix_rm_nm_col.txt", row.names = F, col.names = T,
            quote = F, sep = '\t')

