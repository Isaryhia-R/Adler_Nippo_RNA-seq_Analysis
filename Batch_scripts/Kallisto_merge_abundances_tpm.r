#################### Load metadata #################### 
metadata = read.csv('/central/groups/stathopouloslab/imrodrig/Adler_Nippo/Nippostronglyus_project/scripts/Nippo_metadata.csv', head=TRUE, sep=',')
samples=metadata$Sample
#################### Read in .tsv and make counts matrix across all samples in data folder ###############################
thisMat=as.data.frame(read.table(paste0("/central/groups/stathopouloslab/imrodrig/Adler_Nippo/Nippostronglyus_project/output/2023_02_13/",samples[1],"/abundance.tsv"),sep="",head=T, stringsAsFactors = F))
head(thisMat)
counts=thisMat[,c('target_id','tpm')] # est_counts or tpm or eff_length
head(counts)

for(i in 2:length(samples)){
  thisMat=as.data.frame(read.table(paste0("/central/groups/stathopouloslab/imrodrig/Adler_Nippo/Nippostronglyus_project/output/2023_02_13/",samples[i],"/abundance.tsv"),sep="",head=T, stringsAsFactors = F))
  thisMat=thisMat[na.omit(match(counts$target_id,thisMat$target_id)),]
  counts=as.data.frame(cbind(counts,thisMat[,"tpm"]))   #est_counts  or tpm or eff_length (CHANGE NAME)
}

dim(counts)
rownames(counts) = counts$target_id
counts$target_id = NULL
colnames(counts) = samples
head(counts)

write.csv(counts,"Nippo_Allsamples_2-13-23_tpm_matrix.csv", quote = FALSE)  ##CHANGE NAME IF TPM or counts!! 
~                                
