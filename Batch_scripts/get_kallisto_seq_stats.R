setwd("/central/groups/stathopouloslab/imrodrig/Adler_Nippo/Nippostronglyus_project/scripts/")

library(stringr)

files = list.files("/central/groups/stathopouloslab/imrodrig/Adler_Nippo/Nippostronglyus_project/scripts/kallisto_stats/", pattern=".err$", all.files=TRUE, full.names=FALSE) 

output = as.data.frame(matrix( nrow = length(files), ncol = 4))
colnames(output) = c("sample","total_reads","aligned_reads","percent_aligned")

err = read.delim(files[1], stringsAsFactors = F)

total_reads = gsub("[^0-9.-]", "", sapply(strsplit(as.character(err[9,]), " reads,"), "[[", 1))
aligned_reads = gsub("[^0-9.-]", "", sapply(strsplit(as.character(err[9,]), " reads,"), "[[", 2))
fasta_name = sapply(strsplit(as.character(err[6,]), "fastq"), "[[", 2)
percent_aligned = as.numeric(aligned_reads)/as.numeric(total_reads)

output[1,] = c(fasta_name,total_reads,aligned_reads,percent_aligned)

for (i in 2:length(files)) {
  err = read.delim(files[i], stringsAsFactors = F)
  
  total_reads = gsub("[^0-9.-]", "", sapply(strsplit(as.character(err[9,]), " reads,"), "[[", 1))
  aligned_reads = gsub("[^0-9.-]", "", sapply(strsplit(as.character(err[9,]), " reads,"), "[[", 2))
  fasta_name = sapply(strsplit(as.character(err[6,]), "fastq"), "[[", 2)
  percent_aligned = as.numeric(aligned_reads)/as.numeric(total_reads)
  output[i,] = c(fasta_name,total_reads,aligned_reads,percent_aligned)
}


write.table(output,file="./sequencing_stats.tsv",sep="\t", row.names = F, quote = F)
