samples <- read.table("samplesheet_jobsub_34-62.txt", header = T)
print(samples)
# Step01_runFastqc.sub Step02_runkallisto.sub
numOfRows <- dim(samples)[1]
for (i in 1:numOfRows){
  currentRow <- unlist(samples[i, ])
  system(paste0("sbatch Step01_runFastqc.sub " , currentRow[1], sep=""))}
