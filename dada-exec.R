library(dada2)
path <- "/mnt/research/germs/soil-column2-16S"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[1], x[2], sep="_"))

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
filtFs2 <- filtFs[out[,2] >= 1000]
filtRs2 <- filtRs[out[,2] >= 1000]
errF <- learnErrors(filtFs2, multithread=TRUE)
errR <- learnErrors(filtRs2, multithread=TRUE)
dadaFs <- dada(filtFs2, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs2, err=errR, multithread=TRUE)
save.image("/mnt/research/germs/soil-column2-16S/dada.RData")

mergers <- mergePairs(dadaFs, filtFs2, dadaRs, filtRs2, verbose=TRUE)



head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
save.image("/mnt/research/germs/soil-column2-16S/dada2.RData")