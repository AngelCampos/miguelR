library(txtools)
library(pasillaBamSubset)
library(GenomicRanges)

BED_file <- tx_dm3_geneAnnot()
FASTA_file <- dm3_chr4()
PE_BAM_file <- untreated3_chr4()

dm3_geneAnnot <- tx_load_bed(BED_file)
dm3_genome <- tx_load_genome(FASTA_file)
dm3_PEreads <- tx_load_bam(file = PE_BAM_file, pairedEnd = T, loadSeq = T)


reads_SE <- tx_reads(reads = dm3_PEreads, 
                     geneAnnot = dm3_geneAnnot, 
                     withSeq = TRUE, 
                     nCores = 1, 
                     minReads = 1)
seqnames(dm3_PEreads)
start(dm3_PEreads@last)

library(magrittr)
library(stringr)

rep("chr")

# Artifical genome
mk_genome <- function(gLen, nucs = c("A", "T", "G", "C"), prob = c(0.25, 0.25, 0.25, 0.25)){
    genome <- lapply(gLen, function(x){
        sample(nucs, size = x, replace = TRUE, prob = prob) %>% 
            paste(collapse = "") 
    }) %>% unlist %>% DNAStringSet()
}

testFun <- Vectorize(mk_genome, vectorize.args = "gLen")
# Create artifical gene annotation
genome <- mk_transcriptome(c(300, 200))

nGenes <- 3 
nExons <- c(2, 1, 4)
mkGeneModels <- function(genome, nGenes, exons){
    
}

sample(x = 1:100, 2) %>% sort

rpois(10000, 100) %>% min
split(df, sample(1:N, nrow(df), replace=T))

df <- data.frame(start=c(2:-1, 13:15),
                 width=c(0:3, 2:0))

# produces IRanges
rng <- df %>% as_iranges()
rng

df %>% transform(seqnames = sample(c("chr1", "chr2"), 7, replace = TRUE),
              strand = sample(c("+", "-"), 7, replace = TRUE),
              gc = runif(7)) %>%
    as_granges()

library(plyranges)
plyranges::as_granges()
# Artificial FASTQ data

# Map using subread

# Retrieve using txtools

