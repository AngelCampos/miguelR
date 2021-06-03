library(txtools)
library(pasillaBamSubset)
library(GenomicRanges)
library(magrittr)
library(stringr)
library(Biostrings)
library(data.table)

# BED_file <- tx_dm3_geneAnnot()
# FASTA_file <- dm3_chr4()
# PE_BAM_file <- untreated3_chr4()
# 
# dm3_geneAnnot <- tx_load_bed(BED_file)
# dm3_genome <- tx_load_genome(FASTA_file)
# dm3_PEreads <- tx_load_bam(file = PE_BAM_file, pairedEnd = T, loadSeq = T)
# 
# 
# reads_SE <- tx_reads(reads = dm3_PEreads, 
#                      geneAnnot = dm3_geneAnnot, 
#                      withSeq = TRUE, 
#                      nCores = 1, 
#                      minReads = 1)
# seqnames(dm3_PEreads)
# start(dm3_PEreads@last)



# Artifical genome
mk_genome <- function(gLen, nucs = c("A", "T", "G", "C"), prob = c(0.25, 0.25, 0.25, 0.25)){
    genome <- lapply(gLen, function(x){
        sample(nucs, size = x, replace = TRUE, prob = prob) %>% 
            paste(collapse = "") 
    }) %>% unlist %>% DNAStringSet() %>% set_names(paste("chr", seq_along(gLen), sep = ""))
    
}


# Create artifical gene annotation
genome <- mk_genome(1000)


mk_geneModels <- function(genome, nGenes, nExons = "auto", lambda_exons = 30, lambda_introns = 10){
    if(nExons[1] == "auto"){nExons = rpois(n = nGenes, lambda = 1) + 1}
    nIntrons <- nExons - 1
    l_e <- lapply(nExons, function(x) rpois(n = x, lambda = lambda_exons))
    l_i <- lapply(nIntrons, function(x) c(rpois(n = x, lambda = lambda_introns), 0))
    len_b <- lapply(seq_along(l_e), function(i){
        idx <- order(c(seq_along(l_e[[i]]), seq_along(l_i[[i]])))
        unlist(c(l_e[[i]],l_i[[i]]))[idx]
    })
    totalLen <- lapply(len_b, sum) %>% unlist()
    if(any(totalLen >= width(genome))){stop("Some genes exceeded genome ", 
    "size.\nIncrease genome size or modify lambda parameters, or try your luck again;)")}
    f_tbl <- data.frame(chrom = rep(names(genome), nGenes))
    f_tbl$start <- lapply(seq_along(totalLen), function(i){
        sample(x = 1:(nchar(genome) - totalLen[i] - 1), size = 1)
    }) %>% unlist()
    f_tbl$end <- f_tbl$start + totalLen
    f_tbl$name <- paste("gene", seq_along(totalLen), sep = "_")
    f_tbl$score <- 0
    f_tbl$strand <- sample(x = c("+", "-"), size = nGenes, replace = TRUE)
    f_tbl$thickStart <- f_tbl$start
    f_tbl$thickEnd <- f_tbl$end
    f_tbl$itemRGB <- "0,0,0"
    f_tbl$blockCount <- nExons
    f_tbl$blockSize <- lapply(l_e, function(x) paste0(x, ",", collapse = "")) %>% unlist()
    f_tbl$blockStarts <- lapply(seq_along(l_e), function(i) c(0, l_e[[i]] + l_i[[i]])) %>% 
        lapply(function(x) head(x, -1)) %>% 
        lapply(function(x) paste0(x, ",", collapse = ""))
    tmpF <- tempfile()
    fwrite(x = f_tbl, file = tmpF, col.names = FALSE, sep = "\t")
    tx_load_bed(tmpF)
}

geneAnnot <- mk_geneModels(genome, nGenes = 5, lambda_exons = 100)

# Artificial FASTQ data

get_tx_seq <- function(geneAnnot, genome){
    iGenes <- as.character(geneAnnot$name)
    iChrs <- as.character(seqnames(geneAnnot))
    iStrs <- as.character(strand(geneAnnot))
    iBlocks <- lapply(seq_along(geneAnnot), function(i){
        S4Vectors::mcols(geneAnnot[i])$blocks %>% txtools:::if_IRangesList_Unlist() %>%
            IRanges::shift(IRanges::start(geneAnnot[i]) - 1)
    })
    seqs <- lapply(seq_along(iGenes), function(i){
        stringr::str_sub(genome[[iChrs[i]]], start = IRanges::start(iBlocks[[i]]),
                         end = IRanges::end(iBlocks[[i]])) %>% paste(collapse = "") %>%
            Biostrings::DNAString()
    }) %>% DNAStringSet() 
    lapply(seq_along(iGenes), function(i){
        if (iStrs[i] == "-") {
            Biostrings::reverseComplement(seqs[[i]])
        }else{
            seqs[[i]]
        }
    }) %>% DNAStringSet()
}

txSeqs <- get_tx_seq(geneAnnot, genome)
writeLines(as.character(txSeqs), con = "toySeqs.txt")
writeLines(as.character(genome), con = "toyGenome.txt")
pairedEnd <- TRUE
R1_len <- 30
R2_len <- 30
lambda_reads <- 100
0

mk_fastqReads <- function(txSeqs, pairedEnd = TRUE, lambda_reads){
    nReads <- rpois(length(txSeqs), lambda = lambda_reads)
    # rnbinom(length(txSeqs), size = 100, prob = 0.5)
    lapply(nExons, function(x) rpois(n = x, lambda = lambda_exons))
    R1_len
    R2_len
}

genPEread


library(insect)
?insect::writeFASTQ()
# Map using subread

# Retrieve using txtools

## download and extract example FASTQ file to temporary directory
td <- tempdir()
URL <- "https://www.dropbox.com/s/71ixehy8e51etdd/insect_tutorial1_files.zip?dl=1"
dest <- paste0(td, "/insect_tutorial1_files.zip")
download.file(URL, destfile = dest, mode = "wb")
unzip(dest, exdir = td)
insect::char2dna("CGTCAGTAGTCA")
 
x <- readFASTQ(paste0(td, "/COI_sample2.fastq"))
## trim primers from sequences
mlCOIintF <- "GGWACWGGWTGAACWGTWTAYCCYCC"
jgHCO2198 <- "TAIACYTCIGGRTGICCRAARAAYCA"
x <- trim(x, up = mlCOIintF, down = jgHCO2198)
## quality filter with size selection and singleton removal
x <- qfilter(x, minlength = 250, maxlength = 350)
## output filtered FASTQ file
writeFASTQ(x, file = paste0(td, "/COI_sample2_filtered.fastq"))
writeFASTA(x, file = paste0(td, "/COI_sample2_filtered.fasta"))



## Toy objects for check
set.seed(100)
toyGenome <- mk_genome(1000)
toyGAnnot <- mk_geneModels(toyGenome, 3, nExons = c(1,2,3), lambda_introns = 50, lambda_exons = 100)
toyTxSeqs <- get_tx_seq(genome = toyGenome, geneAnnot = toyGAnnot)
writeLines(as.character(toyGenome), con = "toyGenome.txt")
writeLines(as.character(toyGenome), con = "toyGenome.txt")
writeLines(as.character(toyTxSeqs), con = "toySeqs.txt")
writeLines(as.character(Biostrings::reverseComplement(toyTxSeqs)), con = "toySeqs_revComp.txt")
writeLines(as.character(Biostrings::reverseComplement(toyGenome)), con = "toyGeno_revComp.txt")
