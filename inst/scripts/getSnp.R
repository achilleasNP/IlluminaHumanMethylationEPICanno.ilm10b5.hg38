processUCSCsnp <- function(snpfile) {
    require(GenomicRanges)
    cat("Reading file\n")
    df <- read.delim(gzfile(snpfile), header = FALSE,
                     stringsAsFactors = FALSE)
    names(df) <- c("chr", "start", "end", "name", "strand",
                   "refNCBI", "class", "alleleFreqs")
    print(table(df$chr))
    cat("Only keeping chrs 1-22, X, Y\n")
    df <- df[df$chr %in% paste0("chr", c(1:22, "X", "Y")),]
    print(table(df$class))
    cat("Only keeping class 'single'\n")
    df <- df[df$class == "single",]
    cat("Computing MAF\n")
    df$alleleFreqs <- sub(",$", "", df$alleleFreqs)
    sp <- strsplit(df$alleleFreqs, ",")
    minFreq <- sapply(sp, function(xx) min(as.numeric(xx)))
    cat("Instantiating object\n")
    grSNP <- GRanges(seqnames = df$chr, strand = df$strand,
                     ranges = IRanges(start = df$start + 1, end = df$end),
                     MAF = minFreq, ref = df$refNCBI)
    names(grSNP) <- df$name
    grSNP
}


grSnp132CommonSingle <- processUCSCsnp("files/snp132Common_small.txt.gz")
save(grSnp132CommonSingle, file = "objects/grSnp132CommonSingle.rda")

grSnp135CommonSingle <- processUCSCsnp("files/snp135Common_small.txt.gz")
save(grSnp135CommonSingle, file = "objects/grSnp135CommonSingle.rda")

grSnp137CommonSingle <- processUCSCsnp("files/snp137Common_small.txt.gz")
save(grSnp137CommonSingle, file = "objects/grSnp137CommonSingle.rda")

grSnp138CommonSingle <- processUCSCsnp("files/snp138Common_small.txt.gz")
save(grSnp138CommonSingle, file = "objects/grSnp138CommonSingle.rda")

grSnp141CommonSingle <- processUCSCsnp("files/snp141Common_small.txt.gz")
save(grSnp141CommonSingle, file = "objects/grSnp141CommonSingle.rda")

grSnp142CommonSingle <- processUCSCsnp("files/snp142Common_small.txt.gz")
save(grSnp142CommonSingle, file = "objects/grSnp142CommonSingle.rda")

grSnp144CommonSingle <- processUCSCsnp("files/snp144Common_small.txt.gz")
save(grSnp144CommonSingle, file = "objects/grSnp144CommonSingle.rda")

grSnp146CommonSingle <- processUCSCsnp("files/snp146Common_small.txt.gz")
save(grSnp146CommonSingle, file = "objects/grSnp146CommonSingle.rda")

grSnp147CommonSingle <- processUCSCsnp("files/snp147Common_small.txt.gz")
save(grSnp147CommonSingle, file = "objects/grSnp147CommonSingle.rda")


